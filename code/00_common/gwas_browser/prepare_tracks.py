#!/usr/bin/env python3
"""Prepare GWAS, line and gene tracks for the local GWAS browser.

Reads a track manifest (see README.md, "Source manifest") and writes the binary/JSON layout described in
README.md, "Prepared data contract", to output_data/gwas_browser/ by default. Tracks whose sources and options
are unchanged are skipped.

    prepare_tracks.py MANIFEST [--out DIR] [--only id,id] [--force] [--jobs N]
    prepare_tracks.py --discover > tracks.json
    prepare_tracks.py --selftest
"""

import argparse
import array
import datetime
import glob
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from urllib.parse import unquote

import numpy as np

PREP_VERSION = 1          # bump to force a rebuild of every track after a format/logic change
LOGP_SCALE = 1000
U16_MAX = 65535
BIN_SIZE = 10000
DETAIL_CUTOFF = 2.0
OVERVIEW_CUTOFF = 2.0
OVERVIEW_MAX = 300000
# keep candidate detail rows a little above p = 0.01 so the scaled cutoff (>= 2000) decides
DETAIL_P_MAX = 10 ** -(DETAIL_CUTOFF - 0.001)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_FAI = ("input_data/00_Reference_Genome/"
               "jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai")
DEFAULT_GFF = "input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.sorted.gff"
TRACK_ID_RE = re.compile(r"^[A-Za-z0-9_-]+$")

P_COLS = ["p_wald", "p_lrt", "p_score", "p"]
TEST_BY_PCOL = {"p_wald": "lmm1", "p_lrt": "lmm2", "p_score": "lmm3", "p": "gwas"}
MIN_COL_BY_TEST = {"lmm1": 12, "lmm2": 10}   # p column in min.txt (05_collect_min_values.sh)
DETAIL_SOURCES = [            # (contract field, accepted source header names, type)
    ("rs", ["rs", "snp"], str),
    ("p", None, str),
    ("af", ["af"], float),
    ("beta", ["beta"], float),
    ("se", ["se"], float),
    ("n_miss", ["n_miss"], int),
    ("allele1", ["allele1"], str),
    ("allele0", ["allele0"], str),
]

LINE_PRESETS = {
    "pixy_pi":  {"chr": "chromosome", "start": "window_pos_1", "end": "window_pos_2", "value": "avg_pi",
                 "filter": "pop", "kind": "window"},
    "pixy_fst": {"chr": "chromosome", "start": "window_pos_1", "end": "window_pos_2", "value": "avg_wc_fst",
                 "filter": "pair", "kind": "window", "ymin": 0, "ymax": 1},
    "pixy_dxy": {"chr": "chromosome", "start": "window_pos_1", "end": "window_pos_2", "value": "avg_dxy",
                 "filter": "pair", "kind": "window"},
    "tajd":     {"chr": "CHROM", "start": "BIN_START", "end": None, "value": "TajimaD",
                 "filter": "pop1", "kind": "window", "zero_line": True},
    "dsuite":   {"chr": "chr", "start": "windowStart", "end": "windowEnd", "value": "f_dM",
                 "filter": None, "kind": "window", "ymin": -1, "ymax": 1, "zero_line": True},
    "sweed":    {"chr": "chromosome", "start": "position", "end": None, "value": "CLR",
                 "filter": None, "kind": "point"},
    "sf2":      {"chr": "chromosome", "start": "position", "end": None, "value": "CLR",
                 "filter": None, "kind": "point"},
}


# ----------------------------------------------------------------------------------------------- helpers

def log(msg):
    print(msg, flush=True)


def repo_root():
    try:
        out = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=SCRIPT_DIR,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", ".."))


def rel(path, root):
    try:
        return os.path.relpath(path, root)
    except ValueError:
        return path


def open_text(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb", buffering=1 << 22)


def file_sig(path):
    st = os.stat(path)
    return [rel(path, ROOT), st.st_size, int(st.st_mtime)]


def opt_sig(path):
    return file_sig(path) if path and os.path.exists(path) else None


def options_hash(obj):
    return hashlib.sha1(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()


def atomic_write_bytes(path, data):
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, path)


def atomic_write_json(path, obj, indent=None):
    atomic_write_bytes(path, json.dumps(obj, indent=indent, separators=None if indent else (",", ":"))
                       .encode())


def read_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def swap_dir(tmp_dir, final_dir):
    """Replace final_dir with tmp_dir (removes stale files from an earlier build)."""
    old = final_dir + ".old"
    if os.path.exists(old):
        shutil.rmtree(old)
    if os.path.exists(final_dir):
        os.replace(final_dir, old)
    os.replace(tmp_dir, final_dir)
    if os.path.exists(old):
        shutil.rmtree(old)


def contig_sort_key(name):
    m = re.fullmatch(r"chr(\d+)", name)
    if m:
        return (0, int(m.group(1)), name)
    m = re.fullmatch(r"ups(\d+)", name)
    if m:
        return (1, int(m.group(1)), name)
    return (2, 0, name)


def read_fai(path):
    lengths = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def bed_contig(name):
    """Peak BEDs from gemma.peak.bed() use numeric chr; ups contigs are dropped before peak calling
    (gwas_report.Rmd:48), so numeric N always means chrN."""
    name = name.strip()
    if re.fullmatch(r"\d+", name):
        return "chr" + str(int(name))
    return name


def thresholds_from_pvals(pvals):
    """gwas_report.Rmd:35-39. R: pvals[length - floor(length*0.95)] with 1-based indexing."""
    p = sorted(pvals)
    n = len(p)
    if n == 0:
        return None, None
    i_sig = n - math.floor(n * 0.95)
    i_sug = n - math.floor(n * 0.67)
    if i_sig < 1 or i_sug < 1:
        return None, None
    return -math.log10(p[i_sig - 1]), -math.log10(p[i_sug - 1])


def parse_float(s):
    try:
        return float(s)
    except ValueError:
        return math.nan


def scale_logp(p):
    """p (float64 array) -> (scaled u16 logp, n_capped). p = 0 or logp > 65.535 caps at 65535."""
    with np.errstate(divide="ignore"):
        logp = -np.log10(p)
    scaled = np.rint(logp * LOGP_SCALE)
    capped = int(np.count_nonzero(scaled > U16_MAX))
    scaled = np.clip(scaled, 0, U16_MAX).astype(np.uint16)
    return scaled, capped


# ------------------------------------------------------------------------------------------- GWAS tracks

def resolve_auto(value, path, suffix):
    """'auto' -> <dir>/<NAME><suffix>; None/'' -> None; else a root-relative path."""
    if value is None or value == "" or value is False:
        return None
    if value == "auto":
        d = os.path.dirname(path)
        name = os.path.basename(path).split(".")[0]
        cand = os.path.join(d, name + suffix)
        return cand if os.path.exists(cand) else None
    return os.path.join(ROOT, value)


def read_header(path):
    with open_text(path) as fh:
        return fh.readline().decode().rstrip("\r\n").split("\t")


def gwas_columns(header, p_col):
    low = [h.strip().lower() for h in header]
    idx = {h: i for i, h in enumerate(low)}
    chr_i = idx.get("chr")
    pos_i = idx.get("ps", idx.get("bp"))
    if chr_i is None or pos_i is None:
        raise ValueError("no chr/ps (or CHR/BP) column in header: %s" % header)
    if p_col in (None, "auto"):
        p_name = next((c for c in P_COLS if c in idx), None)
        if p_name is None:
            raise ValueError("no p column (%s) in header" % ", ".join(P_COLS))
    else:
        p_name = p_col.lower()
        if p_name not in idx:
            raise ValueError("p_col %r not in header" % p_col)
    detail = []
    for field, names, typ in DETAIL_SOURCES:
        if field == "p":
            detail.append((field, idx[p_name], typ))
            continue
        col = next((idx[n] for n in names if n in idx), None)
        if col is not None:
            detail.append((field, col, typ))
    return chr_i, pos_i, idx[p_name], header[idx[p_name]], idx.get("af"), detail


def parse_gwas(path, p_col="auto", min_af=0.0):
    """Stream a GEMMA assoc / .GWAS file. Returns per-contig position/p arrays and detail candidate rows."""
    header = read_header(path)
    chr_i, pos_i, p_i, p_name, af_i, detail_cols = gwas_columns(header, p_col)
    if min_af and af_i is None:
        log("  WARN %s: min_af set but no af column; ignored" % rel(path, ROOT))
        min_af = 0.0
    maxcol = max([chr_i, pos_i, p_i] + ([af_i] if min_af else []) + [c for _, c, _ in detail_cols])
    contigs = {}       # name -> [array('I') pos, array('d') p, list detail rows]
    order = []
    n_rows = n_nan = n_af = 0
    cur_name = None
    cur = None
    with open_text(path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip(b"\r\n").split(b"\t")
            if len(parts) <= maxcol:
                if not line.strip():
                    continue
                raise ValueError("short row in %s: %r" % (path, line[:200]))
            n_rows += 1
            try:
                p = float(parts[p_i])
            except ValueError:
                p = math.nan
            if p != p:
                n_nan += 1
                continue
            if min_af and float(parts[af_i]) < min_af:
                n_af += 1
                continue
            c = parts[chr_i]
            if c != cur_name:
                cur_name = c
                cur = contigs.get(c)
                if cur is None:
                    cur = contigs[c] = [array.array("I"), array.array("d"), []]
                    order.append(c)
            pos = int(parts[pos_i])
            cur[0].append(pos)
            cur[1].append(p)
            if p <= DETAIL_P_MAX:
                cur[2].append((pos, [parts[i] for _, i, _ in detail_cols]))
    out = {}
    for c in order:
        pos_a, p_a, det = contigs[c]
        out[c.decode()] = (np.frombuffer(pos_a, dtype=np.uint32).copy(),
                           np.frombuffer(p_a, dtype=np.float64).copy(), det)
    return {"contigs": out, "p_name": p_name, "detail_fields": [f for f, _, _ in detail_cols],
            "detail_types": [t for _, _, t in detail_cols], "n_rows": n_rows, "n_nan": n_nan,
            "n_af_filtered": n_af}


def read_min_txt(path, test, min_col=None):
    col = min_col or MIN_COL_BY_TEST.get(test)
    if col is None:
        log("  WARN %s: no min.txt p column known for test %s; set min_col" % (rel(path, ROOT), test))
        return None, None
    pvals = []
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            expect = MIN_COL_BY_TEST.get(test)
            if expect and min_col is None and len(parts) != expect:
                log("  WARN %s: %d columns, expected %d for %s" % (rel(path, ROOT), len(parts), expect, test))
            if len(parts) >= col:
                v = parse_float(parts[col - 1])
                if v == v:
                    pvals.append(v)
    return thresholds_from_pvals(pvals)


def read_peaks(path, filtered, lengths, label):
    peaks = []
    if not path:
        return peaks
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[1].strip().lstrip("-").isdigit():
                continue
            contig = bed_contig(parts[0])
            if contig not in lengths:
                log("  WARN %s: peak contig %s not in .fai; skipped" % (label, contig))
                continue
            rec = {"contig": contig, "start": int(parts[1]), "end": int(parts[2])}
            if filtered:
                rec["id"] = parts[3].strip() if len(parts) > 3 else str(len(peaks) + 1)
                rec["max_logp"] = round(parse_float(parts[4]), 4) if len(parts) > 4 else None
                rec["lead"] = int(float(parts[5])) if len(parts) > 5 and parts[5].strip() else None
            peaks.append(rec)
    return peaks


def detail_value(raw, typ):
    s = raw.decode()
    if typ is str:
        return s
    if typ is int:
        try:
            return int(s)
        except ValueError:
            return None
    v = parse_float(s)
    return None if v != v or math.isinf(v) else v


def gwas_plan(track, lengths):
    """Resolve paths and the build signature for a GWAS track."""
    path = os.path.join(ROOT, track["path"])
    min_txt = resolve_auto(track.get("min_txt", "auto"), path, "_permution/min.txt")
    peaks = resolve_auto(track.get("peaks", "auto"), path, "_PEAKS.bed")
    peaks_u = resolve_auto(track.get("peaks_unfiltered", "auto"), path, "_PEAKS_unfiltered.bed")
    sig = {"v": PREP_VERSION, "track": track, "source": file_sig(path), "min_txt": opt_sig(min_txt),
           "peaks": opt_sig(peaks), "peaks_u": opt_sig(peaks_u), "fai": FAI_SIG}
    return path, min_txt, peaks, peaks_u, options_hash(sig)


def build_gwas(track, out_dir, lengths, force=False):
    t0 = time.time()
    tid = track["id"]
    final = os.path.join(out_dir, "gwas", tid)
    path, min_txt, peaks_path, peaks_u_path, sig = gwas_plan(track, lengths)
    meta = read_json(os.path.join(final, "meta.json"))
    if not force and meta and meta.get("sig") == sig:
        return {"id": tid, "skipped": True, "entry": meta["entry"], "secs": time.time() - t0}

    parsed = parse_gwas(path, track.get("p_col", "auto"), float(track.get("min_af") or 0))
    p_name = parsed["p_name"]
    test = TEST_BY_PCOL.get(p_name.lower(), "gwas")
    os.makedirs(os.path.join(out_dir, "gwas"), exist_ok=True)
    tmp = tempfile.mkdtemp(prefix=tid + ".", dir=os.path.join(out_dir, "gwas"))
    os.makedirs(os.path.join(tmp, "full"))
    os.makedirs(os.path.join(tmp, "detail"))

    names = sorted((c for c in parsed["contigs"] if c in lengths), key=contig_sort_key)
    for c in parsed["contigs"]:
        if c not in lengths:
            log("  WARN %s: contig %s not in .fai; %d SNPs skipped"
                % (tid, c, len(parsed["contigs"][c][0])))
    per = {}
    n_snps = n_capped = 0
    for c in names:
        pos, p, det = parsed["contigs"][c]
        if len(pos) > 1 and np.any(np.diff(pos.astype(np.int64)) < 0):
            o = np.argsort(pos, kind="stable")
            pos, p = pos[o], p[o]
        scaled, capped = scale_logp(p)
        per[c] = (pos, scaled, det)
        n_snps += len(pos)
        n_capped += capped

    # overview cutoff: start at 2, raise in 0.5 steps while too many SNPs pass
    cutoff = OVERVIEW_CUTOFF
    while True:
        thr = int(round(cutoff * LOGP_SCALE))
        n_ov = sum(int(np.count_nonzero(v[1] >= thr)) for v in per.values())
        if n_ov <= OVERVIEW_MAX or cutoff >= 60:
            break
        cutoff += 0.5
    det_thr = int(round(DETAIL_CUTOFF * LOGP_SCALE))

    ov_pos, ov_lp, bins_all, contig_meta = [], [], [], {}
    ov_off = bin_off = 0
    max_all = 0
    types = parsed["detail_types"]
    fields = parsed["detail_fields"]
    for c in names:
        pos, scaled, det = per[c]
        atomic_write_bytes(os.path.join(tmp, "full", c + ".bin"),
                           pos.astype("<u4").tobytes() + scaled.astype("<u2").tobytes())
        keep = scaled >= thr
        ov_pos.append(pos[keep])
        ov_lp.append(scaled[keep])
        n_bins = math.ceil(lengths[c] / BIN_SIZE)
        bins = np.zeros(n_bins, dtype=np.uint16)
        if len(pos):
            idx = np.minimum((pos.astype(np.int64) - 1) // BIN_SIZE, n_bins - 1)
            idx = np.maximum(idx, 0)
            np.maximum.at(bins, idx, scaled)
        bins_all.append(bins)
        cmax = int(scaled.max()) if len(scaled) else 0
        max_all = max(max_all, cmax)
        contig_meta[c] = {"n": int(len(pos)), "ov_off": ov_off, "ov_n": int(keep.sum()),
                          "bin_off": bin_off, "bin_n": n_bins, "max_logp": round(cmax / LOGP_SCALE, 3)}
        ov_off += int(keep.sum())
        bin_off += n_bins

        # detail rows: candidates kept at parse time, final cut on the scaled value
        if det:
            dpos = np.fromiter((d[0] for d in det), dtype=np.int64, count=len(det))
            dp = np.array([parse_float(d[1][fields.index("p")].decode()) for d in det])
            dscaled, _ = scale_logp(dp)
            sel = np.nonzero(dscaled >= det_thr)[0]
            if len(sel):
                sel = sel[np.argsort(dpos[sel], kind="stable")]
                cols = {"pos": [int(dpos[i]) for i in sel]}
                for j, (f, typ) in enumerate(zip(fields, types)):
                    cols[f] = [detail_value(det[i][1][j], typ) for i in sel]
                atomic_write_json(os.path.join(tmp, "detail", c + ".json"), cols)

    pos_cat = np.concatenate(ov_pos) if ov_pos else np.zeros(0, np.uint32)
    lp_cat = np.concatenate(ov_lp) if ov_lp else np.zeros(0, np.uint16)
    atomic_write_bytes(os.path.join(tmp, "overview.bin"),
                       pos_cat.astype("<u4").tobytes() + lp_cat.astype("<u2").tobytes())
    bins_cat = np.concatenate(bins_all) if bins_all else np.zeros(0, np.uint16)
    atomic_write_bytes(os.path.join(tmp, "bins.bin"), bins_cat.astype("<u2").tobytes())

    sig_v = track.get("sig")
    sug_v = track.get("sug")
    thresh_source = None
    if sig_v is None and sug_v is None and min_txt:
        if test in MIN_COL_BY_TEST or track.get("min_col"):
            sig_v, sug_v = read_min_txt(min_txt, test, track.get("min_col"))
            thresh_source = rel(min_txt, ROOT)
        else:
            log("  WARN %s: .GWAS source; min.txt ignored without min_col" % tid)
    elif sig_v is not None or sug_v is not None:
        thresh_source = "manual"
    elif track.get("min_txt", "auto"):
        log("  WARN %s: no min.txt found; no threshold lines" % tid)

    peaks = read_peaks(peaks_path, True, lengths, tid)
    peaks_u = read_peaks(peaks_u_path, False, lengths, tid)
    entry = {
        "id": tid, "type": "gwas", "label": track.get("label") or tid, "source": track["path"],
        "test": test, "p_col": p_name, "n_rows": parsed["n_rows"], "n_snps": n_snps,
        "n_nan": parsed["n_nan"], "n_capped": n_capped, "min_af": float(track.get("min_af") or 0),
        "n_af_filtered": parsed["n_af_filtered"],
        "overview_cutoff": cutoff, "detail_cutoff": DETAIL_CUTOFF,
        "sig": None if sig_v is None else round(float(sig_v), 4),
        "sug": None if sug_v is None else round(float(sug_v), 4),
        "thresh_source": thresh_source, "detail_fields": fields,
        "peaks": peaks, "peaks_unfiltered": peaks_u,
        "max_logp": round(max_all / LOGP_SCALE, 3), "contigs": contig_meta,
    }
    atomic_write_json(os.path.join(tmp, "meta.json"), {"sig": sig, "entry": entry}, indent=1)
    swap_dir(tmp, final)
    return {"id": tid, "skipped": False, "entry": entry, "secs": time.time() - t0}


# ------------------------------------------------------------------------------------------- line tracks

def line_paths(track):
    pattern = os.path.join(ROOT, track["path"])
    paths = sorted(glob.glob(pattern)) if any(ch in track["path"] for ch in "*?[") else [pattern]
    if not paths:
        raise FileNotFoundError("no files match %s" % track["path"])
    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(p)
    return paths


def line_spec(track):
    preset = track.get("preset")
    if preset not in LINE_PRESETS:
        raise ValueError("track %s: unknown preset %r" % (track["id"], preset))
    spec = dict(LINE_PRESETS[preset])
    if track.get("value_col"):
        spec["value"] = track["value_col"]
    if preset == "dsuite" and spec["value"] not in ("f_dM", "f_d", "D"):
        spec.pop("ymin", None)
        spec.pop("ymax", None)
    return spec


def filter_key(track, spec):
    flt = track.get("filter") or {}
    if spec["filter"] == "pop":
        v = flt.get("pop") or flt.get("pop1")
        if not v:
            raise ValueError("track %s: preset %s needs filter.pop" % (track["id"], track["preset"]))
        return ("pop", v)
    if spec["filter"] == "pop1":
        v = flt.get("pop1") or flt.get("pop")
        if not v:
            raise ValueError("track %s: preset %s needs filter.pop1" % (track["id"], track["preset"]))
        return ("pop", v)
    if spec["filter"] == "pair":
        a, b = flt.get("pop1"), flt.get("pop2")
        if not (a and b):
            raise ValueError("track %s: preset %s needs filter.pop1 and pop2" % (track["id"], track["preset"]))
        return ("pair", frozenset((a, b)))
    return None


def line_sig(track):
    return options_hash({"v": PREP_VERSION, "track": track, "files": [file_sig(p) for p in line_paths(track)],
                         "fai": FAI_SIG})


def build_line_group(tracks, out_dir, lengths, force=False):
    """Build every line track that reads the same source file(s) in one pass."""
    t0 = time.time()
    line_dir = os.path.join(out_dir, "line")
    os.makedirs(line_dir, exist_ok=True)
    results, todo = [], []
    for tr in tracks:
        sig = line_sig(tr)
        meta = read_json(os.path.join(line_dir, tr["id"] + ".meta.json"))
        if not force and meta and meta.get("sig") == sig:
            results.append({"id": tr["id"], "skipped": True, "entry": meta["entry"], "secs": 0.0})
        else:
            todo.append((tr, sig))
    if not todo:
        return results

    spec = line_spec(todo[0][0])
    by_key = {}
    data = {}
    for tr, _ in todo:
        k = filter_key(tr, spec)
        by_key.setdefault(k, []).append(tr["id"])
        data[tr["id"]] = {}     # contig -> (array start, array end, array value)
    window = {tr["id"]: int(tr.get("window", 5000)) for tr, _ in todo}
    tid_any = todo[0][0]["id"]
    skipped_contigs = set()

    for path in line_paths(todo[0][0]):
        with open_text(path) as fh:
            header = fh.readline().decode().rstrip("\r\n").split("\t")
            hidx = {h: i for i, h in enumerate(header)}
            need = [spec["chr"], spec["start"], spec["value"]] + ([spec["end"]] if spec["end"] else [])
            for col in need:
                if col not in hidx:
                    raise ValueError("%s: column %s missing (header %s)" % (rel(path, ROOT), col, header))
            ci, si, vi = hidx[spec["chr"]], hidx[spec["start"]], hidx[spec["value"]]
            ei = hidx[spec["end"]] if spec["end"] else None
            if spec["filter"] == "pair":
                f1, f2 = hidx["pop1"], hidx["pop2"]
            elif spec["filter"] == "pop":
                f1 = hidx["pop"]
            elif spec["filter"] == "pop1":
                f1 = hidx["pop1"]
            # only fully split rows whose filter columns match
            wanted_pops = set()
            for k in by_key:
                if k and k[0] == "pop":
                    wanted_pops.add(k[1].encode())
                elif k and k[0] == "pair":
                    wanted_pops.update(x.encode() for x in k[1])
            for line in fh:
                parts = line.rstrip(b"\r\n").split(b"\t")
                if len(parts) < len(header):
                    continue
                if spec["filter"] is None:
                    ids = by_key[None]
                elif spec["filter"] == "pair":
                    a, b = parts[f1], parts[f2]
                    if a not in wanted_pops or b not in wanted_pops:
                        continue
                    ids = by_key.get(("pair", frozenset((a.decode(), b.decode()))))
                    if not ids:
                        continue
                else:
                    a = parts[f1]
                    if a not in wanted_pops:
                        continue
                    ids = by_key.get(("pop", a.decode()))
                    if not ids:
                        continue
                c = parts[ci].decode()
                if c not in lengths:
                    skipped_contigs.add(c)
                    continue
                v = parse_float(parts[vi].decode())
                if spec["kind"] == "point":
                    s = max(1, int(round(float(parts[si]))))
                    e = s
                else:
                    s = int(float(parts[si]))
                    if ei is None:     # tajd: 0-based BIN_START
                        s += 1
                for tid in ids:
                    if spec["kind"] != "point" and ei is None:
                        e = s + window[tid] - 1
                    elif ei is not None:
                        e = int(float(parts[ei]))
                    d = data[tid].get(c)
                    if d is None:
                        d = data[tid][c] = (array.array("I"), array.array("I"), array.array("f"))
                    d[0].append(s)
                    d[1].append(e)
                    d[2].append(v)
    if skipped_contigs:
        log("  WARN %s: contigs not in .fai skipped: %s" % (tid_any, ", ".join(sorted(skipped_contigs))))

    for tr, sig in todo:
        tid = tr["id"]
        names = sorted(data[tid], key=contig_sort_key)
        starts, ends, vals, cmeta = [], [], [], {}
        off = 0
        for c in names:
            s = np.frombuffer(data[tid][c][0], dtype=np.uint32)
            e = np.frombuffer(data[tid][c][1], dtype=np.uint32)
            v = np.frombuffer(data[tid][c][2], dtype=np.float32)
            o = np.argsort(s, kind="stable")
            starts.append(s[o]); ends.append(e[o]); vals.append(v[o])
            cmeta[c] = {"off": off, "n": int(len(s))}
            off += len(s)
        s_all = np.concatenate(starts) if starts else np.zeros(0, np.uint32)
        e_all = np.concatenate(ends) if ends else np.zeros(0, np.uint32)
        v_all = np.concatenate(vals) if vals else np.zeros(0, np.float32)
        if off == 0:
            log("  WARN %s: no rows matched (check filter/path)" % tid)
        atomic_write_bytes(os.path.join(line_dir, tid + ".bin"),
                           s_all.astype("<u4").tobytes() + e_all.astype("<u4").tobytes()
                           + v_all.astype("<f4").tobytes())
        finite = v_all[np.isfinite(v_all)]
        entry = {
            "id": tid, "type": "line", "label": tr.get("label") or default_line_label(tr, spec),
            "preset": tr["preset"],
            "source": tr["path"], "value_name": spec["value"], "kind": spec["kind"],
            "ymin": tr["ymin"] if "ymin" in tr else spec.get("ymin"),
            "ymax": tr["ymax"] if "ymax" in tr else spec.get("ymax"),
            "zero_line": bool(tr.get("zero_line", spec.get("zero_line", False))),
            "vmin": float(finite.min()) if len(finite) else None,
            "vmax": float(finite.max()) if len(finite) else None,
            "n": int(off), "n_nan": int(len(v_all) - len(finite)), "contigs": cmeta,
        }
        atomic_write_json(os.path.join(line_dir, tid + ".meta.json"), {"sig": sig, "entry": entry}, indent=1)
        results.append({"id": tid, "skipped": False, "entry": entry, "secs": time.time() - t0})
    return results


def default_line_label(track, spec):
    flt = track.get("filter") or {}
    pre = track["preset"]
    if pre == "pixy_pi":
        return "pi %s" % (flt.get("pop") or flt.get("pop1"))
    if pre in ("pixy_fst", "pixy_dxy"):
        return "%s %s-%s" % ("Fst" if pre == "pixy_fst" else "dxy", flt.get("pop1"), flt.get("pop2"))
    if pre == "tajd":
        return "Tajima's D %s" % (flt.get("pop1") or flt.get("pop"))
    if pre == "dsuite":
        return "%s %s" % (spec["value"], os.path.basename(track["path"]).split("_localFstats")[0])
    return "CLR %s %s" % (pre, os.path.basename(track["path"]).split(".")[0])


# ------------------------------------------------------------------------------------------------- genes

def parse_attrs(s):
    out = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def get_attr(s, key):
    i = s.find(key + "=")
    while i > 0 and s[i - 1] != ";":
        i = s.find(key + "=", i + 1)
    if i < 0:
        return None
    j = s.find(";", i)
    return s[i + len(key) + 1: j if j >= 0 else len(s)]


def build_genes(gff_path, out_dir, lengths, force=False):
    t0 = time.time()
    final = os.path.join(out_dir, "genes")
    sig = options_hash({"v": PREP_VERSION, "gff": file_sig(gff_path), "fai": FAI_SIG})
    meta = read_json(os.path.join(final, "meta.json"))
    if not force and meta and meta.get("sig") == sig:
        return {"id": "genes", "skipped": True, "entry": meta["entry"], "secs": time.time() - t0}

    genes = {}          # gene id -> record
    parent_of = {}      # feature id -> parent id (transcripts and nested features)
    exons, cds = {}, {}
    skip_types = {"region", "cDNA_match", "match", "sequence_feature"}
    with open(gff_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype in ("exon", "CDS"):
                par = get_attr(parts[8], "Parent")
                if par:
                    target = exons if ftype == "exon" else cds
                    se = (int(parts[3]), int(parts[4]))
                    for pid in par.split(","):
                        target.setdefault(pid, []).append(se)
                continue
            if ftype in skip_types:
                continue
            a = parse_attrs(parts[8])
            fid = a.get("ID")
            if not fid:
                continue
            if ftype in ("gene", "pseudogene"):
                genes[fid] = {"contig": parts[0], "s": int(parts[3]), "e": int(parts[4]), "st": parts[6],
                              "n": a.get("Name") or a.get("gene") or fid, "id": fid,
                              "b": a.get("gene_biotype") or ftype, "d": a.get("description", "")}
            elif "Parent" in a:
                parent_of[fid] = a["Parent"].split(",")[0]

    def gene_of(fid):
        seen = 0
        while fid in parent_of and seen < 10:
            fid = parent_of[fid]
            seen += 1
        return fid if fid in genes else None

    tx_by_gene = {}
    for fid in parent_of:
        if fid in exons or fid in cds:
            g = gene_of(fid)
            if g:
                tx_by_gene.setdefault(g, []).append(fid)
    for gid in genes:
        if gid in exons or gid in cds:      # pseudogene exons parented directly by the gene
            tx_by_gene.setdefault(gid, []).append(gid)

    per_contig = {}
    n_skipped = 0
    for gid, g in genes.items():
        if g["contig"] not in lengths:
            n_skipped += 1
            continue
        best, best_score = None, (-1, -1)
        for t in tx_by_gene.get(gid, []):
            ex = exons.get(t) or cds.get(t) or []
            score = (sum(e - s + 1 for s, e in cds.get(t, [])), sum(e - s + 1 for s, e in ex))
            if score > best_score:
                best, best_score = t, score
        x, cs, ce = [], None, None
        if best:
            for s, e in sorted(exons.get(best) or cds.get(best) or []):
                x += [s, e]
            if cds.get(best):
                cs = min(s for s, _ in cds[best])
                ce = max(e for _, e in cds[best])
        rec = {"s": g["s"], "e": g["e"], "st": g["st"], "n": g["n"], "id": gid, "b": g["b"], "d": g["d"],
               "t": best if best and best != gid else None, "x": x, "cs": cs, "ce": ce}
        per_contig.setdefault(g["contig"], []).append(rec)
    if n_skipped:
        log("  WARN genes: %d genes on contigs not in .fai skipped" % n_skipped)

    os.makedirs(out_dir, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="genes.", dir=out_dir)
    names = sorted(per_contig, key=contig_sort_key)
    rows = ["name\tcontig\tstart\tend\tstrand\tbiotype\tdescription"]
    clean = lambda s: str(s).replace("\t", " ").replace("\n", " ").replace("\r", " ")
    n_genes = 0
    for c in names:
        recs = sorted(per_contig[c], key=lambda r: (r["s"], r["e"]))
        atomic_write_json(os.path.join(tmp, c + ".json"), recs)
        for r in recs:
            rows.append("\t".join(clean(v) for v in (r["n"], c, r["s"], r["e"], r["st"], r["b"], r["d"])))
        n_genes += len(recs)
    atomic_write_bytes(os.path.join(tmp, "index.tsv"), ("\n".join(rows) + "\n").encode())
    entry = {"index": "genes/index.tsv", "contigs": names}
    atomic_write_json(os.path.join(tmp, "meta.json"), {"sig": sig, "entry": entry, "n_genes": n_genes})
    swap_dir(tmp, final)
    return {"id": "genes", "skipped": False, "entry": entry, "secs": time.time() - t0, "n_genes": n_genes}


# -------------------------------------------------------------------------------------------- orchestration

def _init_worker(root, fai_sig):
    global ROOT, FAI_SIG
    ROOT, FAI_SIG = root, fai_sig


def run_job(kind, payload, out_dir, lengths, force):
    if kind == "gwas":
        return [build_gwas(payload, out_dir, lengths, force)]
    if kind == "line":
        return build_line_group(payload, out_dir, lengths, force)
    if kind == "genes":
        return [build_genes(payload, out_dir, lengths, force)]
    raise ValueError(kind)


def describe(res):
    e = res["entry"]
    if res["id"] == "genes":
        return "genes: %s contigs%s" % (len(e["contigs"]), ", %d genes" % res["n_genes"] if "n_genes" in res else "")
    if e["type"] == "gwas":
        return ("%s: %d rows, %d SNPs, n_nan=%d, n_capped=%d, af_filtered=%d, sig=%s sug=%s (%s), "
                "overview_cutoff=%.1f, peaks=%d/%d, max_logp=%.3f"
                % (e["id"], e["n_rows"], e["n_snps"], e["n_nan"], e["n_capped"], e["n_af_filtered"], e["sig"],
                   e["sug"], e["thresh_source"], e["overview_cutoff"], len(e["peaks"]),
                   len(e["peaks_unfiltered"]), e["max_logp"]))
    return "%s: %d windows, n_nan=%d, %s range %s..%s" % (e["id"], e["n"], e["n_nan"], e["value_name"],
                                                         e["vmin"], e["vmax"])


def prepare(manifest_path, out_dir, only=None, force=False, jobs=None):
    global FAI_SIG
    t_start = time.time()
    with open(manifest_path) as fh:
        man = json.load(fh)
    genome = man.get("genome") or {}
    fai = os.path.join(ROOT, genome.get("fai") or DEFAULT_FAI)
    lengths = read_fai(fai)
    FAI_SIG = file_sig(fai)
    gff = genome.get("gff")
    gff = os.path.join(ROOT, gff) if gff else None
    tracks = man.get("tracks", [])

    ids = [t.get("id") for t in tracks]
    for t in tracks:
        if not t.get("id") or not TRACK_ID_RE.match(t["id"]):
            sys.exit("bad track id %r (use letters, digits, _ and -)" % t.get("id"))
        if t.get("type") not in ("gwas", "line"):
            sys.exit("track %s: type must be gwas or line" % t["id"])
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        sys.exit("duplicate track ids: %s" % ", ".join(sorted(dup)))
    only = set(only) if only else None

    os.makedirs(out_dir, exist_ok=True)
    atomic_write_bytes(os.path.join(out_dir, ".gitignore"), b"*\n")

    jobs_list = []
    if gff and only is None:
        jobs_list.append(("genes", gff))
    line_groups = {}
    for t in tracks:
        if only is not None and t["id"] not in only:
            continue
        if t["type"] == "gwas":
            jobs_list.append(("gwas", t))
        else:
            key = (t.get("preset"), t["path"])
            line_groups.setdefault(key, []).append(t)
    for grp in line_groups.values():
        jobs_list.append(("line", grp))

    results = {}
    njobs = max(1, min(jobs or min(4, os.cpu_count() or 1), len(jobs_list) or 1))
    failed = []
    with ProcessPoolExecutor(max_workers=njobs, initializer=_init_worker, initargs=(ROOT, FAI_SIG)) as ex:
        futs = {ex.submit(run_job, k, p, out_dir, lengths, force): (k, p) for k, p in jobs_list}
        for f in as_completed(futs):
            k, p = futs[f]
            try:
                for res in f.result():
                    results[res["id"]] = res
                    tag = "skip " if res["skipped"] else "built"
                    log("%s %6.1fs  %s" % (tag, res["secs"], describe(res)))
            except Exception as err:   # report and continue with other tracks
                name = p if isinstance(p, str) else (p["id"] if isinstance(p, dict) else
                                                     ",".join(t["id"] for t in p))
                log("FAIL %s %s: %s" % (k, name, err))
                failed.append(name)

    # manifest from all track metas, in manifest order
    entries = []
    for t in tracks:
        if t["id"] in results:
            entries.append(results[t["id"]]["entry"])
            continue
        mpath = (os.path.join(out_dir, "gwas", t["id"], "meta.json") if t["type"] == "gwas"
                 else os.path.join(out_dir, "line", t["id"] + ".meta.json"))
        meta = read_json(mpath)
        if meta:
            entries.append(meta["entry"])
        else:
            log("WARN %s: not built yet; left out of manifest.json" % t["id"])
    genes_entry = None
    if gff:
        if "genes" in results:
            genes_entry = results["genes"]["entry"]
        else:
            meta = read_json(os.path.join(out_dir, "genes", "meta.json"))
            genes_entry = meta["entry"] if meta else None
    manifest = {
        "version": 1,
        "built": datetime.datetime.now().replace(microsecond=0).isoformat(),
        "logp_scale": LOGP_SCALE,
        "bin_size": BIN_SIZE,
        "contigs": [{"name": c, "length": lengths[c]} for c in sorted(lengths, key=contig_sort_key)],
        "genes": genes_entry,
        "tracks": entries,
    }
    atomic_write_json(os.path.join(out_dir, "manifest.json"), manifest)
    log("manifest.json: %d tracks; total %.1fs%s" % (len(entries), time.time() - t_start,
                                                     "; FAILED: " + ", ".join(failed) if failed else ""))
    return 1 if failed else 0


def discover():
    found = sorted(glob.glob(os.path.join(ROOT, "output_data", "**", "*.assoc.txt"), recursive=True))

    def rank(p):
        r = rel(p, ROOT)
        return (2 if "/set_aside/" in r else 1 if "14_Wobble" in r else 0, r)

    tracks = []
    for p in sorted(found, key=rank):
        r = rel(p, ROOT)
        stem = os.path.basename(p)[: -len(".assoc.txt")]
        tid = re.sub(r"[^A-Za-z0-9_-]", "_", stem)
        try:
            header = [h.lower() for h in read_header(p)]
        except OSError:
            continue
        test = "lmm1" if "p_wald" in header else "lmm2" if "p_lrt" in header else "lmm"
        parts = stem.split(".")
        label = parts[0] + (" " + " ".join(x for x in parts[1:] if not x.startswith("lmm")) if len(parts) > 2
                            else "")
        label = "%s (%s)%s" % (label.strip(), test, " [set aside]" if "/set_aside/" in r else "")
        tracks.append({"id": tid, "label": label, "type": "gwas", "path": r, "p_col": "auto",
                       "min_txt": "auto", "peaks": "auto", "peaks_unfiltered": "auto", "min_af": 0,
                       "sig": None, "sug": None})
    return {"genome": {"fai": DEFAULT_FAI, "gff": DEFAULT_GFF}, "tracks": tracks}


# ------------------------------------------------------------------------------------------------ selftest

def selftest():
    ok = True

    def check(name, cond):
        nonlocal ok
        log(("ok    " if cond else "FAIL  ") + name)
        ok = ok and bool(cond)

    # threshold formula: n=100 -> 5th and 33rd smallest p (1-based), as in gwas_report.Rmd
    pv = [10 ** -(i / 10) for i in range(100)]          # 1 .. 1e-9.9
    srt = sorted(pv)
    sig, sug = thresholds_from_pvals(pv)
    check("threshold n=100 sig uses 5th smallest", abs(sig - -math.log10(srt[4])) < 1e-12)
    check("threshold n=100 sug uses 33rd smallest", abs(sug - -math.log10(srt[32])) < 1e-12)
    sig, sug = thresholds_from_pvals(pv[:69])
    # n=69: 69 - floor(65.55) = 4 ; 69 - floor(46.23) = 23
    check("threshold n=69 indices 4/23",
          abs(sig - -math.log10(sorted(pv[:69])[3])) < 1e-12 and abs(sug - -math.log10(sorted(pv[:69])[22])) < 1e-12)
    real = os.path.join(ROOT, "output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/"
                              "PKINGS_ALL_WOB_EXCLUDED_permution/min.txt")
    if os.path.exists(real):
        s, u = read_min_txt(real, "lmm2")
        check("PKINGS_ALL_WOB_EXCLUDED sig 7.2223 sug 6.3572 (got %.4f %.4f)" % (s, u),
              round(s, 4) == 7.2223 and round(u, 4) == 6.3572)

    check("bed chr 24 -> chr24", bed_contig("24") == "chr24")
    check("bed chr chr5 kept", bed_contig("chr5") == "chr5")
    check("bed chr ups3 kept", bed_contig("ups3") == "ups3")
    order = sorted(["ups2", "chr10", "chr2", "chr1", "ups10", "scaffold_x", "chrM"], key=contig_sort_key)
    check("display order %s" % order, order == ["chr1", "chr2", "chr10", "ups2", "ups10", "chrM", "scaffold_x"])

    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "T.lmm2.assoc.txt")
        rows = ["chr\trs\tps\tn_miss\tallele1\tallele0\taf\tlogl_H1\tl_mle\tp_lrt",
                "chr2\ta\t10\t0\tA\tG\t0.2\t0\t0\t0.5",
                "chr2\tb\t20\t0\tA\tG\t0.2\t0\t0\tnan",
                "chr2\tc\t30\t0\tA\tG\t0.01\t0\t0\t0",
                "chr2\td\t40\t0\tA\tG\t0.3\t0\t0\t-nan",
                "chr1\te\t5\t1\tC\tT\t0.4\t0\t0\t1.000000e-05",
                "chr1\tf\t7\t1\tC\tT\t0.4\t0\t0\t0.0099"]
        with open(f, "w") as fh:
            fh.write("\n".join(rows) + "\n")
        r = parse_gwas(f)
        check("NaN rows skipped (n_nan=2)", r["n_nan"] == 2 and r["n_rows"] == 6)
        pos, p, det = r["contigs"]["chr2"]
        check("chr2 kept 2 SNPs", list(pos) == [10, 30])
        scaled, capped = scale_logp(p)
        check("p=0 capped at 65535", int(scaled[1]) == U16_MAX and capped == 1)
        check("p=0.5 scaled 301", int(scaled[0]) == 301)
        check("detail candidates: p=0, 1e-5 and 0.0099 kept",
              len(det) == 1 and len(r["contigs"]["chr1"][2]) == 2)
        r2 = parse_gwas(f, min_af=0.05)
        check("min_af drops af<0.05", r2["n_af_filtered"] == 1 and list(r2["contigs"]["chr2"][0]) == [10])
        check("p column auto = p_lrt", r["p_name"] == "p_lrt")
        check("detail fields lmm2", r["detail_fields"] == ["rs", "p", "af", "n_miss", "allele1", "allele0"])
    log("selftest " + ("passed" if ok else "FAILED"))
    return 0 if ok else 1


ROOT = repo_root()
FAI_SIG = None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", nargs="?", help="track manifest JSON (paths relative to the repo root)")
    ap.add_argument("--out", default=os.path.join(ROOT, "output_data", "gwas_browser"))
    ap.add_argument("--only", help="comma-separated track ids to (re)build; others keep their last build")
    ap.add_argument("--force", action="store_true", help="rebuild even if inputs are unchanged")
    ap.add_argument("--jobs", type=int, default=None, help="parallel workers (default min(4, cpus))")
    ap.add_argument("--discover", action="store_true", help="print a manifest of all output_data/**/*.assoc.txt")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()
    if args.discover:
        json.dump(discover(), sys.stdout, indent=1)
        sys.stdout.write("\n")
        return 0
    if not args.manifest:
        ap.error("manifest required (or --discover / --selftest)")
    only = [x for x in args.only.split(",") if x] if args.only else None
    return prepare(os.path.abspath(args.manifest), os.path.abspath(args.out), only, args.force, args.jobs)


if __name__ == "__main__":
    sys.exit(main())
