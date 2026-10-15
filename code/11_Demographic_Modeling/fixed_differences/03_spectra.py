#!/usr/bin/env python3
"""03_spectra.py
Turn one chunk of the site table into the quantities DESIGN.md sections 3-4 need, in three passes so
that every rule fixed in DESIGN.md is applied in the order it depends on:

  --mode an_hist   AN distributions of masked, polarized SNPs per pair and side (for the projection
                   rule, applied genome-wide by 04_projections.py)
  --mode bins      per 100-kb bin and pair: callable sites and callable SNPs polymorphic in the pair, at
                   the chosen projection sizes (for the call-desert mask, 05_deserts.py)
  --mode spectra   projected unfolded 2D spectra per (1-Mb block, exon class, region), the monomorphic
                   counts that complete them, a fine histogram of raw delta for empirical ranks, the
                   near-fixed SNPs for clustering, and candidate rows for the parallel count

Definitions (DESIGN.md section 3-4)
  derived count d = AC if the outgroup carries REF, AN - AC if it carries ALT; unpolarized SNPs dropped
  raw delta    = d_bp/AN_bp - d_tp/AN_tp
  exon class   c3: 0 inside an exon, 1 within 5 kb of one, 2 farther (headline fit mask: c3 >= 1)
  region       "none", "peak1".."peak6" (extended peak +/-250 kb) or "chr4_inversion"
  paralog filter: SNPs with ExcHet < 0.001 in any whole population are dropped from data and L alike
  callable     AN >= projection size in both populations; invariant sites also need the outgroup called

Pairs: by default the headline whole-population pairs and the TP-TP controls; --pairs_role selects others.

Usage (inside dadi_tools.sif, from the repo root):
  python3 .../03_spectra.py --mode an_hist --chunk_id 17
"""
import argparse
import collections
import gzip
import os
import pickle
import sys

import numpy as np

ROOT = os.getcwd()
INP = os.path.join(ROOT, "input_data/11_Demographic_Modeling")
OUT = os.path.join(ROOT, os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
ST = os.path.join(OUT, "site_table")
EXCHET_MIN = 0.001
DELTA_BINS = np.linspace(-1, 1, 801)        # 0.0025 resolution for empirical ranks


def to_int(x):
    return 0 if x in (".", "") else int(x)


def c3_of(dbin):
    return 0 if dbin == 0 else (1 if dbin <= 2 else 2)


def load_pairs(roles):
    pairs = []
    with open(os.path.join(INP, "pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            pid, bp, tp, role = line.rstrip("\n").split("\t")
            if role in roles:
                pairs.append((pid, bp, tp))
    return pairs


def load_regions(chrom):
    regs = []
    with open(os.path.join(OUT, "masks/exclude.bed")) as fh:
        for line in fh:
            c, s, e, name = line.rstrip("\n").split("\t")
            if c == chrom:
                regs.append((int(s), int(e), name))
    return regs


def region_of(regs, x):
    for s, e, name in regs:
        if s <= x < e:
            return name
    return "none"


_hyp = {}


def hyper(N, d, n):
    """P(k derived among n drawn without replacement from N haplotypes carrying d derived), k = 0..n."""
    key = (N, d, n)
    v = _hyp.get(key)
    if v is None:
        from math import comb
        denom = comb(N, n)
        v = np.array([comb(d, k) * comb(N - d, n - k) / denom for k in range(n + 1)])
        _hyp[key] = v
    return v


def read_projections(suffix=""):
    proj = {}
    with open(os.path.join(OUT, f"sfs/projections{suffix}.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            proj[f[0]] = (int(f[1]), int(f[2]))
    return proj


def read_deserts(suffix=""):
    d = collections.defaultdict(set)
    path = os.path.join(OUT, f"sfs/deserts{suffix}.tsv")
    with open(path) as fh:
        next(fh)
        for line in fh:
            pid, chrom, b = line.rstrip("\n").split("\t")[:3]
            d[pid].add((chrom, int(b)))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["an_hist", "bins", "spectra"])
    ap.add_argument("--chunk_id", required=True)
    ap.add_argument("--pairs_role", default="headline_whole,tp_control_whole",
                    help="comma-separated pair roles from pairs.tsv")
    ap.add_argument("--suffix", default="", help="output suffix for non-default pair sets")
    a = ap.parse_args()

    pairs = load_pairs(set(a.pairs_role.split(",")))
    if not pairs:
        sys.exit("no pairs selected")
    snp_path = os.path.join(ST, f"snps.{a.chunk_id}.tsv.gz")
    fh = gzip.open(snp_path, "rt")
    header = fh.readline().lstrip("#").rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(header)}
    exchet_cols = [ix[c] for c in header if c.startswith("ExcHet_")]
    cols = {pid: (ix["AN_" + bp], ix["AC_" + bp], ix["AN_" + tp], ix["AC_" + tp]) for pid, bp, tp in pairs}

    proj = read_projections(a.suffix) if a.mode != "an_hist" else {}
    deserts = read_deserts(a.suffix) if a.mode == "spectra" else {}

    an_hist = collections.Counter()                      # (pid, side, AN) -> n
    bins = collections.Counter()                         # (pid, chrom, bin, kind) -> n
    raw = {pid: collections.Counter() for pid, _, _ in pairs}   # (block, c3, region, Nb, db, Nt, dt) -> n
    dhist = {pid: collections.defaultdict(lambda: np.zeros(len(DELTA_BINS) - 1)) for pid, _, _ in pairs}
    near_fixed, parallel = [], []
    regs, chrom = None, None
    # local-Ne proxy for E_hetNe (DESIGN.md section 4): per 1-Mb block and whole population, the summed
    # per-site diversity and divergence to the outgroup over polarized SNPs; f_w is their ratio, so the
    # number of callable sites cancels
    whole8 = [g for g in ("wAPA", "wBAM", "wBAVA", "wBEN", "wBIK", "wBIR", "wDOG", "wDOV")]
    w_cols = {g: (ix["AN_" + g], ix["AC_" + g]) for g in whole8}
    blockstats = collections.defaultdict(lambda: [0.0, 0.0])

    # parallel count (P2, P3, P6): the same derived allele high in the southern BP creeks, low in TP2
    par_cols = {g: (ix["AN_" + g], ix["AC_" + g]) for g in ("wBIK", "wBAVA", "wBIR", "wTP2")}

    for line in fh:
        f = line.rstrip("\n").split("\t")
        if regs is None:
            chrom = f[0]
            regs = load_regions(chrom)
        aa = f[4]
        if aa == ".":
            continue
        if any(f[i] not in (".", "") and float(f[i]) < EXCHET_MIN for i in exchet_cols):
            continue
        pos = int(f[1]); x = pos - 1
        c3 = c3_of(int(f[5])); xcode = int(f[6])
        region = "none" if xcode == 0 else region_of(regs, x)
        b100 = x // 100_000
        for pid, (i_an, i_ac, j_an, j_ac) in cols.items():
            Nb, ab, Nt, at = to_int(f[i_an]), to_int(f[i_ac]), to_int(f[j_an]), to_int(f[j_ac])
            if a.mode == "an_hist":
                if region == "none" and c3 >= 1:
                    an_hist[(pid, "bp", Nb)] += 1
                    an_hist[(pid, "tp", Nt)] += 1
                continue
            n1, n2 = proj[pid]
            if Nb < n1 or Nt < n2:
                continue
            db = ab if aa == "R" else Nb - ab
            dt = at if aa == "R" else Nt - at
            if a.mode == "bins":
                bins[(pid, chrom, b100, "callable_snp_rows")] += 1
                if not ((db == 0 and dt == 0) or (db == Nb and dt == Nt)):
                    bins[(pid, chrom, b100, "polymorphic")] += 1
                continue
            if (chrom, b100) in deserts.get(pid, ()):
                continue
            raw[pid][(x // 1_000_000, c3, region, Nb, db, Nt, dt)] += 1
            delta = db / Nb - dt / Nt
            dhist[pid][(c3, region)][min(np.searchsorted(DELTA_BINS, delta, side="right") - 1, len(DELTA_BINS) - 2)] += 1
            if abs(delta) >= 0.85:
                near_fixed.append((pid, chrom, pos, db / Nb, dt / Nt, region, c3))
        if a.mode == "spectra":
            blk = x // 1_000_000
            for g, (ci, cj) in w_cols.items():
                an_g, ac_g = to_int(f[ci]), to_int(f[cj])
                if an_g >= 2:
                    d_g = ac_g if aa == "R" else an_g - ac_g
                    st = blockstats[(blk, g)]
                    st[0] += d_g * (an_g - d_g) / (an_g * (an_g - 1) / 2)
                    st[1] += d_g / an_g
            an_t, ac_t = to_int(f[par_cols["wTP2"][0]]), to_int(f[par_cols["wTP2"][1]])
            if an_t:
                ft = (ac_t if aa == "R" else an_t - ac_t) / an_t
                if ft <= 0.1:
                    fr = {}
                    for g in ("wBIK", "wBAVA", "wBIR"):
                        an_g, ac_g = to_int(f[par_cols[g][0]]), to_int(f[par_cols[g][1]])
                        fr[g] = ((ac_g if aa == "R" else an_g - ac_g) / an_g) if an_g else float("nan")
                    vals = [v for v in fr.values() if v == v]
                    if vals and max(vals) >= 0.5:
                        parallel.append((chrom, pos, fr["wBIK"], fr["wBAVA"], fr["wBIR"], ft, region, c3))

    os.makedirs(os.path.join(OUT, "sfs/chunks"), exist_ok=True)
    tag = f"{a.chunk_id}{a.suffix}"
    if a.mode == "an_hist":
        with open(os.path.join(OUT, f"sfs/chunks/an_hist.{tag}.tsv"), "w") as o:
            o.write("pair_id\tside\tAN\tn\n")
            for (pid, side, an), n in sorted(an_hist.items()):
                o.write(f"{pid}\t{side}\t{an}\t{n}\n")
        return

    # invariant sites from the site table histograms
    inv = collections.Counter()
    hist_path = os.path.join(ST, f"invariant_hist.{a.chunk_id}.tsv")
    with open(hist_path) as h:
        next(h)
        for line in h:
            pid, c, b100, dbin, xcode, out_ok, an_bp, an_tp, n = line.rstrip("\n").split("\t")
            if pid not in proj or out_ok != "1":
                continue
            n1, n2 = proj[pid]
            if int(an_bp) < n1 or int(an_tp) < n2:
                continue
            b100 = int(b100)
            if a.mode == "bins":
                bins[(pid, c, b100, "callable_invariant")] += int(n)
                continue
            if (c, b100) in deserts.get(pid, ()):
                continue
            # invariant tallies carry the exclusion code but not the peak id; peaks are excluded from fits
            # and L_scan counts every non-desert site, so the code is enough
            reg = "none" if xcode == "0" else ("chr4_inversion" if xcode == "2" else "peak")
            inv[(pid, b100 // 10, c3_of(int(dbin)), reg)] += int(n)

    if a.mode == "bins":
        with open(os.path.join(OUT, f"sfs/chunks/bins.{tag}.tsv"), "w") as o:
            o.write("pair_id\tchrom\tbin100k\tkind\tn\n")
            for (pid, c, b, kind), n in sorted(bins.items()):
                o.write(f"{pid}\t{c}\t{b}\t{kind}\t{n}\n")
        return

    # spectra: expand the raw tallies by hypergeometric projection
    spectra = {}
    for pid, _, _ in pairs:
        n1, n2 = proj[pid]
        sp = {}
        for (blk, c3, region, Nb, db, Nt, dt), n in raw[pid].items():
            key = (blk, c3, region)
            if key not in sp:
                sp[key] = np.zeros((n1 + 1, n2 + 1))
            sp[key] += n * np.outer(hyper(Nb, db, n1), hyper(Nt, dt, n2))
        spectra[pid] = sp
    payload = {
        "chrom": chrom, "chunk_id": a.chunk_id, "projections": {p: proj[p] for p, _, _ in pairs},
        "spectra": spectra,
        "invariant": dict(inv),
        "delta_bins": DELTA_BINS,
        "delta_hist": {pid: dict(v) for pid, v in dhist.items()},
        "blockstats": {k: tuple(v) for k, v in blockstats.items()},
    }
    with open(os.path.join(OUT, f"sfs/chunks/spectra.{tag}.pkl"), "wb") as o:
        pickle.dump(payload, o, protocol=pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(OUT, f"sfs/chunks/near_fixed.{tag}.tsv"), "w") as o:
        o.write("pair_id\tchrom\tpos\tf_bp\tf_tp\tregion\tc3\n")
        for r in near_fixed:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]:.4f}\t{r[5]}\t{r[6]}\n")
    with open(os.path.join(OUT, f"sfs/chunks/parallel.{tag}.tsv"), "w") as o:
        o.write("chrom\tpos\tf_BIK\tf_BAVA\tf_BIR\tf_TP2\tregion\tc3\n")
        for r in parallel:
            o.write("\t".join(str(v) if not isinstance(v, float) else f"{v:.4f}" for v in r) + "\n")


if __name__ == "__main__":
    main()
