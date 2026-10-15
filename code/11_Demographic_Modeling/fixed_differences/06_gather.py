#!/usr/bin/env python3
"""06_gather.py
Sum the per-chunk spectra into what the fits and the tail calculations use (DESIGN.md sections 3-5).

Per pair, writes sfs/<pair><suffix>.pkl holding
  data         projected unfolded spectrum over the fitting mask (region none, exon class in --mask_c3)
  L_mask       callable sites under that mask (spectrum total + invariant sites)
  L_scan       callable sites anywhere on chr1-25 outside call deserts (what the GWAS could scan)
  region_spec  spectra summed per region ("none", "peak1".."peak6", "chr4_inversion"), all exon classes
  blocks       1-Mb block keys, the per-block masked spectra and invariant counts (bootstrap)
  dhist        raw-delta histogram per region (empirical ranks)
  f_blocks     per 1-Mb block local-Ne factor f_w for E_hetNe, from the populations not in the pair
and sfs/<pair><suffix>.fs (dadi format, corners masked).

Also writes
  tail/index_units<suffix>.tsv     raw frequencies and delta at the six index SNPs for every pair
  fixed_diff_counts<suffix>.tsv    near-fixed and fixed differences per pair, mask and genome-wide

Usage (inside dadi_tools.sif, from the repo root):
  python3 .../06_gather.py [--suffix _x] [--mask_c3 1,2]
"""
import argparse
import collections
import glob
import os
import pickle
import re

import numpy as np
import dadi

ROOT = os.getcwd()
OUT = os.path.join(ROOT, os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"))
WHOLE8 = ["wAPA", "wBAM", "wBAVA", "wBEN", "wBIK", "wBIR", "wDOG", "wDOV"]
OWN = {"A1": [1, 4, 5], "A2": [2, 3, 6], "A3": [2, 3, 6], "A4": [2, 3]}


def to_int(x):
    return 0 if x in (".", "") else int(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="", help="suffix of the 03_spectra.py run to gather")
    ap.add_argument("--mask_c3", default="1,2", help="exon classes in the fitting mask (headline 1,2; +/-5 kb: 2)")
    ap.add_argument("--out_suffix", default=None, help="defaults to --suffix; give a mask variant its own name")
    a = ap.parse_args()
    out_suffix = a.suffix if a.out_suffix is None else a.out_suffix
    mask_c3 = {int(c) for c in a.mask_c3.split(",")}
    os.makedirs(os.path.join(OUT, "tail"), exist_ok=True)

    pat = re.compile(r"spectra\.\d+" + re.escape(a.suffix) + r"\.pkl$")
    files = sorted(f for f in glob.glob(os.path.join(OUT, "sfs/chunks/spectra.*.pkl")) if pat.search(f))
    n_chunks = sum(1 for _ in open(os.path.join(ROOT, "input_data/11_Demographic_Modeling/chunks.tsv"))) - 1
    if len(files) != n_chunks:
        raise SystemExit(f"{len(files)} spectra files for {n_chunks} chunks; finish 03_spectra.py first")

    # Two variants per pair: "" = polarized by outgroup agreement (unfolded analysis) and "_alt" = every
    # SNP, ALT-oriented, invariant sites whether or not the outgroups are callable (folded, label-free
    # analysis; DESIGN.md Deviation 1)
    VARIANTS = (("", "spectra", "invariant", "delta_hist"),
                ("_alt", "spectra_alt", "invariant_all", "delta_hist_alt"))
    acc = {}
    blockstats = collections.defaultdict(lambda: np.zeros(2))
    groups_of = {}
    with open(os.path.join(ROOT, "input_data/11_Demographic_Modeling/pairs.tsv")) as fh:
        next(fh)
        for line in fh:
            pid, bp, tp, role = line.rstrip("\n").split("\t")
            groups_of[pid] = (bp, tp)

    for path in files:
        with open(path, "rb") as fh:
            p = pickle.load(fh)
        chrom = p["chrom"]
        for (blk, g), (pi_s, d_s) in p["blockstats"].items():
            blockstats[(chrom, blk, g)] += (pi_s, d_s)
        for var, spec_key, inv_key, dh_key in VARIANTS:
            for pid, sp in p[spec_key].items():
                n1, n2 = p["projections"][pid]
                A = acc.setdefault((pid, var), {
                    "n": (n1, n2), "data": np.zeros((n1 + 1, n2 + 1)), "L_mask": 0.0, "L_scan": 0.0,
                    "region_spec": collections.defaultdict(lambda n1=n1, n2=n2: np.zeros((n1 + 1, n2 + 1))),
                    "block_spec": collections.defaultdict(lambda n1=n1, n2=n2: np.zeros((n1 + 1, n2 + 1))),
                    "block_inv": collections.Counter(),
                    "dhist": collections.defaultdict(lambda nb=len(p["delta_bins"]) - 1: np.zeros(nb)),
                    "delta_bins": p["delta_bins"]})
                for (blk, c3, region), arr in sp.items():
                    A["L_scan"] += arr.sum()
                    A["region_spec"][region] += arr
                    if region == "none" and c3 in mask_c3:
                        A["data"] += arr
                        A["L_mask"] += arr.sum()
                        A["block_spec"][(chrom, blk)] += arr
                for (c3, region), h in p[dh_key][pid].items():
                    A["dhist"][region] += h
            for (pid, blk, c3, reg), n in p[inv_key].items():
                if (pid, var) not in acc:
                    continue
                A = acc[(pid, var)]
                A["L_scan"] += n
                if reg == "none" and c3 in mask_c3:
                    A["L_mask"] += n
                    A["block_inv"][(chrom, blk)] += n

    # local-Ne factor per block for each pair: diversity / outgroup divergence in the populations not in
    # the pair, relative to the same ratio genome-wide (DESIGN.md section 4)
    blocks_all = sorted({(c, b) for c, b, g in blockstats}, key=lambda k: (int(k[0][3:]), k[1]))

    def f_blocks_for(pid):
        bp, tp = groups_of[pid]
        members = {bp, tp}
        if tp == "wTP1":
            members |= {"wAPA", "wBEN"}
        if tp == "wTP2":
            members |= {"wDOV", "wDOG"}
        off = [g for g in WHOLE8 if g not in members]
        tot = sum(blockstats[(c, b, g)] for c, b in blocks_all for g in off)
        genome_ratio = tot[0] / tot[1]
        d_med = np.median([sum(blockstats[(c, b, g)][1] for g in off) for c, b in blocks_all])
        f = {}
        for c, b in blocks_all:
            s = sum(blockstats[(c, b, g)] for g in off)
            f[(c, b)] = 1.0 if s[1] < 0.25 * d_med else (s[0] / s[1]) / genome_ratio
        return f, off

    for pid in sorted({k[0] for k in acc}):
        f_blocks, off = f_blocks_for(pid)
        out = {"pair": pid, "groups": groups_of[pid], "mask_c3": sorted(mask_c3),
               "f_blocks": f_blocks, "f_off_groups": off}
        for var, _, _, _ in VARIANTS:
            A = acc[(pid, var)]
            n1, n2 = A["n"]
            keys = sorted(A["block_spec"], key=lambda k: (int(k[0][3:]), k[1]))
            out.update({"n": A["n"], f"data{var}": A["data"], f"L_mask{var}": A["L_mask"],
                        f"L_scan{var}": A["L_scan"], f"region_spec{var}": dict(A["region_spec"]),
                        f"block_keys{var}": keys,
                        f"block_spec{var}": np.array([A["block_spec"][k] for k in keys]),
                        f"block_inv{var}": np.array([A["block_inv"][k] for k in keys]),
                        f"dhist{var}": dict(A["dhist"]), "delta_bins": A["delta_bins"]})
            print(f"{pid}{var or ' (polarized)'}: n=({n1},{n2}) SNPs in mask {A['data'].sum():.0f} "
                  f"L_mask {A['L_mask']:.0f} L_scan {A['L_scan']:.0f} blocks {len(keys)}")
        data = dadi.Spectrum(out["data"], mask_corners=True, data_folded=False, pop_ids=list(groups_of[pid]))
        data.to_file(os.path.join(OUT, f"sfs/{pid}{out_suffix}.fs"))
        dadi.Spectrum(out["data_alt"], mask_corners=True, data_folded=False,
                      pop_ids=list(groups_of[pid])).fold().to_file(os.path.join(OUT, f"sfs/{pid}{out_suffix}.folded.fs"))
        with open(os.path.join(OUT, f"sfs/{pid}{out_suffix}.pkl"), "wb") as fh:
            pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)

    # index SNPs, raw frequencies from the site table
    rows = {}
    for path in glob.glob(os.path.join(OUT, "site_table/index_snps.*.tsv")):
        with open(path) as fh:
            header = fh.readline().lstrip("#").rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(header, line.rstrip("\n").split("\t")))
                rows[(f["chrom"], int(f["pos"]))] = f
    peaks = []
    with open(os.path.join(ROOT, "input_data/10_Introgression/peaks.bed")) as fh:
        for line in fh:
            f = line.split()
            peaks.append((int(f[3]), "chr" + f[0] if not f[0].startswith("chr") else f[0], int(f[5])))
    with open(os.path.join(OUT, f"tail/index_units{out_suffix}.tsv"), "w") as o:
        o.write("pair_id\tpeak\town\tchrom\tpos\taa\tAN_bp\td_bp\tAN_tp\td_tp\tf_bp\tf_tp\tdelta\tdelta_alt\n")
        for pid in sorted({k[0] for k in acc}):
            bp, tp = groups_of[pid]
            for pk, chrom, pos in sorted(peaks):
                f = rows.get((chrom, pos))
                if f is None:
                    o.write(f"{pid}\t{pk}\tNA\t{chrom}\t{pos}\tmissing\t" + "\t".join(["NA"] * 8) + "\n")
                    continue
                anb, acb, ant, act = (to_int(f[f"AN_{bp}"]), to_int(f[f"AC_{bp}"]),
                                      to_int(f[f"AN_{tp}"]), to_int(f[f"AC_{tp}"]))
                aa = f["aa"]
                db = acb if aa == "R" else anb - acb
                dt = act if aa == "R" else ant - act
                own = pk in OWN.get(re.match(r"[A-Z]\d", pid).group(0), [])
                # unpolarized under the agreement rule: no derived-allele delta, but the label-free
                # (ALT-oriented) delta still defines the folded analysis's tail
                if aa == ".":
                    fb = ft = float("nan")
                else:
                    fb, ft = db / anb if anb else float("nan"), dt / ant if ant else float("nan")
                d_alt = (acb / anb - act / ant) if anb and ant else float("nan")
                o.write(f"{pid}\t{pk}\t{own}\t{chrom}\t{pos}\t{aa}\t{anb}\t{db}\t{ant}\t{dt}\t{fb:.4f}\t{ft:.4f}\t"
                        f"{fb - ft:.4f}\t{d_alt:.4f}\n")

    # fixed and near-fixed differences, raw delta, mask and genome-wide (outside peaks)
    with open(os.path.join(OUT, f"fixed_diff_counts{out_suffix}.tsv"), "w") as o:
        o.write("pair_id\tscope\tdirection\tthreshold\tn_snps\n")
        for (pid, var), A in sorted(acc.items()):
            if var:
                continue
            edges = A["delta_bins"]
            for scope, regions in (("nonpeak", [r for r in A["dhist"] if not r.startswith("peak")]),
                                   ("all", list(A["dhist"]))):
                h = sum(A["dhist"][r] for r in regions)
                for thr in (0.95, 1.0):
                    # delta = +/-1 sits in the outermost bins; no other delta reaches them (AN <= 134)
                    lo = thr if thr < 1 else edges[-2]
                    up = h[edges[:-1] >= lo - 1e-9].sum()
                    down = h[edges[1:] <= -lo + 1e-9].sum()
                    o.write(f"{pid}\t{scope}\tBP_derived\t{thr}\t{up:.0f}\n")
                    o.write(f"{pid}\t{scope}\tTP_derived\t{thr}\t{down:.0f}\n")


if __name__ == "__main__":
    main()
