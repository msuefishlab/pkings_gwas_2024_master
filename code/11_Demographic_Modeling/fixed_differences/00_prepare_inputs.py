#!/usr/bin/env python3
"""00_prepare_inputs.py
Write every input the site-table pass needs, so that one pass over the all-sites VCF serves the
headline test, every sensitivity row, the TP-TP controls, and the pixy reproduction check.

  groups.tsv   sample -> comma-separated groups, in the format of `bcftools +fill-tags -S`
  pairs.tsv    the population pairs (headline, single-TP, pooled BP2, TP-TP controls, pixy checks)
  chunks.tsv   chr1-25 in pieces of at most 10 Mb (one SLURM array task each)
  exons.bed    GFF exons (all transcript types), merged; the site table stores distance to them
  exclude.bed  regions never used for fitting: GWAS extended peaks +/-250 kb and the chr4
               inversion-like polymorphism (chr4:40.88-42.77 Mb, 14_Wobble_Analysis)

Groups. Two samplings of each creek, both computed so that DESIGN.md can name either as the
headline without a second pass over the VCF:
  p<POP>  phenotype-pure GWAS lists (input_data/12_Sweep_Detection/split_up/), split by name prefix
  w<POP>  whole populations as pixy saw them (input_data/05_PopGen/popsfile.txt), wobble fish included
plus the pooled references pTP1/pTP2/pBP2 and wTP1/wTP2/wBP2, and one group per outer-clade outgroup
species (OUT_PSZA, OUT_PMAG, OUT_PPAR, OUT_PHOP, OUT_PGAB, OUT_PVAD) from the popsfile. PSZA_455 is
P. hopkinsi (fish_data_2023) and the popsfile lists it under PHOP.

Every listed sample is checked against the VCF header (--vcf_samples, from `bcftools query -l`);
any mismatch aborts.

Usage (from the repo root, inside gwas_tools.sif):
  bcftools query -l "$VCF" > /tmp/vcf_samples.txt
  python3 code/11_Demographic_Modeling/fixed_differences/00_prepare_inputs.py --vcf_samples /tmp/vcf_samples.txt
"""
import argparse
import collections
import os
import sys

POPS = ["BAM", "APA", "BEN", "BIK", "BAVA", "BIR", "DOV", "DOG"]
SPLIT_FILES = ["BP1", "BP2", "BP3", "TP1", "TP2"]
# Outer clade of the Paramormyrops species tree (output_data/04_Phylogenomics/...cf.tree): the species
# used to call the ancestral allele by agreement (DESIGN.md, Deviation 1). The curvifrons group
# (PCUR, PNTO, PLON, PSN2) is closer to P. kingsleyae and shares excess derived alleles with it (D-006),
# so it is not used.
OUTGROUP_SPECIES = ["PSZA", "PMAG", "PPAR", "PHOP", "PGAB", "PVAD"]
POOLED = {"TP1": ["APA", "BEN"], "TP2": ["DOV", "DOG"], "BP2": ["BIK", "BAVA"]}
CHUNK = 10_000_000
PEAK_BUFFER = 250_000
CHR4_INVERSION = ("chr4", 40_880_000, 42_770_000)

# pair_id, bp group, tp group, role. Roles: headline (both samplings), single_tp, pooled_bp2,
# tp_control (no BP selection expected: calibrates the model's tail), check (pixy dxy only).
def build_pairs():
    pairs = []
    for s in ("p", "w"):
        tag = "pure" if s == "p" else "whole"
        for pid, bp, tp in (("A1", "BAM", "TP1"), ("A2", "BIK", "TP2"),
                            ("A3", "BAVA", "TP2"), ("A4", "BIR", "TP2")):
            pairs.append((f"{pid}{s}", s + bp, s + tp, f"headline_{tag}"))
            for single in POOLED[tp]:
                pairs.append((f"{pid}{s}_{single}", s + bp, s + single, f"single_tp_{tag}"))
        pairs.append((f"B2{s}", s + "BP2", s + "TP2", f"pooled_bp2_{tag}"))
        pairs.append((f"C1{s}", s + "APA", s + "BEN", f"tp_control_{tag}"))
        pairs.append((f"C2{s}", s + "DOV", s + "DOG", f"tp_control_{tag}"))
    # pixy dxy reproduction (whole populations, as pixy ran)
    for a, b in (("APA", "BAM"), ("BAM", "BEN"), ("BIK", "BAVA"), ("BIK", "DOV"),
                 ("BAVA", "DOV"), ("BIR", "DOV")):
        pairs.append((f"X_{a}_{b}", "w" + a, "w" + b, "check"))
    return pairs


def read_list(path):
    with open(path) as fh:
        return [l.strip() for l in fh if l.strip()]


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf_samples", required=True, help="output of bcftools query -l on the all-sites VCF")
    ap.add_argument("--root", default=os.getcwd())
    args = ap.parse_args()
    root = args.root
    out_in = os.path.join(root, "input_data/11_Demographic_Modeling")
    out_mask = os.path.join(root, os.environ.get("FD_OUT", "output_data/11_Demographic_Modeling/fixed_differences"), "masks")
    os.makedirs(out_in, exist_ok=True)
    os.makedirs(out_mask, exist_ok=True)

    vcf_samples = set(read_list(args.vcf_samples))
    groups = collections.defaultdict(list)

    # phenotype-pure GWAS lists, split by population prefix
    for f in SPLIT_FILES:
        for s in read_list(os.path.join(root, "input_data/12_Sweep_Detection/split_up", f + ".txt")):
            pop = s.split("_")[0]
            if pop not in POPS:
                sys.exit(f"unexpected population prefix {pop} in {f}.txt")
            groups[s].append("p" + pop)
    # whole populations and the outgroup, from the pixy popsfile
    with open(os.path.join(root, "input_data/05_PopGen/popsfile.txt")) as fh:
        for line in fh:
            if not line.strip():
                continue
            s, pop = line.rstrip("\n").split("\t")[:2]
            if pop in POPS:
                groups[s].append("w" + pop)
            elif pop in OUTGROUP_SPECIES:
                groups[s].append("OUT_" + pop)
    if "OUT_PSZA" in groups.get("PSZA_455", []):
        sys.exit("PSZA_455 is P. hopkinsi; the popsfile should list it as PHOP")
    # pooled references
    for s, gs in groups.items():
        for pooled, members in POOLED.items():
            for prefix in ("p", "w"):
                if any(g == prefix + m for m in members for g in gs):
                    gs.append(prefix + pooled)

    missing = sorted(s for s in groups if s not in vcf_samples)
    if missing:
        sys.exit(f"{len(missing)} listed samples are not in the VCF: {missing[:10]}")

    with open(os.path.join(out_in, "groups.tsv"), "w") as fh:
        for s in sorted(groups):
            fh.write(f"{s}\t{','.join(groups[s])}\n")
    counts = collections.Counter(g for gs in groups.values() for g in gs)
    with open(os.path.join(out_in, "group_sizes.tsv"), "w") as fh:
        fh.write("group\tn_samples\n")
        for g in sorted(counts):
            fh.write(f"{g}\t{counts[g]}\n")

    with open(os.path.join(out_in, "pairs.tsv"), "w") as fh:
        fh.write("pair_id\tbp\ttp\trole\n")
        for p in build_pairs():
            for g in p[1:3]:
                if g not in counts:
                    sys.exit(f"pair {p[0]} uses unknown group {g}")
            fh.write("\t".join(p) + "\n")

    # chunks over chr1-25
    fai = os.path.join(root, "input_data/00_Reference_Genome",
                       "jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai")
    chroms = []
    with open(fai) as fh:
        for line in fh:
            c, length = line.split("\t")[:2]
            if c.startswith("chr") and c[3:].isdigit() and 1 <= int(c[3:]) <= 25:
                chroms.append((c, int(length)))
    with open(os.path.join(out_in, "chunks.tsv"), "w") as fh:
        fh.write("chunk_id\tchrom\tstart\tend\n")
        i = 0
        for c, length in sorted(chroms, key=lambda x: int(x[0][3:])):
            for start in range(1, length + 1, CHUNK):
                end = min(start + CHUNK - 1, length)
                fh.write(f"{i}\t{c}\t{start}\t{end}\n")
                i += 1

    # exons (all transcript types), merged, 0-based half-open BED
    ex = collections.defaultdict(list)
    with open(os.path.join(root, "input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 4 and f[2] == "exon":
                ex[f[0]].append((int(f[3]) - 1, int(f[4])))
    with open(os.path.join(out_mask, "exons.bed"), "w") as fh:
        for c, _ in sorted(chroms, key=lambda x: int(x[0][3:])):
            for s, e in merge(ex.get(c, [])):
                fh.write(f"{c}\t{s}\t{e}\n")

    # exclusions: GWAS extended peaks (D-002) +/-250 kb, and the chr4 inversion-like polymorphism
    excl = collections.defaultdict(list)
    with open(os.path.join(root, "input_data/10_Introgression/peaks.bed")) as fh:
        for line in fh:
            f = line.split()
            c = "chr" + f[0] if not f[0].startswith("chr") else f[0]
            excl[c].append((max(0, int(f[1]) - PEAK_BUFFER), int(f[2]) + PEAK_BUFFER, "peak" + f[3]))
    excl[CHR4_INVERSION[0]].append((CHR4_INVERSION[1], CHR4_INVERSION[2], "chr4_inversion"))
    with open(os.path.join(out_mask, "exclude.bed"), "w") as fh:
        for c in sorted(excl, key=lambda x: int(x[3:])):
            for s, e, name in sorted(excl[c]):
                fh.write(f"{c}\t{s}\t{e}\t{name}\n")

    print("groups:", ", ".join(f"{g}={counts[g]}" for g in sorted(counts)))
    print("pairs:", len(build_pairs()), "| chunks:", i, "| exon intervals:", sum(len(merge(v)) for v in ex.values()))


if __name__ == "__main__":
    main()
