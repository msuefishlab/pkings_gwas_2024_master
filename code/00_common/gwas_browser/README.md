# GWAS browser

A local, browser-based viewer for comparing any number of GEMMA results side by side against the genome and
annotation. It replaces loading `.GWAS` files into IGV, which parses every SNP of every file into memory.

`prepare_tracks.py` converts each result once into compact binary files (a thinned genome-wide layer plus
full-resolution per-contig files fetched only when zoomed in). The viewer is a static page (`index.html`,
`viewer.js`, `core.js`) with no build step and no dependencies, served from the repo root on localhost.

## Usage

```bash
source pkings_gwas.env
cd $root

# 1. Write a track manifest (or copy tracks.example.json and edit it)
python3 code/00_common/gwas_browser/prepare_tracks.py --discover > code/00_common/gwas_browser/tracks.json

# 2. Prepare data (skips tracks whose source and options are unchanged)
python3 code/00_common/gwas_browser/prepare_tracks.py code/00_common/gwas_browser/tracks.json

# 3. Serve on localhost only and open the printed URL
bash code/00_common/gwas_browser/serve.sh
```

The server binds to 127.0.0.1, so the unpublished data is not reachable from the network.

Prepared data goes to `output_data/gwas_browser/`, which carries its own `.gitignore` (`output_data/` itself is
not ignored).

The server sends `Cache-Control: no-cache`, so after re-running the prep a page reload picks up the new data.

## Viewer controls

- Wheel: zoom around the cursor. Drag: pan. Shift-drag: zoom to the selected region. Click the genome bar to jump
  to a chromosome, or drag on it to select a range.
- Keys: `+`/`-` zoom, arrow keys pan, `0` whole genome, `/` search, `Esc` clear the selection.
- Locus box: a gene name, a description word, `chr16`, or `chr16:10,800,000-11,100,000`.
- Click a SNP: rs, p, af, beta/se for every GWAS track at that position, the nearest gene (exon, intron, or distance
  5′/3′), and line-track values there. The position is ringed on every track.
- Click a gene: its description and, for every GWAS track, the strongest SNP within the gene ±10 kb.
- Y scale: per track (genome max), per track in view, or shared; `max` fixes it. Line tracks without a preset
  range scale to the data in view.
- Above 100 Mb the GWAS tracks show the thinned layer: SNPs with −log10 p ≥ 2 as dots and the per-pixel maximum of
  the rest as a shaded mass. Below 100 Mb every SNP is drawn.
- Genes: density above 10 Mb, models below 10 Mb, labels below 2 Mb. Ups scaffolds are off the axis until toggled
  (or until you navigate to one).
- The URL hash holds the locus, selected SNP, hidden tracks, y-scale and height, so a view can be bookmarked.
- `window.gwasBrowser` in the browser console exposes `state`, `setView`, `goText`, `selectSnp`, `setTheme` and
  `drawAll` for debugging and timing.

## Source manifest (`tracks.json`)

Paths are relative to the repo root.

```json
{"genome": {"fai": "input_data/00_Reference_Genome/<assembly>.renamed.fa.fai",
            "gff": "input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.sorted.gff"},
 "tracks": [
   {"id": "wobexcl", "label": "All, wobble excluded (LRT)", "type": "gwas",
    "path": "output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt",
    "p_col": "auto", "min_txt": "auto", "peaks": "auto", "peaks_unfiltered": "auto",
    "min_af": 0, "sig": null, "sug": null},
   {"id": "fst_apa_bam", "type": "line", "preset": "pixy_fst",
    "path": "output_data/05_PopGen/<prefix>_fst.txt", "filter": {"pop1": "APA", "pop2": "BAM"}}]}
```

GWAS track options:

- `path`: a GEMMA `*.assoc.txt` (lmm1 or lmm2) or a `.GWAS` file (`CHR BP SNP P`). Plain or `.gz`.
- `p_col`: `auto` takes the first of `p_wald`, `p_lrt`, `p_score`, `P` (case-insensitive).
- `min_txt`: `auto` means `<dir>/<NAME>_permution/min.txt`, where NAME is the file name up to the first dot.
  Thresholds follow `code/06_Association/gwas_report.Rmd:35-39`: with the permutation min p-values sorted
  ascending and n of them, `sig = -log10(p[n - floor(0.95 n)])` and `sug = -log10(p[n - floor(0.67 n)])`
  (1-based). p is column 12 for lmm1 (`p_wald`) and column 10 for lmm2 (`p_lrt`). Numbers in `sig`/`sug`
  override the computed values; `null` for `min_txt` means no threshold lines.
- `peaks`, `peaks_unfiltered`: `auto` means `<dir>/<NAME>_PEAKS.bed` and `<dir>/<NAME>_PEAKS_unfiltered.bed`.
  Their numeric chromosome N means `chrN` (ups contigs are dropped before peak calling, `gwas_report.Rmd:48`).
  Coordinates are 1-based inclusive SNP positions.
- `min_af`: drop SNPs with `af` below this value. Peaks were called on `af > 0.05`, so use 0.05 to match the
  report.

Line track presets:

| preset  | source                                          | value         | position                                |
|---------|-------------------------------------------------|---------------|-----------------------------------------|
| pixy_pi  | `05_PopGen/*_pi.txt`, `filter.pop`             | `avg_pi`      | `window_pos_1`..`window_pos_2`          |
| pixy_fst | `05_PopGen/*_fst.txt`, `filter.pop1/pop2` (unordered) | `avg_wc_fst` | `window_pos_1`..`window_pos_2`   |
| pixy_dxy | `05_PopGen/*_dxy.txt`, `filter.pop1/pop2` (unordered) | `avg_dxy`    | `window_pos_1`..`window_pos_2`   |
| tajd     | `05_PopGen/*_tajd.txt`, `filter.pop1`          | `TajimaD`     | `BIN_START+1`..`BIN_START+window` (5000) |
| dsuite   | `10_Introgression/<trio>_localFstats_chr*_50_25.txt` (glob) | `value_col`, default `f_dM` | `windowStart`..`windowEnd` |
| sweed, sf2 | `12_Sweep_Detection/{sweed,sweepfinder2}/merged/<POP>.genome_wide_clr.txt` | `CLR` | `round(position)` (point) |

`NA`, `nan` and empty values become gaps. Optional keys on any line track: `label`, `ymin`, `ymax`
(fixed y-range), `zero_line` (default true for dsuite and tajd). Presets fix the range for `pixy_fst` (0 to 1, so
negative Weir-Cockerham estimates draw at 0) and `dsuite` (−1 to 1); other line tracks scale to the data in view.

## Prepared data contract (`output_data/gwas_browser/`)

All binary data is little-endian. Contig order everywhere is the display order in `manifest.json`:
`chr1`..`chr25` numerically, then `upsN` numerically, then any other names alphabetically. Element offsets
below count elements, not bytes. Positions are 1-based.

### `manifest.json`

```json
{
  "version": 1,
  "built": "2026-10-12T15:04:05",
  "logp_scale": 1000,
  "bin_size": 10000,
  "contigs": [{"name": "chr1", "length": 77100000}],
  "genes": {"index": "genes/index.tsv", "contigs": ["chr1"]},
  "tracks": [
    {"id": "wobexcl", "type": "gwas", "label": "...", "source": "output_data/...assoc.txt",
     "test": "lmm2", "p_col": "p_lrt", "n_snps": 5898721, "n_nan": 0, "n_capped": 0, "min_af": 0,
     "overview_cutoff": 2.0, "detail_cutoff": 2.0,
     "sig": 7.2223, "sug": 6.3572, "thresh_source": "output_data/.../min.txt",
     "detail_fields": ["rs", "p", "af", "beta", "se", "n_miss", "allele1", "allele0"],
     "peaks": [{"contig": "chr16", "start": 10861379, "end": 11021498, "id": "4",
                "max_logp": 9.733, "lead": 10916883}],
     "peaks_unfiltered": [{"contig": "chr16", "start": 10861379, "end": 11021498}],
     "max_logp": 9.733,
     "contigs": {"chr16": {"n": 210000, "ov_off": 0, "ov_n": 1200, "bin_off": 0, "bin_n": 3400,
                           "max_logp": 9.733}}},
    {"id": "fst_apa_bam", "type": "line", "label": "...", "preset": "pixy_fst", "source": "...",
     "value_name": "avg_wc_fst", "kind": "window", "ymin": 0, "ymax": 1, "zero_line": false,
     "vmin": -0.1, "vmax": 0.98, "n": 187002,
     "contigs": {"chr16": {"off": 0, "n": 4000}}}
  ]
}
```

- `sig`, `sug`, `thresh_source` may be `null`. `peaks`, `peaks_unfiltered` may be empty lists.
- A GWAS track's `contigs` holds only contigs with SNPs. `bin_n = ceil(length / bin_size)` for that contig.
- `detail_fields` lists the columns present in the detail files; `.GWAS` sources have only `rs` and `p`.
- Line `kind` is `window` (start..end) or `point` (start == end). `ymin`/`ymax` are `null` unless fixed.
- `genes` is `null` when no GFF was given.

### GWAS track files (`gwas/<id>/`)

- `overview.bin`: SNPs with logp >= `overview_cutoff`, all contigs concatenated. Layout `[pos u32 x N][logp u16 x N]`,
  N = sum of `ov_n`; contig slice starts at `ov_off`. Sorted by position within each contig.
  logp is stored as `round(logp * logp_scale)`, capped at 65535 (`n_capped` counts p = 0 or logp > 65.535).
- `bins.bin`: `u16 x sum(bin_n)`, the max scaled logp of SNPs in each `bin_size` bin (bin k covers positions
  `k*bin_size+1`..`(k+1)*bin_size`); 0 = empty. Contig slice starts at `bin_off`.
- `full/<contig>.bin`: every SNP on that contig, `[pos u32 x n][logp u16 x n]`, sorted by position.
- `detail/<contig>.json`: columnar detail for SNPs with logp >= `detail_cutoff`, sorted by position:
  `{"pos": [...], "rs": [...], "p": ["1.23e-09", ...], "af": [...], ...}` with one array per `detail_fields`
  entry plus `pos`. `p` keeps the original string from the source file. Contigs with no such SNP have no file.
- `meta.json`: source size, mtime and an options hash, used to skip unchanged tracks.

### Line track file (`line/<id>.bin`)

`[start u32 x N][end u32 x N][value f32 x N]`, N = track `n`, contig slice starts at `off`, sorted by start.
Missing values are NaN.

### Genes (`genes/`)

- `genes/<contig>.json`: array of gene records sorted by start:
  `{"s": 77168, "e": 100526, "st": "+", "n": "ikbip", "id": "gene-ikbip", "b": "protein_coding",
    "d": "IKBKB interacting protein", "t": "rna-XM_...", "x": [s1, e1, s2, e2, ...], "cs": 77300, "ce": 99000}`.
  `x` is the exon list of one representative transcript (longest CDS, else longest exon total; empty for genes
  without transcripts), `cs`/`ce` are the CDS bounds or `null`. Descriptions are %-decoded.
- `genes/index.tsv`: header `name contig start end strand biotype description`, one line per gene, tab-separated.
