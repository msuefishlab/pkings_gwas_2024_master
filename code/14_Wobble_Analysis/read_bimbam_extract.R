library(data.table)
library(stringr)

# final_long must have columns: IID, peak, chr, pos, dosage
# order_file is the same .order file you used for the BIMBAM
enforce_sample_order <- function(final_long, order_file,
                                 sort_keys = c("peak","chr","pos"),
                                 samples_first = FALSE) {
  x <- as.data.table(final_long)
  ord <- fread(order_file, header = FALSE)[, .(IID = trimws(V1))]
  ord[, sample_index := .I]
  
  # attach index; keep all rows from final_long
  x <- ord[x, on = "IID"]  # left join from ord onto x to bring sample_index
  # If some IIDs weren’t in the order file, put them at the end
  x[is.na(sample_index), sample_index := .Machine$integer.max]
  
  # optional: ensure deterministic ordering of locus keys
  if ("peak" %in% names(x)) x[, peak := suppressWarnings(as.integer(peak)) %>% fifelse(is.na(peak), peak, as.integer(peak))]
  
  if (samples_first) {
    setorder(x, sample_index, !!!lapply(sort_keys, as.name))
  } else {
    setorder(x, !!!lapply(sort_keys, as.name), sample_index)
  }
  
  x[, sample_index := NULL]
  x[]
}



library(data.table)

# df.bed must contain: maxp_pos (integer), peak (integer/char)
# If df.bed also has a chromosome column (e.g., 'chr'), we will join on (chr, pos).
# Otherwise we join on pos only.
make_dosage_long <- function(res, df.bed) {
  # res$genotypes: matrix [variants x samples]
  # res$meta: data.table(variant_id, chr, pos, geno_variant_id, A1, A2)
  
  # 1) Melt genotypes to long: variant_id, IID, dosage
  Gdt <- as.data.table(res$genotypes, keep.rownames = "variant_id")
  longG <- data.table::melt(
    Gdt,
    id.vars = "variant_id",
    variable.name = "IID",
    value.name = "dosage"
  )
  longG[, IID := as.character(IID)]
  
  # 2) Attach chr/pos from meta
  meta_small <- res$meta[, .(variant_id, chr, pos)]
  longG <- meta_small[longG, on = "variant_id"]
  
  # 3) Prepare df.bed (ensure types and column names)
  bed <- as.data.table(df.bed)
  # normalize column names
  if (!"maxp_pos" %in% names(bed)) stop("df.bed must have a 'maxp_pos' column.")
  setnames(bed, "maxp_pos", "pos")
  
  # If df.bed has a chromosome column, use it; else join on pos only.
  has_chr <- any(c("chr", "chrom", "chromosome") %in% names(bed))
  if (!"chr" %in% names(bed) && "chrom" %in% names(bed)) setnames(bed, "chrom", "chr")
  if (!"chr" %in% names(bed) && "chromosome" %in% names(bed)) setnames(bed, "chromosome", "chr")
  
  bed[, pos := as.integer(pos)]  # be strict
  if ("chr" %in% names(bed)) bed[, chr := as.character(chr)]
  
  # 4) Join to add 'peak'
  if ("chr" %in% names(bed)) {
    setkey(bed, chr, pos)
    setkey(longG, chr, pos)
    longG <- bed[longG, nomatch = 0L][, .(IID, peak, chr, pos, dosage)]
  } else {
    setkey(bed, pos)
    setkey(longG, pos)
    longG <- bed[longG, nomatch = 0L][, .(IID, peak, chr, pos, dosage)]
  }
  
  # 5) Optional: order nicely
  setorder(longG, peak, chr, pos, IID)
  longG[]
}


# Read selected BIMBAM genotypes using data.table only
read_bimbam_extract <- function(
    geno_file,
    map_file,
    order_file,
    positions,          # integer vector (base-pair positions)
    chroms = NULL       # optional character vector (e.g. "chr16")
) {
  # ---- 1) MAP: parse positions and chromosome; record row indices ----
  # Example line: "chr16_232504_A_T:chr16:232504, 232504, chr16"
  map <- fread(
    map_file,
    header = FALSE,
    sep = ",",
    strip.white = TRUE,
    col.names = c("variant_id", "pos", "chr")
  )
  map[, pos := as.integer(pos)]
  map[, chr := as.character(chr)]
  map[, row_idx := .I]  # 1-based row index (matches .geno row order)
  
  # Optional chr filter + required pos filter
  if (!is.null(chroms)) {
    keep_map <- map[pos %in% positions & chr %in% chroms]
  } else {
    keep_map <- map[pos %in% positions]
  }
  if (nrow(keep_map) == 0L) stop("No requested positions found in the map.")
  
  # ---- 2) SAMPLE ORDER ----
  samples <- trimws(readLines(order_file))
  n_samp  <- length(samples)
  
  # ---- 3) GENO: read once, set column names, then subset by row_idx ----
  # BIMBAM .geno: 1=id, 2=A1, 3=A2, then one dosage per sample
  geno <- fread(geno_file, header = FALSE, na.strings = c("NA"))
  # Safety: check width
  expected_cols <- 3 + n_samp
  if (ncol(geno) < expected_cols) {
    stop(sprintf("Expected at least %d columns (3 + %d samples), got %d.",
                 expected_cols, n_samp, ncol(geno)))
  }
  setnames(geno, c("variant_id", "A1", "A2", samples))
  
  # Keep only the rows we want (use row index from the map to be order-safe)
  geno_keep <- geno[keep_map$row_idx]
  
  # ---- 4) Build outputs ----
  # Minimal metadata from map (align to the kept order)
  meta <- keep_map[, .(variant_id, chr, pos)]
  # Ensure variant_id consistency (qctool keeps the same row order/id)
  # but keep both, in case you want to check:
  meta[, geno_variant_id := geno_keep$variant_id]
  meta[, A1 := geno_keep$A1]
  meta[, A2 := geno_keep$A2]
  
  # Genotype matrix: variants x samples
  G <- as.matrix(geno_keep[, ..samples])
  storage.mode(G) <- "numeric"
  rownames(G) <- meta$variant_id
  colnames(G) <- samples
  
  list(
    meta = meta,             # variant_id, chr, pos, geno_variant_id, A1, A2
    genotypes = G,           # numeric matrix [variants x samples]
    table = cbind(meta, as.data.table(G))  # convenient long table
  )
}
