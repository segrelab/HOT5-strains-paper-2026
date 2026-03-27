# Figure2.R
# RSCU-based cluster dendrogram of all replicons from five HOT5 strains
# Osborne et al. 2026
#
# Strains: B08, C03, B10, E06, F03
# Replicons: chromosome + all extra-chromosomal replicons (chromids/plasmids)
# Tip labels: name | length (kbp) | GC% | replication origin | T4SS +/-

library(seqinr)
library(ape)

setwd("/usr2/postdoc/melosbor/HOT5-strains-paper-2026/Figure2-RSCU-dendrogram")

# ============================================================
# PATHS
# ============================================================

base_dir <- "/projectnb/npsegre/Rhodobacters_sequencing/Genomics_analysis/Final_annotated_assemblies"
btc <- function(strain, file) {
  file.path(base_dir, strain, "broken_to_contigs", file)
}

# ============================================================
# REPLICON TABLE
# strain, display_name, type, fasta_path, gff_path
# ============================================================

replicons <- list(
  # --- B08 ---
  list(strain="B08", short="chromosome", type="chromosome",
       fasta=btc("B08", "B08_genome.fasta"), gff=btc("B08", "B08_genome.gff")),
  list(strain="B08", short="chromid",    type="chromid",
       fasta=btc("B08", "B08_chromid.fasta"), gff=btc("B08", "B08_chromid.gff")),
  # --- C03 ---
  list(strain="C03", short="chromosome", type="chromosome",
       fasta=btc("C03", "C03_genome.fasta"), gff=btc("C03", "C03_genome.gff")),
  list(strain="C03", short="chromid",    type="chromid",
       fasta=btc("C03", "C03_chromid.fasta"), gff=btc("C03", "C03_chromid.gff")),
  # --- B10 ---
  list(strain="B10", short="chromosome", type="chromosome",
       fasta=btc("B10", "B10_genome.fasta"),   gff=btc("B10", "B10_genome.gff")),
  list(strain="B10", short="p53kbp",     type="plasmid",
       fasta=btc("B10", "B10_p53kbp.fasta"),   gff=btc("B10", "B10_p53kbp.gff")),
  list(strain="B10", short="p86kbp",     type="plasmid",
       fasta=btc("B10", "B10_p86kbp.fasta"),   gff=btc("B10", "B10_p86kbp.gff")),
  list(strain="B10", short="p90kbp",     type="plasmid",
       fasta=btc("B10", "B10_p90kbp.fasta"),   gff=btc("B10", "B10_p90kbp.gff")),
  list(strain="B10", short="p102kbp",    type="plasmid",
       fasta=btc("B10", "B10_p102kbp.fasta"),  gff=btc("B10", "B10_p102kbp.gff")),
  list(strain="B10", short="p111kbp",    type="plasmid",
       fasta=btc("B10", "B10_p111kbp.fasta"),  gff=btc("B10", "B10_p111kbp.gff")),
  list(strain="B10", short="p125kbp",    type="plasmid",
       fasta=btc("B10", "B10_p125kbp.fasta"),  gff=btc("B10", "B10_p125kbp.gff")),
  list(strain="B10", short="p143kbp",    type="plasmid",
       fasta=btc("B10", "B10_p143kbp.fasta"),  gff=btc("B10", "B10_p143kbp.gff")),
  list(strain="B10", short="p198kbp",    type="plasmid",
       fasta=btc("B10", "B10_p198kbp.fasta"),  gff=btc("B10", "B10_p198kbp.gff")),
  list(strain="B10", short="p293kbp",    type="plasmid",
       fasta=btc("B10", "B10_p293kbp.fasta"),  gff=btc("B10", "B10_p293kbp.gff")),
  list(strain="B10", short="p360kbp",    type="plasmid",
       fasta=btc("B10", "B10_p360kbp.fasta"),  gff=btc("B10", "B10_p360kbp.gff")),
  list(strain="B10", short="p373kbp",    type="plasmid",
       fasta=btc("B10", "B10_p373kbp.fasta"),  gff=btc("B10", "B10_p373kbp.gff")),
  # --- E06 ---
  list(strain="E06", short="chromosome", type="chromosome",
       fasta=btc("E06", "E06_genome.fasta"),   gff=btc("E06", "E06_genome.gff")),
  list(strain="E06", short="p63kbp",     type="plasmid",
       fasta=btc("E06", "E06_p63kbp.fasta"),   gff=btc("E06", "E06_p63kbp.gff")),
  list(strain="E06", short="p87kbp",     type="plasmid",
       fasta=btc("E06", "E06_p87kbp.fasta"),   gff=btc("E06", "E06_p87kbp.gff")),
  list(strain="E06", short="p98kbp",     type="plasmid",
       fasta=btc("E06", "E06_p98kbp.fasta"),   gff=btc("E06", "E06_p98kbp.gff")),
  list(strain="E06", short="p99kbp",     type="plasmid",
       fasta=btc("E06", "E06_p99kbp.fasta"),   gff=btc("E06", "E06_p99kbp.gff")),
  list(strain="E06", short="p111kbp",    type="plasmid",
       fasta=btc("E06", "E06_p111kbp.fasta"),  gff=btc("E06", "E06_p111kbp.gff")),
  list(strain="E06", short="p143kbp",    type="plasmid",
       fasta=btc("E06", "E06_p143kbp.fasta"),  gff=btc("E06", "E06_p143kbp.gff")),
  list(strain="E06", short="p183kbp",    type="plasmid",
       fasta=btc("E06", "E06_p183kbp.fasta"),  gff=btc("E06", "E06_p183kbp.gff")),
  list(strain="E06", short="p188kbp",    type="plasmid",
       fasta=btc("E06", "E06_p188kbp.fasta"),  gff=btc("E06", "E06_p188kbp.gff")),
  list(strain="E06", short="p195kbp",    type="plasmid",
       fasta=btc("E06", "E06_p195kbp.fasta"),  gff=btc("E06", "E06_p195kbp.gff")),
  list(strain="E06", short="p303kbp",    type="plasmid",
       fasta=btc("E06", "E06_p303kbp.fasta"),  gff=btc("E06", "E06_p303kbp.gff")),
  list(strain="E06", short="p360kbp",    type="plasmid",
       fasta=btc("E06", "E06_p360kbp.fasta"),  gff=btc("E06", "E06_p360kbp.gff")),
  # --- F03 ---
  list(strain="F03", short="chromosome", type="chromosome",
       fasta=btc("F03", "F03_genome.fasta"),   gff=btc("F03", "F03_genome.gff")),
  list(strain="F03", short="p53kbp",     type="plasmid",
       fasta=btc("F03", "F03_p53kbp.fasta"),   gff=btc("F03", "F03_p53kbp.gff")),
  list(strain="F03", short="p63kbp",     type="plasmid",
       fasta=btc("F03", "F03_p63kbp.fasta"),   gff=btc("F03", "F03_p63kbp.gff")),
  list(strain="F03", short="p87kbp",     type="plasmid",
       fasta=btc("F03", "F03_p87kbp.fasta"),   gff=btc("F03", "F03_p87kbp.gff")),
  list(strain="F03", short="p98kbp",     type="plasmid",
       fasta=btc("F03", "F03_p98kbp.fasta"),   gff=btc("F03", "F03_p98kbp.gff")),
  list(strain="F03", short="p99kbp",     type="plasmid",
       fasta=btc("F03", "F03_p99kbp.fasta"),   gff=btc("F03", "F03_p99kbp.gff")),
  list(strain="F03", short="p111kbp",    type="plasmid",
       fasta=btc("F03", "F03_p111kbp.fasta"),  gff=btc("F03", "F03_p111kbp.gff")),
  list(strain="F03", short="p117kbp",    type="plasmid",
       fasta=btc("F03", "F03_p117kbp.fasta"),  gff=btc("F03", "F03_p117kbp.gff")),
  list(strain="F03", short="p143kbp",    type="plasmid",
       fasta=btc("F03", "F03_p143kbp.fasta"),  gff=btc("F03", "F03_p143kbp.gff")),
  list(strain="F03", short="p183kbp",    type="plasmid",
       fasta=btc("F03", "F03_p183kbp.fasta"),  gff=btc("F03", "F03_p183kbp.gff")),
  list(strain="F03", short="p188kbp",    type="plasmid",
       fasta=btc("F03", "F03_p188kbp.fasta"),  gff=btc("F03", "F03_p188kbp.gff")),
  list(strain="F03", short="p195kbp",    type="plasmid",
       fasta=btc("F03", "F03_p195kbp.fasta"),  gff=btc("F03", "F03_p195kbp.gff")),
  list(strain="F03", short="p303kbp",    type="plasmid",
       fasta=btc("F03", "F03_p303kbp.fasta"),  gff=btc("F03", "F03_p303kbp.gff")),
  list(strain="F03", short="p360kbp",    type="plasmid",
       fasta=btc("F03", "F03_p360kbp.fasta"),  gff=btc("F03", "F03_p360kbp.gff"))
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Parse Geneious GFF3, return data.frame with CDS rows
read_gff_cds <- function(gff_file) {
  lines <- readLines(gff_file, warn=FALSE)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(lines) == 0) return(data.frame())

  parsed <- lapply(lines, function(l) {
    x <- strsplit(l, "\t")[[1]]
    if (length(x) < 9) return(NULL)
    list(seqname=x[1], feature=x[3], start=as.integer(x[4]),
         end=as.integer(x[5]), strand=x[7], attrs=x[9])
  })
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) return(data.frame())

  df <- data.frame(
    seqname = sapply(parsed, `[[`, "seqname"),
    feature = sapply(parsed, `[[`, "feature"),
    start   = sapply(parsed, `[[`, "start"),
    end     = sapply(parsed, `[[`, "end"),
    strand  = sapply(parsed, `[[`, "strand"),
    attrs   = sapply(parsed, `[[`, "attrs"),
    stringsAsFactors = FALSE
  )
  df[df$feature == "CDS", ]
}

# Compute reverse complement of a nucleotide character vector
revcomp <- function(seq) {
  rev(comp(seq, ambiguous=TRUE))
}

# All 64 sense codons in alphabetical order (a/c/g/t).
# Pre-computed once so count_codons() doesn't rebuild it on every call.
.CODONS_64 <- apply(
  expand.grid(c("a","c","g","t"), c("a","c","g","t"), c("a","c","g","t"),
              stringsAsFactors = FALSE),
  1, paste0, collapse = "")

# Count raw codons from a seqinr-style character vector (single lowercase bases).
# Returns a named integer vector of length 64.
count_codons <- function(dna_vec) {
  str    <- paste(dna_vec, collapse = "")
  n      <- nchar(str)
  trim   <- n - (n %% 3L)
  starts <- seq(1L, trim - 2L, by = 3L)
  codons <- substring(str, starts, starts + 2L)
  tbl    <- table(factor(codons, levels = .CODONS_64))
  setNames(as.integer(tbl), .CODONS_64)
}

# Compute RSCU from a named vector of raw codon counts.
# Excludes stop codons (TAA, TAG, TGA).
compute_rscu <- function(raw) {
  stopcodons <- c("taa", "tag", "tga")
  raw_sense <- raw[!names(raw) %in% stopcodons]

  # Map each codon to its amino acid using the standard genetic code
  aa_map <- sapply(names(raw_sense), function(codon) {
    tryCatch(translate(s2c(codon)), error=function(e) NA_character_)
  })

  rscu_vals <- numeric(length(raw_sense))
  names(rscu_vals) <- names(raw_sense)

  for (aa in unique(na.omit(aa_map))) {
    syn <- names(aa_map)[!is.na(aa_map) & aa_map == aa]
    n   <- length(syn)
    tot <- sum(raw_sense[syn])
    for (codon in syn) {
      rscu_vals[codon] <- if (tot == 0) 1.0 else raw_sense[codon] / (tot / n)
    }
  }
  rscu_vals
}

# Extract all CDS sequences from a fasta+gff pair, return summed raw codon counts
replicon_raw_counts <- function(fasta_file, gff_file) {
  fa   <- read.fasta(fasta_file, seqtype="DNA", forceDNAtolower=TRUE, whole.header=FALSE)
  seq  <- fa[[1]]
  seqlen <- length(seq)

  cds <- read_gff_cds(gff_file)
  if (nrow(cds) == 0) return(NULL)

  raw_total <- setNames(rep(0L, 64L), .CODONS_64)

  for (i in seq_len(nrow(cds))) {
    s <- cds$start[i]
    e <- cds$end[i]
    if (is.na(s) || is.na(e) || s < 1 || e > seqlen || e < s) next

    subseq <- seq[s:e]
    if (cds$strand[i] == "-") subseq <- revcomp(subseq)

    cds_len <- length(subseq)
    # Trim to complete codons
    trim <- cds_len - (cds_len %% 3L)
    if (trim < 3L) next
    subseq <- subseq[seq_len(trim)]

    raw_total <- raw_total + count_codons(subseq)
  }
  raw_total
}

# Extract metadata: length (kbp), GC%, replication origin, T4SS presence
replicon_meta <- function(fasta_file, gff_file) {
  fa  <- read.fasta(fasta_file, seqtype="DNA", forceDNAtolower=TRUE, whole.header=FALSE)
  seq <- fa[[1]]

  len_kbp <- round(length(seq) / 1e3, 1)
  gc_pct  <- round(GC(seq) * 100, 1)

  # Read raw GFF text for keyword searches
  raw_gff  <- tolower(paste(readLines(gff_file, warn=FALSE), collapse=" "))

  # Replication origin classification
  has_dnaa <- grepl("dnaa|chromosomal replication initiator", raw_gff)
  has_repc <- grepl("repc cds|replication initiation protein repc|replication initiator repc", raw_gff)
  has_repa <- grepl("name=repa cds|replication protein a|replication initiator protein repa", raw_gff) ||
              grepl("\\brepa\\b", raw_gff)  # fallback: any standalone "repA"

  rep_origin <- if (has_dnaa && !has_repc && !has_repa) {
    "oriC"
  } else if (has_dnaa && (has_repa || has_repc)) {
    "oriC"   # chromosomes annotated with DnaA take priority
  } else if (has_repc) {
    "RepC"
  } else if (has_repa) {
    "RepA"
  } else {
    "Rep?"
  }

  # T4SS: VirB proteins or "type iv secretion" (not type IV pili)
  has_t4ss <- grepl("virb[2-9]|virb1[0-1]|vird4|type iv secretion system|t4ss", raw_gff)
  t4ss <- if (has_t4ss) "T4SS+" else "T4SS-"

  list(len_kbp=len_kbp, gc_pct=gc_pct, rep_origin=rep_origin, t4ss=t4ss)
}

# ============================================================
# COMPUTE RSCU AND METADATA FOR ALL REPLICONS
# ============================================================

message("Computing RSCU for all replicons...")
n_rep   <- length(replicons)
rscu_list <- vector("list", n_rep)
meta_list <- vector("list", n_rep)
ids       <- character(n_rep)

for (i in seq_len(n_rep)) {
  r    <- replicons[[i]]
  id   <- paste(r$strain, r$short, sep="_")
  ids[i] <- id
  message(sprintf("  [%d/%d] %s", i, n_rep, id))

  raw <- tryCatch(
    replicon_raw_counts(r$fasta, r$gff),
    error = function(e) { message("    ERROR: ", e$message); NULL }
  )

  meta <- tryCatch(
    replicon_meta(r$fasta, r$gff),
    error = function(e) { message("    ERROR: ", e$message); NULL }
  )

  rscu_list[[i]] <- if (!is.null(raw)) compute_rscu(raw) else NULL
  meta_list[[i]] <- meta
}

# Drop any replicon whose RSCU computation failed
ok <- !sapply(rscu_list, is.null) & !sapply(meta_list, is.null)
if (any(!ok)) {
  warning("Dropping replicons with failed computation: ", paste(ids[!ok], collapse=", "))
}
rscu_list  <- rscu_list[ok]
meta_list  <- meta_list[ok]
ids        <- ids[ok]
rep_types  <- sapply(replicons[ok], `[[`, "type")
rep_strains <- sapply(replicons[ok], `[[`, "strain")

# ============================================================
# BUILD RSCU MATRIX AND DISTANCE MATRIX
# ============================================================

# Ensure all RSCU vectors share the same codon names in the same order
codon_names <- names(rscu_list[[1]])
rscu_mat <- do.call(rbind, lapply(rscu_list, function(v) v[codon_names]))
rownames(rscu_mat) <- ids

# Euclidean distance on RSCU values (sense codons, stop codons excluded)
dist_mat <- dist(rscu_mat, method="euclidean")

# Ward's D2 hierarchical clustering
hc <- hclust(dist_mat, method="ward.D2")

# ============================================================
# BUILD TIP LABELS
# ============================================================

# Format: "B10 p53kbp  |  53.7 kbp  |  62.1% GC  |  RepA  |  T4SS-"
tip_labels <- mapply(function(id, meta, r) {
  display_name <- paste(r$strain, r$short)
  sprintf("%-20s  %7.1f kbp  |  %4.1f%% GC  |  %-5s  |  %s",
          display_name, meta$len_kbp, meta$gc_pct, meta$rep_origin, meta$t4ss)
}, ids, meta_list, replicons[ok], SIMPLIFY=TRUE)

names(tip_labels) <- ids

# ============================================================
# COLOR SCHEME
# ============================================================

strain_colors <- c(
  "B08" = "#E41A1C",  # red
  "C03" = "#377EB8",  # blue
  "B10" = "#4DAF4A",  # green
  "E06" = "#FF7F00",  # orange
  "F03" = "#984EA3"   # purple
)

type_lty <- c(
  "chromosome" = 1,   # solid
  "chromid"    = 2,   # dashed
  "plasmid"    = 3    # dotted
)

tip_cols <- strain_colors[rep_strains]
tip_lty  <- type_lty[rep_types]
names(tip_cols) <- ids
names(tip_lty)  <- ids

# ============================================================
# PLOT DENDROGRAM
# ============================================================

phylo_tree <- as.phylo(hc)

# Map hclust tip order to phylo tip labels
# as.phylo preserves tip labels from hclust$labels
phylo_tree$tip.label <- tip_labels[phylo_tree$tip.label]
tip_cols_ordered  <- tip_cols[hc$labels[hc$order]]
tip_lty_ordered   <- tip_lty[hc$labels[hc$order]]

# Figure dimensions: tall enough for all tips
n_tips  <- length(phylo_tree$tip.label)
fig_h   <- max(12, n_tips * 0.35)

pdf("Figure2_RSCU_dendrogram.pdf", width=18, height=fig_h)

par(mar=c(3, 1, 3, 1))

plot.phylo(
  phylo_tree,
  type      = "phylogram",
  direction = "rightwards",
  show.tip.label = TRUE,
  tip.color = tip_cols[names(tip_labels)[match(phylo_tree$tip.label, tip_labels)]],
  cex       = 0.75,
  font      = 1,
  label.offset = 0.01,
  no.margin = FALSE
)

title(
  main = "RSCU-based cluster dendrogram - HOT5 strains replicons",
  sub  = "Ward's D2 linkage on Euclidean distance of RSCU profiles",
  cex.main = 1.1, cex.sub = 0.85
)

# Legend: strains (color) and replicon types (line type)
legend(
  "bottomleft",
  legend = c(names(strain_colors), "", "chromosome", "chromid", "plasmid"),
  col    = c(strain_colors, NA, "black", "black", "black"),
  lty    = c(rep(1, length(strain_colors)), NA, 1, 2, 3),
  lwd    = c(rep(2, length(strain_colors)), NA, 1, 1, 1),
  pch    = c(rep(NA, length(strain_colors)+1), NA, NA, NA),
  bty    = "n",
  cex    = 0.8,
  title  = "Strain / Replicon type"
)

dev.off()

# ============================================================
# ALSO SAVE RSCU MATRIX AND METADATA AS CSV
# ============================================================

write.csv(rscu_mat, "rscu_matrix.csv")

meta_df <- data.frame(
  replicon   = ids,
  strain     = rep_strains,
  type       = rep_types,
  len_kbp    = sapply(meta_list, `[[`, "len_kbp"),
  gc_pct     = sapply(meta_list, `[[`, "gc_pct"),
  rep_origin = sapply(meta_list, `[[`, "rep_origin"),
  t4ss       = sapply(meta_list, `[[`, "t4ss"),
  stringsAsFactors = FALSE
)
write.csv(meta_df, "replicon_metadata.csv", row.names=FALSE)

message("Done. Output files:")
message("  Figure2_RSCU_dendrogram.pdf")
message("  rscu_matrix.csv")
message("  replicon_metadata.csv")

