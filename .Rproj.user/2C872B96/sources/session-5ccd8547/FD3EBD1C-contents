# scripts/04_deg_mrna_limma.R
# Stage 04: mRNA DEG per dataset (limma) on GPL570
#
# Inputs:
#   - data_raw/GEO/<GSE>_GSEMatrix_list.rds
#   - meta/<GSE>_sample_sheet.csv  (final_group column must be correct)
#
# Outputs:
#   - data_processed/<GSE>_expr_gene.rds
#   - results/tables/<GSE>_limma_all.tsv
#   - results/tables/<GSE>_limma_DEG.tsv

source("scripts/01_params.R")

suppressPackageStartupMessages({
  library(Biobase)
  library(limma)
  library(matrixStats)
  library(AnnotationDbi)
  library(hgu133plus2.db)
})

# -------- helpers --------
needs_log2 <- function(x) {
  qs <- quantile(x, probs = c(0.99, 1), na.rm = TRUE)
  (qs[1] > 100) || (qs[2] > 1000)
}
safe_log2 <- function(x) { x[x <= 0] <- NA; log2(x) }

collapse_probe_to_gene <- function(expr_mat) {
  probes <- rownames(expr_mat)
  symbols <- AnnotationDbi::mapIds(
    hgu133plus2.db,
    keys = probes,
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  keep <- !is.na(symbols)
  expr <- expr_mat[keep, , drop = FALSE]
  sym  <- symbols[keep]
  
  iqr <- matrixStats::rowIQRs(expr, na.rm = TRUE)
  df <- data.frame(probe = rownames(expr), gene = sym, iqr = iqr)
  
  df <- df[order(df$gene, -df$iqr), ]
  best <- df[!duplicated(df$gene), ]
  
  out <- expr[best$probe, , drop = FALSE]
  rownames(out) <- best$gene
  out
}

run_limma <- function(expr_gene, group_vec) {
  group <- factor(group_vec, levels = c("Normal", "EOC"))
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  fit <- lmFit(expr_gene, design)
  cont <- makeContrasts(EOCvsN = EOC - Normal, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))
  
  tt <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
  tt$gene <- rownames(tt)
  tt
}

process_one_mrna <- function(gse_id) {
  message("\n==============================")
  message("Stage 04 (limma DEG): ", gse_id)
  
  # GEO data
  gse_list <- readRDS(file.path(DIR_GEO, paste0(gse_id, "_GSEMatrix_list.rds")))
  eset <- gse_list[[1]]
  
  # sample sheet
  sheet <- read.csv(file.path(DIR_META, paste0(gse_id, "_sample_sheet.csv")))
  sheet$final_group[sheet$final_group == ""] <- NA
  
  keep <- !is.na(sheet$final_group) & sheet$final_group %in% c("EOC","Normal")
  sheet2 <- sheet[keep, , drop = FALSE]
  
  # verify counts vs expected
  exp <- expected_counts[[gse_id]]
  obs_eoc <- sum(sheet2$final_group == "EOC")
  obs_nor <- sum(sheet2$final_group == "Normal")
  message("Observed: EOC=", obs_eoc, " Normal=", obs_nor,
          " | Expected: EOC=", exp["EOC"], " Normal=", exp["Normal"])
  
  if (!(obs_eoc == exp["EOC"] && obs_nor == exp["Normal"])) {
    stop("Counts mismatch for ", gse_id, ". Fix final_group in meta sample sheet.")
  }
  
  # expression
  expr0 <- exprs(eset)
  
  # reorder columns to sample sheet order
  if (!all(sheet2$sample_id %in% colnames(expr0))) {
    missing <- sheet2$sample_id[!sheet2$sample_id %in% colnames(expr0)]
    stop("Sample IDs missing in expression for ", gse_id, ": ", paste(missing, collapse = ", "))
  }
  expr0 <- expr0[, sheet2$sample_id, drop = FALSE]
  
  # log2 if needed
  if (needs_log2(expr0)) {
    message("Detected non-log2 intensity -> applying log2 transform.")
    expr0 <- safe_log2(expr0)
  } else {
    message("Expression already log2-like -> no log2 transform.")
  }
  
  # probe -> gene collapse
  expr_gene <- collapse_probe_to_gene(expr0)
  message("Gene-level matrix: ", nrow(expr_gene), " genes x ", ncol(expr_gene), " samples.")
  
  # save processed gene-level matrix
  saveRDS(expr_gene, file = file.path(DIR_PROC, paste0(gse_id, "_expr_gene.rds")))
  
  # limma
  tt <- run_limma(expr_gene, sheet2$final_group)
  
  # DEG flag
  tt$DEG_flag <- (tt$adj.P.Val < DEG_FDR) & (abs(tt$logFC) > DEG_abs_logFC)
  
  # save tables
  out_all <- file.path(DIR_TAB, paste0(gse_id, "_limma_all.tsv"))
  out_deg <- file.path(DIR_TAB, paste0(gse_id, "_limma_DEG.tsv"))
  
  write.table(tt, out_all, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(tt[tt$DEG_flag, ], out_deg, sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Saved: ", out_all)
  message("Saved: ", out_deg)
  message("N_DEG = ", sum(tt$DEG_flag, na.rm = TRUE))
  
  invisible(tt)
}

# -------- run for all mRNA datasets --------
for (id in GSE_mrna) {
  process_one_mrna(id)
}

message("\nStage 04 complete.")
