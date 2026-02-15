# scripts/03_make_sample_sheets.R
# Purpose:
#   Create or validate sample sheets for GEO datasets.
#   - If sample sheet already exists: DO NOT overwrite; just validate counts.
#   - If not exists (or force_recreate=TRUE): create from GEO metadata, write, then validate.
#
# Run:
#   source("scripts/03_make_sample_sheets.R")

source("scripts/01_params.R")

suppressPackageStartupMessages({
  library(Biobase)
})

# ---- user control ----
force_recreate <- FALSE  # IMPORTANT: keep FALSE so your manual edits won't be overwritten

# ---------- helpers ----------
collapse_chr <- function(pheno_row) paste(na.omit(as.character(pheno_row)), collapse = " | ")

infer_group_from_pheno <- function(pheno) {
  candidate_cols <- c(
    "title","source_name_ch1",
    "characteristics_ch1","characteristics_ch1.1","characteristics_ch1.2",
    "characteristics_ch1.3","characteristics_ch1.4","characteristics_ch1.5"
  )
  cols_try <- intersect(candidate_cols, colnames(pheno))
  
  txt <- apply(pheno[, cols_try, drop = FALSE], 1, collapse_chr)
  x <- tolower(txt)
  
  tumor  <- grepl("tumou?r|cancer|carcinoma|serous|hgsc|high[- ]grade|eoc|malignan|ovarian carcinoma", x)
  normal <- grepl("normal|control|healthy|non[- ]tumou?r|benign|noncancer", x)
  
  grp <- ifelse(tumor & !normal, "EOC",
                ifelse(normal & !tumor, "Normal", NA))
  grp
}

normalize_sheet <- function(sheet) {
  # convert blank strings to NA (critical for proper counting)
  if ("final_group" %in% colnames(sheet)) sheet$final_group[sheet$final_group == ""] <- NA
  if ("inferred_group" %in% colnames(sheet)) sheet$inferred_group[sheet$inferred_group == ""] <- NA
  sheet
}

validate_counts_or_stop <- function(sheet, gse_id) {
  sheet <- normalize_sheet(sheet)
  
  exp <- expected_counts[[gse_id]]
  if (is.null(exp)) {
    warning("No expected_counts defined for ", gse_id, ". Skipping strict validation.")
    return(invisible(TRUE))
  }
  
  obs_eoc <- sum(sheet$final_group == "EOC", na.rm = TRUE)
  obs_nor <- sum(sheet$final_group == "Normal", na.rm = TRUE)
  obs_na  <- sum(is.na(sheet$final_group))
  
  message("\n[GSE] ", gse_id,
          "\nObserved: EOC=", obs_eoc, " Normal=", obs_nor,
          "\nExpected: EOC=", exp["EOC"], " Normal=", exp["Normal"],
          "\nNA/Excluded: ", obs_na)
  
  ok <- (obs_eoc == exp["EOC"]) && (obs_nor == exp["Normal"])
  if (!ok) {
    stop(
      "\nSample counts do NOT match the manuscript for ", gse_id,
      "\nAction required:",
      "\n1) Open: meta/", gse_id, "_sample_sheet.csv",
      "\n2) Edit ONLY 'final_group' to exactly match expected counts.",
      "\n   Use only: EOC or Normal; leave non-used samples blank.",
      "\n3) Save and re-run script 03."
    )
  }
  invisible(TRUE)
}

create_sheet_from_geo <- function(gse_id) {
  rds_path <- file.path(DIR_GEO, paste0(gse_id, "_GSEMatrix_list.rds"))
  if (!file.exists(rds_path)) stop("Missing GEO RDS: ", rds_path)
  
  gse_list <- readRDS(rds_path)
  eset1 <- gse_list[[1]]
  pheno <- pData(eset1)
  
  inferred <- infer_group_from_pheno(pheno)
  
  # optional characteristics column for manual review
  char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  char_text <- if (length(char_cols) > 0) {
    apply(pheno[, char_cols, drop = FALSE], 1, collapse_chr)
  } else {
    rep(NA, nrow(pheno))
  }
  
  sheet <- data.frame(
    sample_id       = rownames(pheno),
    inferred_group  = inferred,
    final_group     = inferred,  # editable
    title           = if ("title" %in% colnames(pheno)) as.character(pheno$title) else NA,
    source_name_ch1 = if ("source_name_ch1" %in% colnames(pheno)) as.character(pheno$source_name_ch1) else NA,
    characteristics = char_text,
    stringsAsFactors = FALSE
  )
  
  sheet
}

# ---------- main ----------
all_gses <- unique(c(GSE_mrna, GSE_mirna))

for (gse_id in all_gses) {
  message("\n==============================")
  message("Sample sheet step for: ", gse_id)
  
  out_csv <- file.path(DIR_META, paste0(gse_id, "_sample_sheet.csv"))
  
  if (file.exists(out_csv) && !force_recreate) {
    message("Found existing sample sheet (will NOT overwrite): ", out_csv)
    sheet <- read.csv(out_csv)
    sheet <- normalize_sheet(sheet)
    # re-save normalized (optional but useful)
    write.csv(sheet, out_csv, row.names = FALSE)
    validate_counts_or_stop(sheet, gse_id)
    next
  }
  
  message("Creating NEW sample sheet from GEO metadata: ", gse_id)
  sheet <- create_sheet_from_geo(gse_id)
  sheet <- normalize_sheet(sheet)
  write.csv(sheet, out_csv, row.names = FALSE)
  message("Wrote: ", out_csv)
  
  validate_counts_or_stop(sheet, gse_id)
}

message("\nStage 03 complete: sample sheets are ready and validated.")
