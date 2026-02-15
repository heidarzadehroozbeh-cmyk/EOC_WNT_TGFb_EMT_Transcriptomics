# scripts/00_setup.R
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c(
  "tidyverse","data.table","matrixStats","pROC","glmnet",
  "rsample","yardstick","doParallel","foreach"
)
to_install <- cran_pkgs[!cran_pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

bioc_pkgs <- c(
  "GEOquery","Biobase","limma","sva","metafor",
  "AnnotationDbi","hgu133plus2.db"
)
to_install_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install_bioc) > 0) BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)

cat("Setup complete.\n")
