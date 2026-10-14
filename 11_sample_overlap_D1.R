suppressPackageStartupMessages(library(tidyverse))

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir <- get_script_dir()
tempdata   <- file.path(dirname(script_dir), "data/tempdata")
source(file.path(script_dir, "helpers.R"))

nc <- prepare_df(readRDS(file.path(tempdata, "ELSA_bmonth-year.rds")))
sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

report <- function(a, b, label) {
  ia <- unique(a$idauniq); ib <- unique(b$idauniq)
  both <- intersect(ia, ib)
  cat("\n===", label, "===\n")
  cat(sprintf("NatCen N            : %6d\n", length(ia)))
  cat(sprintf("UKDS N              : %6d\n", length(ib)))
  cat(sprintf("In both             : %6d  (%.1f%% of NatCen, %.1f%% of UKDS)\n",
              length(both), 100 * length(both) / length(ia),
              100 * length(both) / length(ib)))
  cat(sprintf("NatCen only         : %6d  (%.1f%% of NatCen)\n",
              length(setdiff(ia, ib)), 100 * length(setdiff(ia, ib)) / length(ia)))
  cat(sprintf("UKDS only           : %6d  (%.1f%% of UKDS)\n",
              length(setdiff(ib, ia)), 100 * length(setdiff(ib, ia)) / length(ib)))
  cat(sprintf("Union               : %6d\n", length(union(ia, ib))))
}

cat("\nSupplement D.1: sample overlap\n")
report(nc, sn, "Prepared datasets")

nc_desc <- make_analytical(nc)
sn_desc <- make_analytical(sn)
report(nc_desc, sn_desc, "Analytic samples")

ia <- unique(nc_desc$idauniq); ib <- unique(sn_desc$idauniq)
grp <- list(
  both        = nc_desc %>% filter(idauniq %in% intersect(ia, ib)),
  natcen_only = nc_desc %>% filter(!idauniq %in% ib),
  ukds_only   = sn_desc %>% filter(!idauniq %in% ia)
)
cat("\n=== Mean year of birth by group (analytic samples) ===\n")
for (g in names(grp)) {
  cat(sprintf("%-12s N = %6d   mean byear = %7.1f\n",
              g, nrow(grp[[g]]), mean(grp[[g]]$byear, na.rm = TRUE)))
}
