get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
pkg_root <- get_script_dir()
code_dir <- file.path(pkg_root, "code")
log_dir  <- file.path(pkg_root, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, "replication_log.txt")
file.create(log_file)
setwd(tempdir())

rscript_path <- file.path(R.home("bin"), "Rscript")

log_append <- function(text) {
  for (attempt in 1:5) {
    ok <- tryCatch({
      cat(text, file = log_file, append = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(invisible())
    Sys.sleep(0.5)
  }
  warning("Could not write to log file after 5 attempts")
}

run_r <- function(r_file, target) {
  line   <- strrep("=", 80)
  banner <- paste0(line, "\n", r_file, ": ", target, "\n", line, "\n")
  cat(banner)
  log_append(banner)

  out <- system2(rscript_path, args = shQuote(file.path(code_dir, r_file)),
                 stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  combined <- paste0(paste(out, collapse = "\n"), "\n\n")
  cat(combined)
  log_append(combined)

  if (status != 0) quit(status = 1)
}

steps <- list(
  list(file = "00_prepare_ukds_data.R",              target = "Data preparation, UK Data Service"),
  list(file = "01_prepare_ukds_analysis_file.R",     target = "Data preparation, UK Data Service analysis file"),
  list(file = "02_prepare_natcen_data.R",            target = "Data preparation, NatCen"),
  list(file = "03_figure1.R",                        target = "Figure 1"),
  list(file = "04_table2_figures2to4.R",             target = "Table 2 (Panels A, C); Tables A4, A5 (Panels A, C), A6 (UK Data Service), A8 (IV column); Figures 2, 3, 4"),
  list(file = "05_table2_panelB.R",                  target = "Table 2 (Panel B); Table A5 (Panel B)"),
  list(file = "06_tablesA2_A3.R",                    target = "Tables A2, A3"),
  list(file = "07_figuresA1toA4.R",                  target = "Figures A1, A2, A3, A4"),
  list(file = "08_proportional_hazards_C1.R",        target = "Supplement C.1 text statistics"),
  list(file = "09_missing_data_weights_figureA5.R",  target = "Figure A5"),
  list(file = "10_tableA6_natcen.R",                 target = "Table A6 (NatCen)"),
  list(file = "11_sample_overlap_D1.R",              target = "Supplement D.1 text statistics"),
  list(file = "12_tableA7.R",                        target = "Table A7"),
  list(file = "13_tableA8.R",                        target = "Table A8")
)

for (s in steps) run_r(s$file, s$target)
