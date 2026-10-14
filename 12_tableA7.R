suppressPackageStartupMessages({library(tidyverse); library(survival)})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  getwd()
}
script_dir  <- get_script_dir()
pkg_root    <- dirname(script_dir)
tempdata    <- file.path(pkg_root, "data/tempdata")
output_path <- file.path(pkg_root, "Outputs")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
source(file.path(script_dir, "helpers.R"))

nc <- prepare_df(readRDS(file.path(tempdata, "ELSA_bmonth-year.rds")))
sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

covars <- c("fem", "ngor_cat", "nsibs_cat", "nonwhite", "par_class2")

specs <- list(
  `schlage only`   = "factor(schlage)",
  `preprint Eq. 2` = paste("factor(schlage)", paste(covars, collapse = " + "),
                           sep = " + ")
)

parents <- list(
  list(parent = "Father", time = "flongevity2", event = "fdead"),
  list(parent = "Mother", time = "mlongevity2", event = "mdead")
)

rows <- list()
for (nm in c("nc", "sn")) {
  ds  <- get(nm)
  lab <- if (nm == "sn") "UK Data Service" else "NatCen"
  for (cl in parents) {
    keep <- stats::complete.cases(ds[, c("schlage", covars, cl$time, cl$event)])
    d    <- ds[keep, ]

    fits <- lapply(specs, function(sp) {
      f <- as.formula(paste0("Surv(", cl$time, ", ", cl$event, ") ~ ", sp))
      coxph(f, data = d, id = idauniq, x = TRUE, model = TRUE)
    })

    full <- fits[["preprint Eq. 2"]]
    lp   <- predict(full, type = "lp")
    hr   <- predict(full, type = "risk")
    sla  <- if (is.null(names(lp))) d$schlage else d$schlage[match(names(lp), rownames(d))]
    stopifnot(!anyNA(sla), length(sla) == length(lp))
    rho  <- cor(hr, sla, method = "spearman")
    r2   <- summary(lm(lp ~ factor(sla)))$r.squared

    rows[[paste(nm, cl$parent)]] <- tibble(
      dataset = lab, parent = cl$parent,
      C_edu   = round(unname(summary(fits[["schlage only"]])$concordance[1]), 4),
      C_full  = round(unname(summary(full)$concordance[1]), 4),
      se_full = round(unname(summary(full)$concordance[2]), 4),
      rho     = round(rho, 4),
      R2      = round(r2, 4),
      deaths  = full$nevent, N = full$n
    )
  }
}

res <- bind_rows(rows) %>% mutate(delta = round(C_full - C_edu, 4), .after = C_full)
cat("\nTable A7\n")
print(as.data.frame(res), row.names = FALSE)
write_csv(res, file.path(output_path, "tableA7.csv"))
