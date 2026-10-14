suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(patchwork); library(mice)
})

set.seed(20261008)

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
fig_path    <- file.path(pkg_root, "figures")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_path,    showWarnings = FALSE, recursive = TRUE)
source(file.path(script_dir, "helpers.R"))

M    <- 20
B_MI <- 100
B_W  <- 500

sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

# First positive cross-sectional weight
w <- as.matrix(sn[, paste0("wgt", 1:11)]); storage.mode(w) <- "numeric"
first_w <- apply(w, 1, function(x) { i <- which(!is.na(x) & x > 0); if (length(i)) i[1] else NA })
sn$w_first <- w[cbind(seq_len(nrow(w)), first_w)]

# Two-stage estimator and bootstrap
two_stage <- function(d, cl, weighted = FALSE) {
  fs_f <- first_stage_formula(cl$iv, cl$rv)
  ss_f <- as.formula(paste0("Surv(time_, event_) ~ schlage + presid + ", cl$rv,
                            " + fem + par_class2 + nsibs_cat"))
  fs <- if (weighted) lm(fs_f, data = d, weights = w_first, na.action = na.omit)
        else          lm(fs_f, data = d, na.action = na.omit)
  omit <- attr(model.frame(fs), "na.action")
  used <- if (length(omit)) seq_len(nrow(d))[-as.integer(omit)] else seq_len(nrow(d))
  d$presid <- NA_real_
  d$presid[used] <- residuals(fs)
  fit <- if (weighted) coxph(ss_f, data = d, weights = w_first, robust = TRUE)
         else          coxph(ss_f, data = d)
  unname(coef(fit)["schlage"])
}

boot_draws <- function(d, cl, B, weighted = FALSE) {
  est <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    bd <- d[sample(nrow(d), replace = TRUE), ]
    est[b] <- tryCatch(two_stage(bd, cl, weighted), error = function(e) NA_real_)
  }
  est[!is.na(est)]
}

make_window <- function(cl) {
  d <- sn %>% filter(.data[[cl$cohort]]) %>%
    mutate(time_ = .data[[cl$time]], event_ = .data[[cl$event]], iv_ = .data[[cl$iv]])
  as.data.frame(d)
}

# Panel A: multiple imputation
cat("\nFigure A5, Panel A: multiple imputation\n")
mi_rows <- list()
for (cl in cells) {
  d       <- make_window(cl)
  has_out <- !is.na(d$time_) & !is.na(d$event_)
  cc_n <- sum(has_out & complete.cases(d[, c("par_class2", "nsibs_cat", "fem", "schlage", cl$rv)]))

  dd   <- d[has_out, c("time_", "event_", "iv_", cl$rv, "fem", "schlage", "par_class2", "nsibs_cat")]
  dd$H <- nelsonaalen(dd, time_, event_)
  meth <- make.method(dd); meth[] <- ""
  meth["par_class2"] <- "polyreg"; meth["nsibs_cat"] <- "polyreg"
  pred <- make.predictorMatrix(dd); pred[, ] <- 0
  pred["par_class2", c("event_", "iv_", cl$rv, "fem", "schlage", "H", "nsibs_cat")] <- 1
  pred["nsibs_cat",  c("event_", "iv_", cl$rv, "fem", "schlage", "H", "par_class2")] <- 1
  imp <- mice(dd, m = M, maxit = 10, method = meth, predictorMatrix = pred,
              printFlag = FALSE, seed = 1)

  Q <- numeric(M); draws <- c()
  for (i in seq_len(M)) {
    di <- d
    ci <- complete(imp, i)
    di$par_class2[has_out] <- ci$par_class2
    di$nsibs_cat[has_out]  <- ci$nsibs_cat
    Q[i]  <- two_stage(di, cl)
    draws <- c(draws, boot_draws(di, cl, B_MI))
  }
  qd <- quantile(draws, c(0.025, 0.975))

  mi_rows[[length(mi_rows) + 1]] <- tibble(
    parent = cl$parent, reform = cl$reform, N_window = nrow(d), N_complete_case = cc_n,
    miss_class = sum(is.na(d$par_class2[has_out])), miss_sibs = sum(is.na(d$nsibs_cat[has_out])),
    mi_2sri_HR = exp(mean(Q)), mi_2sri_lo = exp(qd[[1]]), mi_2sri_hi = exp(qd[[2]]),
    n_draws = length(draws))
}
mi_res <- bind_rows(mi_rows)
print(as.data.frame(mi_res))
write_csv(mi_res, file.path(output_path, "mice.csv"))

# Panel B: survey weights
cat("\nFigure A5, Panel B: survey weights\n")
w_rows <- list()
for (cl in cells) {
  d_all <- make_window(cl) %>% filter(!is.na(par_class2), !is.na(nsibs_cat))
  d_w   <- d_all[!is.na(d_all$w_first), ]
  est   <- exp(two_stage(d_w, cl, TRUE))
  ci    <- exp(unname(quantile(boot_draws(d_w, cl, B_W, TRUE), c(0.025, 0.975))))
  w_rows[[length(w_rows) + 1]] <- tibble(
    parent = cl$parent, reform = cl$reform,
    N = sum(!is.na(d_w$time_) & !is.na(d_w$event_)),
    iv_HR = est, iv_lo = ci[1], iv_hi = ci[2])
}
w_res <- bind_rows(w_rows)
print(w_res %>% mutate(across(where(is.numeric), ~ round(.x, 2))), width = 200)
write_csv(w_res, file.path(output_path, "weights.csv"))

# Figure A5
facet_lab <- function(parent, reform) paste0(parent, "s, ", reform, " reform")

cc <- read_csv(file.path(output_path, "table2.csv"), show_col_types = FALSE) %>%
  filter(dataset == "SN 5050/8346", model == "IV (2SRI Cox)") %>%
  transmute(parent, reform = as.character(reform), estimate = HR,
            conf.low = boot_lo, conf.high = boot_hi)
mi <- mi_res %>%
  transmute(parent, reform, estimate = mi_2sri_HR, conf.low = mi_2sri_lo, conf.high = mi_2sri_hi)
wt <- w_res %>%
  transmute(parent, reform, estimate = iv_HR, conf.low = iv_lo, conf.high = iv_hi)

panelA <- bind_rows(cc %>% mutate(x_label = "Complete\ncase"),
                    mi %>% mutate(x_label = "Multiple\nimputation")) %>%
  mutate(facet_title = facet_lab(parent, reform))
panelB <- bind_rows(cc %>% mutate(x_label = "Unweighted"),
                    wt %>% mutate(x_label = "Weighted")) %>%
  mutate(facet_title = facet_lab(parent, reform))

write_csv(bind_rows(panelA %>% mutate(panel = "A"), panelB %>% mutate(panel = "B")),
          file.path(output_path, "figureA5_estimates.csv"))

plot_panel <- function(df, ylab_text) {
  df <- df %>% mutate(facet_title = factor(facet_title, levels = unique(facet_title)),
                      x_label     = factor(x_label, levels = unique(x_label)))
  ggplot(df, aes(x = x_label, y = estimate)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08,
                  linewidth = 0.7, colour = "#F8766D") +
    geom_point(shape = 1, size = 2.6, colour = "grey30") +
    facet_wrap(~ facet_title, nrow = 1) +
    scale_y_log10(breaks = c(0.1, 0.25, 0.5, 1, 2.5, 5)) +
    coord_cartesian(ylim = c(0.08, 6)) +
    labs(x = NULL, y = ylab_text) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5, lineheight = 0.95))
}
fig <- (plot_panel(panelA, "Hazard ratio\n(IV)") + ggtitle("A") +
          theme(plot.title = element_text(face = "bold", size = 18))) /
       (plot_panel(panelB, "Hazard ratio\n(IV)") + ggtitle("B") +
          theme(plot.title = element_text(face = "bold", size = 18)))

ggsave(file.path(fig_path, "figureA5.pdf"), fig, width = 14, height = 9, device = cairo_pdf)
