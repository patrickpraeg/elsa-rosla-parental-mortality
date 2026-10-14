suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(rddensity); library(patchwork)
})

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
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
source(file.path(script_dir, "helpers.R"))

save_fig <- function(plot, name, width, height) {
  ggsave(file.path(fig_path, paste0(name, ".pdf")), plot,
         width = width, height = height, device = cairo_pdf)
}

sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

sn <- sn %>%
  mutate(
    par_class2_int = as.integer(as.character(par_class2)),
    par_class2_d1  = as.numeric(par_class2_int == 1),
    par_class2_d3  = as.numeric(par_class2_int == 3)
  )

cohort1 <- fit_first_stage(sn %>% filter(cohort1), first_stage_formula("post47", "R_months47"))$data
cohort2 <- fit_first_stage(sn %>% filter(cohort2), first_stage_formula("post72", "R_months72"))$data

specs <- list(
  list(lab = "Father, 1947", d = cohort1, t = "flongevity2", e = "fdead", rv = "R_months47"),
  list(lab = "Father, 1972", d = cohort2, t = "flongevity2", e = "fdead", rv = "R_months72"),
  list(lab = "Mother, 1947", d = cohort1, t = "mlongevity2", e = "mdead", rv = "R_months47"),
  list(lab = "Mother, 1972", d = cohort2, t = "mlongevity2", e = "mdead", rv = "R_months72")
)

# Figure A1: proportional hazards
cat("\nFigure A1: Schoenfeld tests\n")

zph_rows <- list()
zph_fits <- list()

for (s in specs) {
  f <- as.formula(paste0("Surv(", s$t, ", ", s$e, ") ~ schlage + presid + ", s$rv,
                         " + fem + par_class2 + nsibs_cat"))
  m <- coxph(f, data = s$d)
  z <- cox.zph(m)
  zph_fits[[s$lab]] <- z

  tb <- as.data.frame(z$table)
  tb$term  <- rownames(tb)
  tb$model <- s$lab
  zph_rows[[s$lab]] <- tb %>% select(model, term, chisq, df, p)
}

zph <- bind_rows(zph_rows) %>% mutate(chisq = round(chisq, 3), p = signif(p, 3))
print(as.data.frame(zph))
write_csv(zph, file.path(output_path, "schoenfeld.csv"))

a1_titles <- c("Father, 1947" = "Fathers, 1947 reform",
               "Father, 1972" = "Fathers, 1972 reform",
               "Mother, 1947" = "Mothers, 1947 reform",
               "Mother, 1972" = "Mothers, 1972 reform")

cairo_pdf(file.path(fig_path, "figureA1.pdf"), width = 9.333, height = 6.667)
par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
for (s in specs) {
  z <- zph_fits[[s$lab]]
  plot(z[1], resid = FALSE, se = TRUE, lwd = 2,
       col = "#D7301F", lty = 1:2,
       main = sprintf("%s  (%s)", a1_titles[[s$lab]],
                      if (z$table["schlage", "p"] < 0.001) "p < 0.001"
                      else sprintf("p = %.3f", z$table["schlage", "p"])),
       xlab = "Age at death", ylab = expression(beta(t)))
  abline(h = 0, lty = 2, col = "grey40")
}
invisible(dev.off())

# Figure A2: covariate balance
cat("\nFigure A2: covariate balance\n")

balance_rows <- list()

covs <- list(
  list(v = "fem",           lab = "Woman"),
  list(v = "nsibs",         lab = "Number of siblings"),
  list(v = "par_class2_d1", lab = "Upper class origin"),
  list(v = "par_class2_d3", lab = "Working class origin"),
  list(v = "nonwhite",      lab = "Non-white")
)

for (r in list(list(d = cohort1, iv = "post47", rv = "R_months47", ref = "1947"),
               list(d = cohort2, iv = "post72", rv = "R_months72", ref = "1972"))) {
  b <- r$d %>% filter(abs(.data[[r$rv]]) <= 36)
  for (cv in covs) {
    m <- lm(as.formula(paste(cv$v, "~", r$iv, "+", r$rv)), data = b)
    cf <- summary(m)$coefficients
    if (!r$iv %in% rownames(cf)) next
    balance_rows[[paste(r$ref, cv$v)]] <- tibble(
      reform     = r$ref,
      covariate  = cv$lab,
      disc       = round(cf[r$iv, "Estimate"], 4),
      se         = round(cf[r$iv, "Std. Error"], 4),
      p          = signif(cf[r$iv, "Pr(>|t|)"], 3),
      N          = nobs(m)
    )
  }
}

balance <- bind_rows(balance_rows)
print(as.data.frame(balance))
write_csv(balance, file.path(output_path, "balance.csv"))
cat("\nSignificant discontinuities at p<.05:", sum(balance$p < 0.05), "of", nrow(balance), "\n")

plot_binned <- function(data, xvar, yvar, xlim, title, xlab, ylab, bins = 60) {
  df <- data %>% select(x = all_of(xvar), y = all_of(yvar)) %>%
    filter(!is.na(x), !is.na(y), x >= xlim[1], x <= xlim[2])
  if (!nrow(df)) return(ggplot() + labs(title = title))

  bw <- diff(range(df$x)) / bins
  if (!is.finite(bw) || bw <= 0) bw <- 1
  binned <- df %>%
    mutate(x_mid = floor(x / bw) * bw + bw / 2) %>%
    group_by(x_mid) %>% summarise(y = mean(y), .groups = "drop")

  ggplot() +
    geom_col(data = binned, aes(x_mid, y), fill = "#A7C7E7", width = bw * 0.9) +
    geom_smooth(data = df %>% filter(x < 0), aes(x, y),
                method = "lm", formula = y ~ poly(x, 2, raw = TRUE),
                colour = "#CC3300", se = TRUE, linewidth = 0.7) +
    geom_smooth(data = df %>% filter(x >= 0), aes(x, y),
                method = "lm", formula = y ~ poly(x, 2, raw = TRUE),
                colour = "#CC3300", se = TRUE, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 10)
}

x47 <- c(-200, 500); x72 <- c(-500, 200)
g1 <- plot_binned(sn, "R_months47", "fem",           x47, "1947, Gender",       "Running var - month 1947", "Proportion women")
g2 <- plot_binned(sn, "R_months47", "nsibs",         x47, "1947, N. siblings",  "Running var - month 1947", "Mean N. siblings")
g3 <- plot_binned(sn, "R_months47", "par_class2_d1", x47, "1947, Social class", "Running var - month 1947", "Proportion upper class")
g4 <- plot_binned(sn, "R_months47", "nonwhite",      x47, "1947, Ethnicity",    "Running var - month 1947", "Proportion non-white")
g5 <- plot_binned(sn, "R_months72", "fem",           x72, "1972, Gender",       "Running var - month 1972", "Proportion women")
g6 <- plot_binned(sn, "R_months72", "nsibs",         x72, "1972, N. siblings",  "Running var - month 1972", "Mean N. siblings")
g7 <- plot_binned(sn, "R_months72", "par_class2_d1", x72, "1972, Social class", "Running var - month 1972", "Proportion upper class")
g8 <- plot_binned(sn, "R_months72", "nonwhite",      x72, "1972, Ethnicity",    "Running var - month 1972", "Proportion non-white")

figA2 <- (g1 | g2 | g3 | g4) / (g5 | g6 | g7 | g8)
save_fig(figA2, "figureA2", width = 14, height = 7)

# Figure A3: density test
cat("\nFigure A3: density tests\n")

get_bw <- function(obj) { h <- as.numeric(obj$h); h <- h[!is.na(h)]; h[1:2] }

run_density <- function(x, plotRange, title, xlab, outfile) {
  bw  <- rdbwdensity(x, c = 0, p = 3)
  h   <- get_bw(bw)
  rdd <- rddensity(X = x, c = 0, p = 3, q = 3, h = c(h[1], h[2]), all = TRUE)
  capture.output(summary(rdd), file = file.path(output_path, outfile))
  invisible(summary(rdd))
  cat(sprintf("  %s: robust p = %s\n", title, signif(rdd$test$p_jk, 3)))

  p <- rdplotdensity(rdd = rdd, X = x, plotRange = plotRange, title = title,
                     xlabel = xlab, ylabel = "Density", legendTitle = "",
                     plotN = c(20, 20), plotGrid = "es", CItype = "region",
                     CIuniform = FALSE, hist = TRUE, histFillCol = "#A7C7E7",
                     histFillShade = 0.25, histLineCol = "white")
  gg <- if (inherits(p, "ggplot")) p else p$Estplot
  gg + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7)
}

d47 <- run_density(sn$R_months47, c(-200, 500), "1947 Reform",
                   "Running variable: months around April 1933", "rddensity_1947.txt")
d72 <- run_density(sn$R_months72, c(-500, 200), "1972 Reform",
                   "Running variable: months around September 1957", "rddensity_1972.txt")

figA3 <- d47 + d72
save_fig(figA3, "figureA3", width = 12, height = 5)

# Figure A4: bandwidths, placebo cut-offs, estimators
cat("\nFigure A4: bandwidth, placebo and estimator comparison\n")

extract_term <- function(model, term, panel, facet_title, x_label) {
  s  <- summary(model)
  if (inherits(model, "coxph")) {
    ci <- s$conf.int[term, ]
    tibble(panel = panel, facet_title = facet_title, x_label = x_label,
           estimate = unname(ci[1]), conf.low = unname(ci[3]), conf.high = unname(ci[4]))
  } else {
    cf <- s$coefficients[term, ]
    est <- exp(cf["Estimate"]); se <- cf["Std. Error"]
    tibble(panel = panel, facet_title = facet_title, x_label = x_label,
           estimate = unname(est),
           conf.low  = unname(exp(cf["Estimate"] - 1.96 * se)),
           conf.high = unname(exp(cf["Estimate"] + 1.96 * se)))
  }
}

fit_itt <- function(d, t, e, rv, iv, bw) {
  dd <- d %>% filter(abs(.data[[rv]]) <= bw)
  coxph(as.formula(paste0("Surv(", t, ", ", e, ") ~ ", iv, " + ", rv,
                          " + fem + par_class2 + nsibs_cat")), data = dd)
}

make_placebo <- function(d, rv, offset) {
  d %>%
    mutate(newrv = .data[[rv]] - offset,
           pv    = as.integer(.data[[rv]] >= offset)) %>%
    filter(abs(newrv) <= 60)
}
placebo_offsets <- c(-18, -12, -6, 12, 18)

panelA <- list()
panelC <- list()
for (s in specs) {
  d  <- s$d; rv <- s$rv
  iv <- if (grepl("1947", s$lab)) "post47" else "post72"
  ft <- paste0(sub(",.*", "", s$lab), "s, ", sub(".*, ", "", s$lab), " reform")

  for (bw in c(12, 24, 36, 60)) {
    panelA[[paste(s$lab, bw)]] <- extract_term(
      fit_itt(d, s$t, s$e, rv, iv, bw), iv, "A", ft, paste0(bw, " mo."))

    dd_c     <- d %>% filter(abs(.data[[rv]]) <= bw)
    m_slopes <- coxph(as.formula(paste0("Surv(", s$t, ", ", s$e, ") ~ ", iv, " * ", rv,
                                        " + fem + par_class2 + nsibs_cat")), data = dd_c)
    panelC[[paste(s$lab, bw, "same")]] <- extract_term(
      fit_itt(d, s$t, s$e, rv, iv, bw), iv, "C", ft, paste0(bw, " mo.")) %>%
      mutate(spec = "Common slope")
    panelC[[paste(s$lab, bw, "diff")]] <- extract_term(
      m_slopes, iv, "C", ft, paste0(bw, " mo.")) %>%
      mutate(spec = "Different slopes")
  }

  for (off in placebo_offsets) {
    dp  <- make_placebo(d, rv, off)
    mp  <- coxph(as.formula(paste0("Surv(", s$t, ", ", s$e, ") ~ pv + newrv",
                                   " + fem + par_class2 + nsibs_cat")), data = dp)
    lbl <- paste0("Placebo ", ifelse(off > 0, "+", ""), off, "mo.")
    panelA[[paste(s$lab, "placebo", off)]] <- extract_term(mp, "pv", "A", ft, lbl)
  }
}
panelA_df <- bind_rows(panelA)
panelC_df <- bind_rows(panelC) %>%
  mutate(spec = factor(spec, levels = c("Common slope", "Different slopes")))

reshape_long <- function(data, tv, ev, rv) {
  data %>%
    filter(!is.na(.data[[tv]]), !is.na(.data[[ev]]), !is.na(presid)) %>%
    select(all_of(c(tv, ev, rv, "fem", "par_class2", "nsibs_cat", "schlage", "presid"))) %>%
    mutate(id = row_number()) %>%
    uncount(weights = .data[[tv]], .remove = FALSE, .id = "time") %>%
    group_by(id) %>%
    mutate(event = if_else(time == .data[[tv]] & .data[[ev]] == 1, 1L, 0L)) %>%
    ungroup() %>%
    filter(time >= 20) %>%
    rename(running = all_of(rv))
}

panelB <- list()
for (s in specs) {
  ft <- paste0(sub(",.*", "", s$lab), "s, ", sub(".*, ", "", s$lab), " reform")

  mcox <- coxph(as.formula(paste0("Surv(", s$t, ", ", s$e, ") ~ schlage + presid + ", s$rv,
                                  " + fem + par_class2 + nsibs_cat")), data = s$d)
  panelB[[paste(s$lab, "cox")]] <- extract_term(mcox, "schlage", "B", ft, "Cox")

  lng <- reshape_long(s$d, s$t, s$e, s$rv)

  for (lk in c("cloglog", "logit")) {
    m <- glm(event ~ schlage + presid + running + fem + par_class2 + nsibs_cat + factor(time),
             family = binomial(lk), data = lng)
    panelB[[paste(s$lab, lk)]] <- extract_term(
      m, "schlage", "B", ft, if (lk == "cloglog") "Cloglog" else "Logit")
  }
  rm(lng); invisible(gc())
}
panelB_df <- bind_rows(panelB)

figA4_df <- bind_rows(panelA_df, panelB_df, panelC_df)
print(as.data.frame(figA4_df))
write_csv(figA4_df, file.path(output_path, "figureA4_estimates.csv"))

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
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
}

plot_panelC <- function(df, ylab_text) {
  df <- df %>% mutate(facet_title = factor(facet_title, levels = unique(facet_title)),
                      x_label     = factor(x_label, levels = unique(x_label)))
  ggplot(df, aes(x = x_label, y = estimate, colour = spec, shape = spec)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.5) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.08,
                  linewidth = 0.7, position = position_dodge(width = 0.4)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.4)) +
    facet_wrap(~ facet_title, nrow = 1) +
    scale_y_log10(breaks = c(0.1, 0.25, 0.5, 1, 2.5, 5)) +
    coord_cartesian(ylim = c(0.08, 6)) +
    scale_colour_manual(values = c("Common slope" = "grey30", "Different slopes" = "#F8766D")) +
    labs(x = NULL, y = ylab_text, colour = NULL, shape = NULL) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(size = 11, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          legend.position = "bottom")
}

figA4 <- (plot_panel(panelA_df, "Hazard ratio\n(ITT)") + ggtitle("A") +
            theme(plot.title = element_text(face = "bold", size = 18))) /
         (plot_panel(panelB_df, "Hazard ratio\n(IV)") + ggtitle("B") +
            theme(plot.title = element_text(face = "bold", size = 18))) /
         (plot_panelC(panelC_df, "Hazard ratio\n(ITT)") + ggtitle("C") +
            theme(plot.title = element_text(face = "bold", size = 18)))

save_fig(figA4, "figureA4", width = 14, height = 13)
