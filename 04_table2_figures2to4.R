suppressPackageStartupMessages({
  library(tidyverse); library(survival); library(survminer); library(patchwork); library(car)
})

set.seed(321)

B_MAIN <- 1000
B_SUB  <- 1000

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

nc <- prepare_df(readRDS(file.path(tempdata, "ELSA_bmonth-year.rds")))
sn <- prepare_df(readRDS(file.path(tempdata, "ELSA_sn8346_bmonth-year.rds")))

# Analytic sample sizes
parent_n <- function(ds, ds_name) {
  a <- make_analytical(ds)
  tibble(
    dataset = ds_name,
    n_father = sum(!is.na(a$flongevity2) & !is.na(a$fdead)),
    n_mother = sum(!is.na(a$mlongevity2) & !is.na(a$mdead))
  )
}
cat("\nAnalytic sample sizes (Table A2, Table A3)\n")
print(as.data.frame(bind_rows(parent_n(sn, "SN 5050/8346"), parent_n(nc, "NatCen"))))

# First stage (Table 2 Panel C, Table A4)
fs_spec <- list(
  list(cohort = "cohort1", iv = "post47", rv = "R_months47", reform = "1947"),
  list(cohort = "cohort2", iv = "post72", rv = "R_months72", reform = "1972")
)

first_stage_rows <- list()
for (nm in c("sn", "nc")) {
  ds       <- get(nm)
  ds_label <- if (nm == "sn") "SN 5050/8346" else "NatCen"

  for (s in fs_spec) {
    fs <- fit_first_stage(ds %>% filter(.data[[s$cohort]]), first_stage_formula(s$iv, s$rv))
    assign(paste0(nm, "_", s$cohort), fs$data)

    lh <- linearHypothesis(fs$model, paste(s$iv, "= 0"))
    cf <- summary(fs$model)$coefficients[s$iv, ]
    ci <- confint(fs$model, s$iv, level = 0.95)

    first_stage_rows[[paste(nm, s$cohort)]] <- tibble(
      dataset = ds_label,
      reform  = s$reform,
      coef    = round(cf[["Estimate"]], 4),
      se      = round(cf[["Std. Error"]], 4),
      ci_lo   = round(ci[1], 4),
      ci_hi   = round(ci[2], 4),
      p       = signif(cf[["Pr(>|t|)"]], 3),
      F_stat  = round(lh$F[2], 2),
      N       = nobs(fs$model)
    )
  }
}

first_stage <- bind_rows(first_stage_rows)
cat("\nFirst stage (Table 2 Panel C, Table A4)\n")
print(as.data.frame(first_stage))
write_csv(first_stage, file.path(output_path, "first_stage.csv"))

# Table 2 Panels A and C; Table A5 Panels A and C; Table A8
hr_row <- function(fit, term) {
  s  <- summary(fit)
  cf <- s$coefficients[term, ]
  ci <- s$conf.int[term, ]
  tibble(
    HR     = round(unname(cf["exp(coef)"]), 4),
    lo     = round(unname(ci["lower .95"]), 4),
    hi     = round(unname(ci["upper .95"]), 4),
    p      = signif(unname(cf["Pr(>|z|)"]), 3),
    N      = fit$n,
    events = fit$nevent,
    AIC    = round(AIC(fit), 1),
    BIC    = round(BIC(fit), 1)
  )
}

tab2_rows <- list()
for (nm in c("sn", "nc")) {
  ds_label <- if (nm == "sn") "SN 5050/8346" else "NatCen"

  for (cl in cells) {
    d <- get(paste0(nm, "_", cl$cohort))
    covs <- "fem + par_class2 + nsibs_cat"

    m_naive <- coxph(as.formula(paste0("Surv(", cl$time, ", ", cl$event,
                                       ") ~ schlage + ", cl$rv, " + ", covs)), data = d)
    m_itt   <- coxph(as.formula(paste0("Surv(", cl$time, ", ", cl$event,
                                       ") ~ ", cl$iv, " + ", cl$rv, " + ", covs)), data = d)
    m_2sri  <- coxph(as.formula(paste0("Surv(", cl$time, ", ", cl$event,
                                       ") ~ schlage + presid + ", cl$rv, " + ", covs)), data = d)

    bs <- run_bootstrap(d, cl$time, cl$event, cl$iv, cl$rv, B = B_MAIN)

    tab2_rows[[paste(nm, cl$parent, cl$reform)]] <- bind_rows(
      hr_row(m_naive, "schlage") %>%
        mutate(model = "Observational (naive Cox)", boot_lo = NA, boot_hi = NA),
      hr_row(m_itt, cl$iv) %>%
        mutate(model = "ITT (reduced form)", boot_lo = NA, boot_hi = NA),
      hr_row(m_2sri, "schlage") %>%
        mutate(model = "IV (2SRI Cox)",
               boot_lo = unname(bs["lo"]), boot_hi = unname(bs["hi"]))
    ) %>%
      mutate(dataset = ds_label, parent = cl$parent, reform = cl$reform,
             .before = 1)
  }
}

tab2 <- bind_rows(tab2_rows) %>%
  select(dataset, parent, reform, model, HR, lo, hi, boot_lo, boot_hi, p, N, events, AIC, BIC)

cat("\nTable 2, Panels A and C (SN 5050/8346); Table A5, Panels A and C (both datasets); Table A8\n")
print(as.data.frame(tab2))
write_csv(tab2, file.path(output_path, "table2.csv"))

# Figure 2
fig2_data <- sn %>%
  filter(!is.na(MoB)) %>%
  group_by(MoB) %>%
  summarise(
    fgd = mean(fagedied, na.rm = TRUE),
    mgd = mean(magedied, na.rm = TRUE),
    fle = mean(falive,   na.rm = TRUE),
    mle = mean(malive,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(between(as.integer(format(MoB, "%Y")), 1920, 1980))

COL_AGE   <- "#4472C4"
COL_ALIVE <- "#CC3300"
LAB_AGE   <- "Mean age at death (left axis)"
LAB_ALIVE <- "Proportion still alive (right axis)"
ALIVE_SCALE <- 100

make_fig2_panel <- function(df, gd, le, title_label) {
  ggplot(df, aes(x = MoB)) +
    geom_point(aes(y = .data[[gd]], colour = LAB_AGE), alpha = 0.25, size = 0.8) +
    geom_point(aes(y = .data[[le]] * ALIVE_SCALE, colour = LAB_ALIVE),
               alpha = 0.25, size = 0.8) +
    geom_smooth(aes(y = .data[[gd]], colour = LAB_AGE),
                linewidth = 0.8, se = FALSE) +
    geom_smooth(aes(y = .data[[le]] * ALIVE_SCALE, colour = LAB_ALIVE),
                linewidth = 0.8, se = FALSE) +
    scale_colour_manual(values = setNames(c(COL_AGE, COL_ALIVE),
                                          c(LAB_AGE, LAB_ALIVE)),
                        name = NULL, limits = c(LAB_AGE, LAB_ALIVE)) +
    scale_y_continuous(
      name     = "Mean age at death (years)",
      breaks   = c(0, 25, 50, 75, 100),
      limits   = c(0, 100),
      sec.axis = sec_axis(~ . / ALIVE_SCALE, name = "Proportion still alive",
                          breaks = c(0, 0.25, 0.5, 0.75, 1),
                          labels = c("0.00", "0.25", "0.50", "0.75", "1.00"))
    ) +
    scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
    labs(title = title_label, x = "Child's month and year of birth") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(size = 11),
      panel.grid.minor = element_blank(),
      axis.title.y       = element_text(colour = COL_AGE),
      axis.text.y        = element_text(colour = COL_AGE),
      axis.title.y.right = element_text(colour = COL_ALIVE),
      axis.text.y.right  = element_text(colour = COL_ALIVE)
    )
}

fig2 <- (make_fig2_panel(fig2_data, "fgd", "fle", "Father") +
         make_fig2_panel(fig2_data, "mgd", "mle", "Mother")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig2 <- fig2 & guides(colour = guide_legend(
  override.aes = list(alpha = 1, size = 1.8, linewidth = 1)))

save_fig(fig2, "figure2", width = 10, height = 4.4)

# Figure 3
km_plots <- list()

for (P in c("f", "m")) {
  title_label   <- if (P == "f") "Father" else "Mother"
  longevity_var <- paste0(P, "longevity2")
  dead_var      <- paste0(P, "dead")

  df_km <- sn %>%
    filter(!is.na(.data[[longevity_var]]),
           !is.na(.data[[dead_var]]),
           !is.na(schlage4)) %>%
    mutate(schlage4 = factor(schlage4,
                             labels = c("14 y. or younger", "15 y.", "16 y.",
                                        "17 y. or older")))

  surv_obj <- Surv(time = df_km[[longevity_var]], event = df_km[[dead_var]])
  km_fit   <- survfit(surv_obj ~ schlage4, data = df_km)

  km_plots[[P]] <- ggsurvplot(
    km_fit,
    data         = df_km,
    conf.int     = FALSE,
    palette      = c("#7B9CB8", "#CC0000", "#336699", "#CC6633"),
    linetype     = c("dashed", "solid", "dotdash", "solid"),
    legend.title = "",
    legend.labs  = c("14 y. or younger", "15 y.", "16 y.", "17 y. or older"),
    xlim         = c(20, 110),
    break.x.by   = 10,
    xlab         = "Longevity",
    ylab         = "Proportion of people still alive",
    title        = title_label,
    ggtheme      = theme_minimal()
  )$plot
}

fig3 <- km_plots[["f"]] + km_plots[["m"]]
save_fig(fig3, "figure3", width = 10, height = 5)

# Figure 4; Table A6 (UK Data Service columns)
fig4_rows <- list()

for (cl in cells) {
  d <- get(paste0("sn_", cl$cohort))

  for (sg in subgroups) {
    sub <- d[!is.na(d[[sg$var]]) &
               as.character(d[[sg$var]]) == as.character(sg$value), , drop = FALSE]

    rhs  <- setdiff(c("schlage", "presid", cl$rv, "fem", "par_class2", "nsibs_cat"), sg$var)
    f_sg <- as.formula(paste0("Surv(", cl$time, ", ", cl$event, ") ~ ",
                              paste(rhs, collapse = " + ")))

    m <- tryCatch(coxph(f_sg, data = sub), error = function(e) NULL)
    if (is.null(m)) next

    bs <- run_bootstrap(d, cl$time, cl$event, cl$iv, cl$rv, B = B_SUB,
                        subgroup_var = sg$var, subgroup_value = sg$value)

    fig4_rows[[paste(cl$parent, cl$reform, sg$label)]] <- tibble(
      panel    = sg$panel,
      outcome  = paste0(cl$parent, "s, ", cl$reform, " reform"),
      subgroup = sg$label,
      HR       = round(unname(summary(m)$coefficients["schlage", "exp(coef)"]), 4),
      lo       = unname(bs["lo"]),
      hi       = unname(bs["hi"]),
      N        = m$n,
      events   = m$nevent
    )
  }
}

fig4_df <- bind_rows(fig4_rows)
cat("\nFigure 4; Table A6, UK Data Service columns\n")
print(as.data.frame(fig4_df))
write_csv(fig4_df, file.path(output_path, "figure4_estimates.csv"))

plot_panel <- function(df, title_txt) {
  ggplot(df, aes(x = outcome, y = HR, colour = subgroup, shape = subgroup)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(ymin = lo, ymax = hi),
                    position = position_dodge(width = 0.5), size = 0.6) +
    scale_colour_manual(values = c("#336699", "#CC3300")) +
    labs(title = title_txt, x = NULL, y = "Hazard ratio\n(IV)",
         colour = NULL, shape = NULL) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 0),
          legend.position = "bottom")
}

fig4 <- plot_panel(fig4_df %>% filter(panel == "A"), "A") /
        plot_panel(fig4_df %>% filter(panel == "B"), "B")

save_fig(fig4, "figure4", width = 12, height = 9)
