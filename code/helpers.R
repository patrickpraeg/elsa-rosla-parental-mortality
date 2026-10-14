# Variable construction
prepare_df <- function(df) {

  invalid_to_na <- function(x) {
    x <- as.numeric(x)
    x[x %in% c(-1, -2, -8, -9)] <- NA_real_
    x
  }

  df <- df %>%
    dplyr::mutate(dplyr::across(c(fage, mage, falive, malive, fagedied, magedied,
                                  byear, sex, nonwhite, nsibs, edend),
                                invalid_to_na))

  df <- df %>% dplyr::mutate(fem = dplyr::if_else(sex == 2, 1L, 0L))

  df <- df %>%
    dplyr::mutate(
      schlage = dplyr::case_when(
        edend %in% c(1, 2) ~ 0L,
        edend == 3          ~ 14L,
        edend == 4          ~ 15L,
        edend == 5          ~ 16L,
        edend == 6          ~ 17L,
        edend == 7          ~ 18L,
        edend == 8          ~ 19L,
        TRUE                ~ NA_integer_
      ),
      schlage15 = dplyr::case_when(
        schlage %in% c(0, 14)    ~ 0L,
        schlage %in% 15:19       ~ 1L,
        TRUE                     ~ NA_integer_
      ),
      schlage16 = dplyr::case_when(
        schlage %in% c(0, 14, 15) ~ 0L,
        schlage %in% 16:19        ~ 1L,
        TRUE                      ~ NA_integer_
      ),
      schlage4 = dplyr::case_when(
        schlage %in% c(0, 14) ~ 1L,
        schlage == 15          ~ 2L,
        schlage == 16          ~ 3L,
        schlage %in% 17:19     ~ 4L,
        TRUE                   ~ NA_integer_
      )
    )

  for (P in c("f", "m")) {
    alive   <- paste0(P, "alive")
    dead    <- paste0(P, "dead")
    longev  <- paste0(P, "longevity")
    longev2 <- paste0(P, "longevity2")
    aged    <- paste0(P, "agedied")
    age     <- paste0(P, "age")

    df[[alive]]  <- dplyr::if_else(as.numeric(df[[alive]]) == 2, 0L,
                                   as.integer(df[[alive]]))
    df[[dead]]   <- dplyr::case_when(
      df[[alive]] == 1 ~ 0L,
      df[[alive]] == 0 ~ 1L,
      TRUE             ~ NA_integer_
    )
    df[[longev]]  <- dplyr::coalesce(as.numeric(df[[aged]]), as.numeric(df[[age]]))
    df[[longev]]  <- dplyr::if_else(!is.na(df[[longev]]) & df[[longev]] < 20, 19,
                                    df[[longev]])
    df[[longev2]] <- dplyr::if_else(!is.na(df[[longev]]) & df[[longev]] > 99, 100,
                                    df[[longev]])
  }

  df$ngor     <- as.integer(df$ngor)
  df$ngor_cat <- addNA(factor(df$ngor, levels = c(8, 1, 2, 3, 4, 5, 6, 7, 9)))

  df <- df %>% dplyr::mutate(nsibs_cat = dplyr::if_else(!is.na(nsibs) & nsibs >= 6,
                                                        5, nsibs))
  df$nsibs_cat <- factor(df$nsibs_cat, levels = c(1, 0, 2, 3, 4, 5))

  df$par_class2 <- factor(as.integer(df$par_class2), levels = c(3, 1, 2, 4))

  df <- df %>%
    dplyr::mutate(
      nmalive_d = as.integer(nmalive_d),
      nfalive_d = as.integer(nfalive_d)
    )

  df %>%
    dplyr::filter(!is.na(schlage), schlage > 0) %>%
    dplyr::mutate(
      cohort1 = !is.na(R_months47) & R_months47 >= -60 & R_months47 <= 60,
      cohort2 = !is.na(R_months72) & R_months72 >= -60 & R_months72 <= 60
    )
}

make_analytical <- function(ds) {
  ds %>% dplyr::filter(
    !is.na(schlage),
    !is.na(sex), !is.na(par_class2), !is.na(nonwhite),
    !is.na(nsibs), !is.na(byear),
    (!is.na(flongevity2) & !is.na(fdead)) | (!is.na(mlongevity2) & !is.na(mdead))
  )
}

# Estimation
cells <- list(
  list(parent = "Father", reform = "1947", cohort = "cohort1", time = "flongevity2",
       event = "fdead", iv = "post47", rv = "R_months47"),
  list(parent = "Father", reform = "1972", cohort = "cohort2", time = "flongevity2",
       event = "fdead", iv = "post72", rv = "R_months72"),
  list(parent = "Mother", reform = "1947", cohort = "cohort1", time = "mlongevity2",
       event = "mdead", iv = "post47", rv = "R_months47"),
  list(parent = "Mother", reform = "1972", cohort = "cohort2", time = "mlongevity2",
       event = "mdead", iv = "post72", rv = "R_months72")
)

subgroups <- list(
  list(panel = "A", var = "fem",        value = 0,   label = "Sons"),
  list(panel = "A", var = "fem",        value = 1,   label = "Daughters"),
  list(panel = "B", var = "par_class2", value = "1", label = "Upper"),
  list(panel = "B", var = "par_class2", value = "3", label = "Working")
)

first_stage_formula <- function(iv, rv) {
  as.formula(paste("schlage ~", iv, "+", rv, "+ fem + par_class2 + nsibs_cat"))
}

fit_first_stage <- function(data, formula) {
  mod  <- lm(formula, data = data, na.action = na.omit)
  used <- as.integer(rownames(model.frame(mod)))
  data$presid       <- NA_real_
  data$presid[used] <- residuals(mod)
  list(model = mod, data = data)
}

run_bootstrap <- function(bdf, surv_time, surv_event, iv, running_var, B,
                          subgroup_var = NULL, subgroup_value = NULL) {
  fs_f <- first_stage_formula(iv, running_var)
  rhs  <- c("schlage", "presid", running_var, "fem", "par_class2", "nsibs_cat")
  if (!is.null(subgroup_var)) rhs <- setdiff(rhs, subgroup_var)
  ss_f <- as.formula(
    paste0("Surv(", surv_time, ", ", surv_event, ") ~ ", paste(rhs, collapse = " + "))
  )

  coefs <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    bd   <- bdf[sample(nrow(bdf), replace = TRUE), ]
    fs   <- try(lm(fs_f, data = bd, na.action = na.omit), silent = TRUE)
    if (inherits(fs, "try-error")) next
    mf   <- model.frame(fs)
    omit <- as.integer(attr(mf, "na.action"))
    used <- if (length(omit)) seq_len(nrow(bd))[-omit] else seq_len(nrow(bd))
    bd$presid       <- NA_real_
    bd$presid[used] <- residuals(fs)

    sd <- if (!is.null(subgroup_var)) {
      bd[!is.na(bd[[subgroup_var]]) &
           as.character(bd[[subgroup_var]]) == as.character(subgroup_value), , drop = FALSE]
    } else bd

    if (nrow(sd) < 50) next

    ss <- tryCatch(coxph(ss_f, data = sd), error = function(e) NULL)
    if (!is.null(ss)) {
      cf <- coef(ss)["schlage"]
      if (is.finite(cf)) coefs[b] <- cf
    }
  }

  n_ok <- sum(!is.na(coefs))
  if (n_ok < B * 0.9)
    warning(sprintf("Only %d/%d bootstrap iterations succeeded", n_ok, B))

  ci <- quantile(coefs, c(0.025, 0.975), na.rm = TRUE)
  c(lo = round(exp(ci[[1]]), 4), hi = round(exp(ci[[2]]), 4), n_ok = n_ok)
}
