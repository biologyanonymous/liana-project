# ══════════════════════════════════════════════════════════════════════════════
# Michaelis-Menten vs log-log WAIC comparison — exploratory script
# Runs for: H1 height, H2 TLS height (2019), H2 ground height (2011)
#
# Category predictors match main analysis exactly:
#   H1:      H1_category (factor; reference = Never infested)
#   H2 TLS:  H2_category (factor; reference = remained uninfested)
#   H2 2011: h2_2011_category (factor; reference = remained uninfested)
#
# Random effects:
#   H1 and H2 TLS: (1|species) + (1|sbpltnm)
#   H2 2011:       (1|species) only (most trees lack subplot records)
#
# Outputs saved to:
#   models/mm_waic_comparison/
#   tables/mm_waic_comparison/
# ══════════════════════════════════════════════════════════════════════════════

library(brms)
library(tidyverse)
library(writexl)

# ── Paths ─────────────────────────────────────────────────────────────────────
data_dir   <- "data"
tables_dir <- "tables/mm_waic_comparison"
model_dir  <- "models/mm_waic_comparison"
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir,  recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
dat_h1      <- read_csv(file.path(data_dir, "dat_h1.csv"),      show_col_types = FALSE) |>
  mutate(H1_category = factor(H1_category,
                              levels = c("Never infested", "Lost lianas",
                                         "Persistently infested")),
         species  = factor(species),
         sbpltnm  = factor(sbpltnm))

dat_h2      <- read_csv(file.path(data_dir, "dat_h2.csv"),      show_col_types = FALSE) |>
  mutate(H2_category = factor(H2_category,
                              levels = c("remained uninfested", "gained infestation")),
         species  = factor(species),
         sbpltnm  = factor(sbpltnm))

dat_h2_2011 <- read_csv(file.path(data_dir, "dat_h2_2011.csv"), show_col_types = FALSE) |>
  mutate(h2_2011_category = factor(h2_2011_category,
                                   levels = c("remained uninfested", "gained infestation")),
         species = factor(species))

# ── Shared brms settings ──────────────────────────────────────────────────────
brm_ctrl <- list(adapt_delta = 0.95)
brm_args <- list(chains = 4, cores = 2, iter = 4000,
                 warmup = 1000, seed = 42, silent = 2,
                 control = brm_ctrl)

# ── Priors for MM nonlinear parameters ───────────────────────────────────────
mm_priors <- c(
  prior(normal(3.8, 0.5), nlpar = "Vm"),
  prior(normal(0.5, 0.2), nlpar = "x", lb = 0.01),
  prior(normal(30,  20),  nlpar = "K", lb = 0.01)
)

# ══════════════════════════════════════════════════════════════════════════════
# DATA PREPARATION
# MM models use raw height (m) and DBH (cm) — response not pre-log-transformed
# Category variables kept as factors so reference level is set correctly
# ══════════════════════════════════════════════════════════════════════════════

prep_h1 <- function(d) {
  d |>
    filter(!is.na(height), !is.na(dbh_cm), !is.na(H1_category)) |>
    mutate(height = as.numeric(height),
           dbh    = as.numeric(dbh_cm))
}

prep_h2 <- function(d) {
  d |>
    filter(!is.na(height), !is.na(dbh_cm), !is.na(H2_category)) |>
    mutate(height = as.numeric(height),
           dbh    = as.numeric(dbh_cm))
}

prep_h2_2011 <- function(d) {
  d |>
    filter(!is.na(height_2011), !is.na(dbh_2011), !is.na(h2_2011_category)) |>
    mutate(height = as.numeric(height_2011),
           dbh    = as.numeric(dbh_2011 / 10))  # mm -> cm
}

# ══════════════════════════════════════════════════════════════════════════════
# MM FORMULA FACTORIES
# H1: category = H1_category, random = (1|species) + (1|sbpltnm)
# H2 TLS: category = H2_category, random = (1|species) + (1|sbpltnm)
# H2 2011: category = h2_2011_category, random = (1|species) only
# ══════════════════════════════════════════════════════════════════════════════

mm_resp <- "log(height) ~ Vm + x * log(dbh) - log(K + (dbh^x))"

make_mm_forms <- function(cat_var, re_str) {
  # cat_var: e.g. "H1_category"
  # re_str:  e.g. "(1|species) + (1|sbpltnm)" or "(1|species)"
  list(
    mm2 = brmsformula(as.formula(mm_resp),
      as.formula(paste("Vm ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("x  ~ 1 +", re_str)),
      as.formula(paste("K  ~ 1 +", re_str)),
      nl = TRUE),

    mm3 = brmsformula(as.formula(mm_resp),
      as.formula(paste("x  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("Vm ~ 1 +", re_str)),
      as.formula(paste("K  ~ 1 +", re_str)),
      nl = TRUE),

    mm4 = brmsformula(as.formula(mm_resp),
      as.formula(paste("K  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("Vm ~ 1 +", re_str)),
      as.formula(paste("x  ~ 1 +", re_str)),
      nl = TRUE),

    mm5 = brmsformula(as.formula(mm_resp),
      as.formula(paste("x  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("K  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("Vm ~ 1 +", re_str)),
      nl = TRUE),

    mm6 = brmsformula(as.formula(mm_resp),
      as.formula(paste("Vm ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("K  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("x  ~ 1 +", re_str)),
      nl = TRUE),

    mm7 = brmsformula(as.formula(mm_resp),
      as.formula(paste("Vm ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("x  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("K  ~ 1 +", re_str)),
      nl = TRUE),

    mm8 = brmsformula(as.formula(mm_resp),
      as.formula(paste("Vm ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("x  ~ 1 +", cat_var, "+", re_str)),
      as.formula(paste("K  ~ 1 +", cat_var, "+", re_str)),
      nl = TRUE)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# HELPER: fit one model, cache to disk, report convergence
# ══════════════════════════════════════════════════════════════════════════════

fit_mm <- function(form, dat, label) {
  rds_path <- file.path(model_dir, paste0(label, ".rds"))

  if (file.exists(rds_path)) {
    cat("Loading cached:", label, "\n")
    m         <- readRDS(rds_path)
    rhat_max  <- max(rhat(m), na.rm = TRUE)
    div_trans <- sum(nuts_params(m)$Value[nuts_params(m)$Parameter == "divergent__"])
    converged <- rhat_max < 1.01 & div_trans == 0
    return(list(model = m, converged = converged, label = label))
  }

  cat("Fitting:", label, "\n")
  m <- tryCatch(
    do.call(brm, c(
      list(formula = form, data = dat,
           family  = gaussian(),
           prior   = mm_priors),
      brm_args
    )),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )

  if (is.null(m)) return(list(model = NULL, converged = FALSE, label = label))

  rhat_max  <- max(rhat(m), na.rm = TRUE)
  div_trans <- sum(nuts_params(m)$Value[nuts_params(m)$Parameter == "divergent__"])
  converged <- rhat_max < 1.01 & div_trans == 0

  cat("  Rhat max:", round(rhat_max, 4),
      "| Divergent:", div_trans,
      "| Converged:", converged, "\n")

  saveRDS(m, rds_path)
  list(model = m, converged = converged, label = label)
}

# ══════════════════════════════════════════════════════════════════════════════
# WAIC TABLE BUILDER
# ══════════════════════════════════════════════════════════════════════════════

make_desc <- function(cat_label, params_list) {
  # Generate readable descriptions for each MM form
  c(
    mm2    = paste0("Vm (asymptote) ~ ", cat_label),
    mm3    = paste0("x (scaling) ~ ",    cat_label),
    mm4    = paste0("K (half-sat.) ~ ",  cat_label),
    mm5    = paste0("x + K ~ ",          cat_label),
    mm6    = paste0("Vm + K ~ ",         cat_label),
    mm7    = paste0("Vm + x ~ ",         cat_label),
    mm8    = paste0("Vm + x + K ~ ",     cat_label, " (full MM)"),
    loglog = paste0("log(H) ~ log(DBH) + ", cat_label, " [log-log reference]")
  )
}

build_waic_table <- function(fit_list, loglog_model, hypothesis, descriptions) {
  converged_fits <- keep(fit_list, ~ !is.null(.x$model) && .x$converged)
  failed_fits    <- keep(fit_list, ~ is.null(.x$model) || !.x$converged)

  waic_list <- map(converged_fits, function(f) {
    tryCatch(waic(f$model), error = function(e) NULL)
  }) |> compact()

  waic_loglog <- tryCatch(waic(loglog_model), error = function(e) NULL)
  if (!is.null(waic_loglog)) waic_list[["loglog"]] <- waic_loglog

  if (length(waic_list) < 2) {
    cat("Not enough converged models to compare for", hypothesis, "\n")
    return(NULL)
  }

  cmp    <- loo_compare(waic_list)
  cmp_df <- as.data.frame(cmp)

  converged_tbl <- tibble(
    Hypothesis  = hypothesis,
    model_id    = rownames(cmp_df),
    Description = descriptions[rownames(cmp_df)],
    DWAIC       = signif(as.numeric(cmp_df$elpd_diff) * -2, 3),
    SE_diff     = signif(as.numeric(cmp_df$se_diff)   * -2, 3),
    Converged   = "Yes"
  )

  failed_tbl <- map_dfr(failed_fits, function(f) {
    short_id <- gsub("^(H1|H2_TLS|H2_2011)_", "", f$label)
    tibble(
      Hypothesis  = hypothesis,
      model_id    = f$label,
      Description = descriptions[short_id],
      DWAIC       = NA_real_,
      SE_diff     = NA_real_,
      Converged   = if_else(is.null(f$model),
                            "Error during fitting",
                            "No (Rhat > 1.01 or divergent transitions)")
    )
  })

  bind_rows(converged_tbl, failed_tbl)
}

# ══════════════════════════════════════════════════════════════════════════════
# H1 HEIGHT
# category: H1_category | RE: (1|species) + (1|sbpltnm)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== H1 HEIGHT ======\n")
df_h1       <- prep_h1(dat_h1)
mm_forms_h1 <- make_mm_forms("H1_category", "(1|species) + (1|sbpltnm)")
fits_h1     <- imap(mm_forms_h1, ~ fit_mm(.x, df_h1, paste0("H1_", .y)))

rds_h1_ll <- file.path(model_dir, "H1_loglog.rds")
if (file.exists(rds_h1_ll)) {
  m_h1_loglog <- readRDS(rds_h1_ll)
} else {
  m_h1_loglog <- brm(
    log(height) ~ log(dbh) + H1_category + (1|species) + (1|sbpltnm),
    data = df_h1, family = student(), seed = 42,
    chains = 4, cores = 2, iter = 4000, warmup = 1000,
    control = brm_ctrl, silent = 2)
  saveRDS(m_h1_loglog, rds_h1_ll)
}

desc_h1 <- make_desc("H1_category (Never infested | Lost lianas | Persistently infested)",
                     NULL)
tbl_h1  <- build_waic_table(fits_h1, m_h1_loglog,
                             "H1 — Tree height (reference: Never infested)",
                             desc_h1)

# ══════════════════════════════════════════════════════════════════════════════
# H2 TLS HEIGHT (2019)
# category: H2_category | RE: (1|species) + (1|sbpltnm)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== H2 TLS HEIGHT (2019) ======\n")
df_h2        <- prep_h2(dat_h2)
mm_forms_h2  <- make_mm_forms("H2_category", "(1|species) + (1|sbpltnm)")
fits_h2      <- imap(mm_forms_h2, ~ fit_mm(.x, df_h2, paste0("H2_TLS_", .y)))

rds_h2_ll <- file.path(model_dir, "H2_TLS_loglog.rds")
if (file.exists(rds_h2_ll)) {
  m_h2_loglog <- readRDS(rds_h2_ll)
} else {
  m_h2_loglog <- brm(
    log(height) ~ log(dbh) + H2_category + (1|species) + (1|sbpltnm),
    data = df_h2, family = student(), seed = 42,
    chains = 4, cores = 2, iter = 4000, warmup = 1000,
    control = brm_ctrl, silent = 2)
  saveRDS(m_h2_loglog, rds_h2_ll)
}

desc_h2 <- make_desc("H2_category (remained uninfested | gained infestation)", NULL)
tbl_h2  <- build_waic_table(fits_h2, m_h2_loglog,
                             "H2 TLS — Tree height 2019 (reference: remained uninfested)",
                             desc_h2)

# ══════════════════════════════════════════════════════════════════════════════
# H2 GROUND HEIGHT (2011)
# category: h2_2011_category | RE: (1|species) only
# ══════════════════════════════════════════════════════════════════════════════

cat("\n====== H2 GROUND HEIGHT (2011) ======\n")
df_h2_2011       <- prep_h2_2011(dat_h2_2011)
mm_forms_h2_2011 <- make_mm_forms("h2_2011_category", "(1|species)")
fits_h2_2011     <- imap(mm_forms_h2_2011,
                         ~ fit_mm(.x, df_h2_2011, paste0("H2_2011_", .y)))

rds_h2_2011_ll <- file.path(model_dir, "H2_2011_loglog.rds")
if (file.exists(rds_h2_2011_ll)) {
  m_h2_2011_loglog <- readRDS(rds_h2_2011_ll)
} else {
  m_h2_2011_loglog <- brm(
    log(height) ~ log(dbh) + h2_2011_category + (1|species),
    data = df_h2_2011, family = student(), seed = 42,
    chains = 4, cores = 2, iter = 4000, warmup = 1000,
    control = brm_ctrl, silent = 2)
  saveRDS(m_h2_2011_loglog, rds_h2_2011_ll)
}

desc_h2_2011 <- make_desc("h2_2011_category (remained uninfested | gained infestation)",
                           NULL)
tbl_h2_2011  <- build_waic_table(
  fits_h2_2011, m_h2_2011_loglog,
  "H2 Ground — Tree height 2011 (reference: remained uninfested)",
  desc_h2_2011)

# ══════════════════════════════════════════════════════════════════════════════
# SAVE OUTPUTS
# ══════════════════════════════════════════════════════════════════════════════

if (!is.null(tbl_h1)) {
  write_csv(tbl_h1, file.path(tables_dir, "H1_height_WAIC_comparison.csv"))
  cat("\n-- H1 WAIC comparison --\n"); print(tbl_h1)
}
if (!is.null(tbl_h2)) {
  write_csv(tbl_h2, file.path(tables_dir, "H2_TLS_height_WAIC_comparison.csv"))
  cat("\n-- H2 TLS WAIC comparison --\n"); print(tbl_h2)
}
if (!is.null(tbl_h2_2011)) {
  write_csv(tbl_h2_2011, file.path(tables_dir, "H2_2011_height_WAIC_comparison.csv"))
  cat("\n-- H2 2011 WAIC comparison --\n"); print(tbl_h2_2011)
}

all_tbls <- compact(list(
  "H1_height"      = tbl_h1,
  "H2_TLS_height"  = tbl_h2,
  "H2_2011_height" = tbl_h2_2011
))
if (length(all_tbls) > 0) {
  write_xlsx(all_tbls,
             path = file.path(tables_dir, "all_height_WAIC_comparisons.xlsx"))
  cat("\nAll outputs saved to:", tables_dir, "\n")
}

