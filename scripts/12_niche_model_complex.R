#12_niche_model_complex.R
#Find out where accidental epiphytes are likely to occur 
#based on environmental predictors identified in step 4
#

# 0) SETUP & GLOBAL OPTIONS ----------------------------------


suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(projpred)
  library(forcats)
  library(data.table)
  library(dplyr)
})

set.seed(1)
options(mc.cores = parallel::detectCores())
options(contrasts = c("contr.sum", "contr.poly"))  # sum-to-zero contrasts

out_dir <- "outputs/11_niche_model"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

env_joined <- fread("outputs/11_niche_model/model_data.csv")
dat <- fread("outputs/11_niche_model/model_data.csv")


# 1) Data prep -------------------------

# Read and convert to a plain tibble (avoids data.table/dplyr pitfalls)
dat_raw <- fread("outputs/11_niche_model/model_data.csv") %>%
  as_tibble()

# Global contrasts (sum-to-zero helpful with factors)
options(contrasts = c("contr.sum", "contr.poly"))

# Minimal sanitation + computed fields
dat <- dat_raw %>%
  mutate(
    # Source factor levels (explicit)
    source_type = factor(source_type, levels = c("field","inat","background")),
    # 0/1 numerics with NA -> 0
    is_field    = as.numeric(replace_na(is_field, 0)),
    is_forest   = as.numeric(replace_na(is_forest, 0)),
    
    
    # Island factor (baseline = mainland)
    is_island   = fct_relevel(as.factor(is_island), "0"),
    
    # Ensure forest_age_fill is 0 outside forests
    forest_age_fill = ifelse(is_forest == 1, forest_age_fill, 0),
    
    # Detection offset: more effort => less non-detection (negate log_effort).
    # If log_effort missing, use 0 (log(1)=0)
    log_effort = if_else(is.finite(log_effort), log_effort, 0),
    det_offset = -log_effort
  ) %>%  
  mutate(
    site              = fct_na_value_to_level(as.factor(site), level = "no_site") %>% fct_drop(),
    user_login_lumped = fct_na_value_to_level(as.factor(user_login_lumped), level = "no_user") %>% fct_drop(),
    host_tree_lumped  = fct_na_value_to_level(as.factor(host_tree_lumped), level = "other_host") %>% fct_drop()
  ) %>% 
  
  # Collapse top users to "top" vs "other" (keep levels even if one is absent)
  mutate(
    user_top2 = if_else(user_top %in% c("gz_uol", "marie-ho", "nadjakrebs", "spencer_"),
                        "top", "other"),
    user_top2 = factor(user_top2, levels = c("other", "top"))
  )

# 5) Quick NA check (optional)
message("NA counts (core fields):")
print(colSums(is.na(dat[, c("presence","lon","lat","forest_age_fill","is_forest",
                            "access_min","log_effort","det_offset",
                            "site","user_login_lumped","host_tree_lumped",
                            "source_type","is_island"), drop = FALSE])))

# Ecological linear candidates (as per your variables)
ecol_terms <- c("bio7","npp","bio3","bio6","vpd_max","sfcWind_max","bio15",
                "LAI","lat_cdeg","is_island")

# Safe within-source standardization for numeric ecological covariates 
ws_scale <- function(x) {
  if (!is.numeric(x)) return(x)
  if (sum(is.finite(x)) < 2) return(rep(NA_real_, length(x)))   # avoid all-NA groups
  as.numeric(scale(x, center = TRUE, scale = TRUE))
}

dat <- dat %>%
  group_by(source_type) %>%
  mutate(across(all_of(intersect(ecol_terms, names(dat))),
                ws_scale, .names = "{.col}_ws")) %>%
  ungroup()

# Use *_ws for numeric; keep factor is_island as-is
ecol_terms_ws <- paste0(intersect(ecol_terms, names(dat)), "_ws")
ecol_terms_ws <- unique(c(setdiff(ecol_terms_ws, "is_island_ws"), "is_island"))

# Quick NA check on key vars (info only) 
vars_check <- unique(c("presence","lon","lat","forest_age_fill","is_forest",
                       "access_min","log_effort","det_offset", ecol_terms_ws))
na_info <- tibble(variable = vars_check[vars_check %in% names(dat)],
                  n_na = map_int(variable, ~ sum(is.na(dat[[.x]]))))
message("NA counts (selected variables):")
print(na_info %>% arrange(desc(n_na)))


str(dat)
# 2) Weights -----------------
# field rows get more influence

w_mult <- tibble(
  source_type = factor(c("field","inat","background"),
                       levels = c("field","inat","background")),
  mult        = c(3.0, 1.0, 0.5)   # tune: 3x field, 1x inat, 0.5x background
)

dat <- dat %>%
  left_join(w_mult, by = "source_type") %>%
  mutate(mult = ifelse(is.na(mult), 1.0, mult),
         w_raw = mult,
         w     = w_raw / mean(w_raw, na.rm = TRUE))

# Expose as vector for brms;
w <- dat$w

# Helper: fit brms safely with weights; if the local stack errors on `weights`,
# refit without weights (and message).
fit_brms_safe <- function(...) {
  dots <- list(...)
  # Try with weights if provided
  try_fit <- try(do.call(brm, dots), silent = TRUE)
  if (!inherits(try_fit, "try-error")) return(try_fit)
  
  msg <- as.character(try_fit)
  if (grepl("unknown arguments: weights", msg, ignore.case = TRUE)) {
    message("WARNING: Your brms/parallel setup rejected `weights`. Refitting WITHOUT weights...")
    dots$weights <- NULL
    try_fit2 <- do.call(brm, dots)
    attr(try_fit2, "weight_warning") <- TRUE
    return(try_fit2)
  } else {
    stop(try_fit)
  }
}


# 3) Grouped K-FOLD folds -------
### Create grouped K-fold cross-validation folds:
# Observations are grouped by source:
#   - Field data: grouped by site
#   - iNaturalist data: grouped by user
#   - Background points: grouped into pseudo-groups
# Why? If records from the same site or user appear in both training and test folds,
# the model could overestimate performance due to shared context.
# Grouped K-fold ensures all records from the same group stay in the same fold,
# making validation more realistic for spatial/ecological data.
# This affects the model by producing less optimistic accuracy estimates and
# better simulating prediction on unseen sites/users.
dat <- dat %>%
  mutate(
    cv_group = case_when(
      source_type == "field"      ~ paste0("site_", as.character(site)),
      source_type == "inat"       ~ paste0("user_", as.character(user_login_lumped)),
      source_type == "background" ~ paste0("bg_", as.character(row_number() %% 100L)),
      TRUE                        ~ paste0("misc_", row_number() %% 100L)
    )
  )

K <- 5
set.seed(1)
grp_levels <- sample(unique(dat$cv_group))
fold_id_map <- tibble(cv_group = grp_levels, fold = rep(1:K, length.out = length(grp_levels)))
dat <- dat %>% left_join(fold_id_map, by = "cv_group")
folds_vec <- dat$fold
attr(folds_vec, "K") <- K

# 
# # 4) PROJPRED REFERENCE MODEL (BERNOULLI) -------------------
# #     - Ecological linear terms are candidates for selection.
# #     - Forest-age spline, is_forest, and spatial smooth remain fixed.
# ### Define reference model for projection predictive selection (projpred):
# # - projpred starts with a full Bayesian model (reference model) and then
# #   projects its posterior onto smaller submodels to select the most predictive variables.
# # - Use Bernoulli because the response 'presence' is binary (presence/absence).
# # - All ecological covariates are candidates for selection, while forest-age spline,
# #   is_forest, and spatial smooth remain fixed as structural terms.
# 
# form_ref <- bf(
#   as.formula(paste(
#     "presence ~ 0 + Intercept +",
#     paste(ecol_terms_ws, collapse = " + "),
#     "+ s(forest_age_fill, by = is_forest, bs = 'tp') + is_forest + s(lon, lat, bs = 'tp', k = 40)"
#   ))
# )
# 
# # Priors: note Intercept is class "b", coef "Intercept" due to 0 + Intercept
# pri_ref <- c(
#   prior(normal(0, 1.0), class = "b"),                      # ecological linear
#   prior(normal(0, 2.0), class = "b", coef = "Intercept"),  # intercept as population-level b
#   prior(exponential(3), class = "sds")                     # generic spline SD shrinkage
# )
# 
# # 5) Run --------------------
# # Inspect mapping if needed:
# # print(get_prior(form_ref, data = dat, family = bernoulli()))
# 
# fit_ref <- fit_brms_safe(
#   formula = form_ref,
#   data    = dat,
#   family  = bernoulli(link = "logit"),
#   prior   = pri_ref,
#   chains  = 4, iter = 3000, warmup = 1000,
#   weights = w,  # may be removed automatically if your stack errors on weights
#   control = list(adapt_delta = 0.95, max_treedepth = 13),
#   seed    = 1
# )
# 
# saveRDS(fit_ref, file.path(out_dir, "ref_ecology_only_bernoulli.rds"))
# 
# ## 5b) Projpred selection over ecological linear terms only -------
# 
# refmod <- init_refmodel(fit_ref)
# cvsel <- cv_varsel(refmod, formula = search_formula, method = "forward")
# 
# search_formula <- ~ 1 + ecol_terms_ws + s(forest_age_fill) + s(lon, lat)
# 
# 
# # Reference model stays as is
# search_formula <- ~ 1 + ecol_terms_ws
# refmod <- get_refmodel(fit_ref, formula = search_formula)
# 
# # Define selection formula with only ecological linear terms
# sel_formula <- paste("Intercept +", paste(ecol_terms_ws, collapse = " + "))
# 
# # Run projpred selection
# vs <- varsel(refmod, formula = sel_formula)
# 
# set.seed(1)
# vs <- cv_varsel(
#   refmod,
#   method     = "forward",
#   cv_method  = "kfold",
#   folds      = folds_vec,
#   stats      = "elpd",
#   nterms_max = length(ecol_terms_ws),
#   verbose    = TRUE
# )
# 
# size_sel  <- suggest_size(vs, rule = "se")
# sel_terms <- solution_terms(vs, k = size_sel)
# cat("Selected ecological terms:", paste(sel_terms, collapse = ", "), "\n")
# 
# writeLines(paste(sel_terms, collapse = ","), file.path(out_dir, "selected_ecological_terms.txt"))
# saveRDS(vs, file.path(out_dir, "projpred_varsel.rds"))

# _________________________ #####################
# 6) FINAL MODEL --------
# ============================================================
# BRMS: ZERO-INFLATED BINOMIAL WITH SMOOTHS AND PRIORS ------
# ============================================================
rhs_mu <- paste(
  c(ecol_terms,
    "s(forest_age_fill, by = is_forest, bs = 'tp')",
    "is_forest",
    "s(lon, lat, bs = 'tp', k = 40)"),
  collapse = " + "
)

rhs_zi <- paste(
  c("1",
    "source_type",
    "access_min",
    "user_top2",
    "s(site, bs = 're', by = is_field)",
    "s(user_login_lumped, bs = 're', by = is_inat_nontop)",
    "offset(det_offset)"),
  collapse = " + "
)

# --- Build formulas from strings ---
form_mu <- as.formula(paste("presence | trials(1) ~", rhs_mu))
form_zi <- as.formula(paste("~", rhs_zi))

form_final <- bf(form_mu, zi = form_zi)

# --- Define priors (adjust coef names after get_prior output) ---
pri_final <- c(
  # mu (occurrence)
  prior(normal(0, 1.0), class = "b"),
  prior(normal(0, 2.0), class = "Intercept"),
  prior(exponential(3), class = "sds"),
  
  # zi (non-detection)
  prior(normal(0, 0.4), class = "b", dpar = "zi", coef = "source_typeinat"),
  prior(normal(0, 0.4), class = "b", dpar = "zi", coef = "source_typebackground"),
  prior(normal(0.4, 0.4), class = "b", dpar = "zi", coef = "access_min"),
  prior(normal(0, 0.5), class = "b", dpar = "zi", coef = "user_top2top"), # <-- check this name!
  prior(normal(0, 1.5), class = "Intercept", dpar = "zi"),
  prior(exponential(4), class = "sds", dpar = "zi")
)

# --- Inspect prior mapping ---
pri_map <- get_prior(form_final, data = dat, family = zero_inflated_binomial(link = "logit"))
print(pri_map)


pri_final <- c(
  prior(normal(0, 1.0), class = "b"),  # ecological
  prior(normal(0, 2.0), class = "Intercept"),
  prior(exponential(3), class = "sds"),
  prior(normal(0, 0.2), class = "b", dpar = "zi"),  # stronger shrinkage for detection
  prior(normal(0, 1.0), class = "Intercept", dpar = "zi"),
  prior(exponential(4), class = "sds", dpar = "zi")
)


# --- Fit model (optional) ---
fit <- brm(
  formula = form_final,
  data = dat,
  family = zero_inflated_binomial(link = "logit", link_zi="logit"),
  prior = pri_final,
  chains = 4, cores = 4, iter = 3000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)


#########################
# ZERO-INFLATED BERNOULLI (OCCUPANCY SURROGATE)
#     - mu (occurrence): selected ecological terms + forest-age spline + is_forest + spatial
#     - zi (non-detection): source_type, access_min, effort offset, user/site REs


form_final <- bf(
  presence ~ rhs_mu,
  zi ~ rhs_zi
)

pri_final <- c(
  prior(normal(0, 1.0), class = "b"),
  prior(normal(0, 2.0), class = "b", coef = "Intercept"),
  prior(exponential(3), class = "sds"),
  prior(normal(0, 0.4), class = "b", dpar = "zi", coef = "source_typeinat"),
  prior(normal(0, 0.4), class = "b", dpar = "zi", coef = "source_typebackground"),
  prior(normal(0.4, 0.4), class = "b", dpar = "zi", coef = "access_min"),
  prior(normal(0, 0.5), class = "b", dpar = "zi", coef = "user_top2top"),
  prior(normal(0, 1.5), class = "Intercept", dpar = "zi"),
  prior(exponential(4), class = "sds", dpar = "zi")
)

get_prior(form_final, data = dat, family = zero_inflated_binomial(link = "logit"))


fit_final <- fit_brms_safe(
  form_final,
  data   = dat,
  family = bernoulli(link = "logit"),
  prior  = pri_final,
  chains = 4, iter = 4000, warmup = 1000,
  #weights = w,  # may be removed automatically if your stack errors on weights
  control = list(adapt_delta = 0.97, max_treedepth = 13),
  seed   = 1
)

saveRDS(fit_final, file.path(out_dir, "final_ZI_occurrence_detection.rds"))

sink(file.path(out_dir, "final_summary.txt")); print(summary(fit_final)); sink()



## 7) PPC & DIAGNOSTICS --------------------------------------

## 7a) Convergence ----------------
print(rhat(fit_final), digits = 2)
print(neff_ratio(fit_final), digits = 2)

## 7b) Source-stratified PPC ----------------------------------
png(file.path(out_dir, "ppc_by_source.png"), width = 1800, height = 600, res = 150)
par(mfrow = c(1,3))
for (lev in levels(dat$source_type)) {
  idx <- which(dat$source_type == lev & !is.na(dat$source_type))
  if (length(idx) > 50) {
    suppressMessages(pp_check(fit_final, ndraws = 200, type = "bars", resp = "presence", newdata = dat[idx, ]))
    title(main = paste("PPC —", lev))
  }
}
dev.off()

## 7c) Conditional effects on mu for selected ecology--------------------
if (length(sel_terms) > 0) {
  png(file.path(out_dir, "conditional_effects_mu.png"), width = 1600, height = 1200, res = 160)
  plot(conditional_effects(fit_final, effects = sel_terms, dpar = "mu"), ask = FALSE)
  dev.off()
}

## 7d) Marginals on zi-----------------
png(file.path(out_dir, "marginal_effects_zi.png"), width = 1600, height = 900, res = 160)
plot(marginal_effects(fit_final, effects = c("source_type","access_min"), dpar = "zi"), ask = FALSE)
dev.off()


# 8) GROUPED K-FOLD FOR FINAL MODEL -------------------------

folds_list <- lapply(1:K, function(k) which(dat$fold == k))
(kf <- kfold(fit_final, folds = folds_list))
sink(file.path(out_dir, "kfold_results.txt")); print(kf); sink()


# 9) EXPORT COEFFICIENTS ------------------------------------

coef_mu <- fixef(fit_final, dpar = "mu") %>% as.data.frame() %>% rownames_to_column("term")
coef_zi <- fixef(fit_final, dpar = "zi") %>% as.data.frame() %>% rownames_to_column("term")
write.csv(coef_mu, file.path(out_dir, "coef_mu.csv"), row.names = FALSE)
write.csv(coef_zi, file.path(out_dir, "coef_zi.csv"), row.names = FALSE)


# 10) QUICK NA AUDIT  -----------------------

# Rows with ANY NA in variables actually used by the final model:
final_mu_vars <- unique(c("presence","forest_age_fill","is_forest","lon","lat", sel_terms))
final_zi_vars <- c("source_type","access_min","user_top2","det_offset",
                   "is_field","site","user_login_lumped","is_inat_nontop")
final_mu_vars <- final_mu_vars[final_mu_vars %in% names(dat)]
final_zi_vars <- final_zi_vars[final_zi_vars %in% names(dat)]

message("# rows MU would drop: ", sum(!complete.cases(dat[, final_mu_vars, drop = FALSE])))
message("# rows ZI would drop: ", sum(!complete.cases(dat[, final_zi_vars, drop = FALSE])))

dat %>% filter(if_any(all_of(final_mu_vars), is.na))          # rows with NA impacting mu
dat %>% filter(if_any(all_of(final_zi_vars), is.na))          # rows with NA impacting zi
