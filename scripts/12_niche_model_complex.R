#12_niche_model_complex.R
#Find out where accidental epiphytes are likely to occur 
#based on environmental predictors identified in step 4
#
# Regularization (Version A): Forces detection terms to share variance, 
# making ecological signals more visible.
#
# Horseshoe (Version B): Lets model decide which detection terms matter, 
# reducing overfitting and dominance.
#
# LOO: Compares predictive performance, not just fit.
#
# Post-hoc plots: Show ecological smooths and detection shrinkage visually
# for interpretation.


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

# Ecological linear candidates 
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


#Be double sure there are no NAs in important columns
dat_clean <- dat %>%
  dplyr::filter(!is.na(presence), !is.na(source_type), !is.na(user_top2))


# 1) Build formulas ---
rhs_mu <- paste(
  c("bio7","npp","bio3","bio6","vpd_max","sfcWind_max","bio15","LAI",
    "lat_cdeg","is_island",
    "s(forest_age_fill, by = is_forest, bs = 'tp')",
    "is_forest",
    "s(lon, lat, bs = 'tp', k = 40)"),
  collapse = " + "
)

rhs_zi <- paste(
  c("1","source_type","access_min","user_top2",
    "s(site, bs = 're', by = is_field)",
    "s(user_login_lumped, bs = 're', by = is_inat_nontop)",
    "offset(det_offset)"),
  collapse = " + "
)

form_mu <- as.formula(paste("presence | trials(1) ~", rhs_mu))
form_zi <- as.formula(paste("~", rhs_zi))
form_final <- bf(form_mu, zi = form_zi)

# 2) Version A: Regularization priors ------------
# Strong shrinkage on detection terms (zi) to prevent dominance
pri_regularized <- c(
  prior(normal(0, 1.0), class = "b"),                # ecological terms
  prior(normal(0, 2.0), class = "Intercept"),
  prior(exponential(3), class = "sds"),
  prior(normal(0, 0.2), class = "b", dpar = "zi"),   # strong shrinkage for detection
  prior(normal(0, 1.0), class = "Intercept", dpar = "zi"),
  prior(exponential(4), class = "sds", dpar = "zi")
)

# 3) Version B: Horseshoe priors for detection terms --------
# Allows automatic variable selection for zi predictors
pri_horseshoe <- c(
  prior(normal(0, 1.0), class = "b"),                # ecological terms
  prior(normal(0, 2.0), class = "Intercept"),
  prior(exponential(3), class = "sds"),
  prior(horseshoe(1), class = "b", dpar = "zi"),     # horseshoe for detection
  prior(normal(0, 1.0), class = "Intercept", dpar = "zi"),
  prior(exponential(4), class = "sds", dpar = "zi")
)

# 4) Fit both models -------------
fit_reg <- brm(
  formula = form_final,
  data = dat_clean,
  family = zero_inflated_binomial(link = "logit", link_zi = "logit"),
  prior = pri_regularized,
  chains = 4, cores = 4, iter = 3000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

fit_hs <- brm(
  formula = form_final,
  data = dat_clean,
  family = zero_inflated_binomial(link = "logit", link_zi = "logit"),
  prior = pri_horseshoe,
  chains = 4, cores = 4, iter = 3000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# 5) Compare models using LOO -----------
loo_reg <- loo(fit_reg)
loo_hs <- loo(fit_hs)
loo_compare(loo_reg, loo_hs)

# Interpretation:
# - Lower LOOIC = better predictive performance.
# - If difference > 4, consider the better model strongly preferred.

# 6) Post-hoc visualization ----------
# Ecological effects: conditional smooths and marginal effects
plot(conditional_effects(fit_reg), points = TRUE)
plot(conditional_effects(fit_hs), points = TRUE)

# 7) Detection effects: check shrinkage or selection----------------
mcmc_intervals(as.matrix(fit_reg), pars = grep("^zi_", fit_reg$fit@sim$fnames_oi, value = TRUE))
mcmc_intervals(as.matrix(fit_hs), pars = grep("^zi_", fit_hs$fit@sim$fnames_oi, value = TRUE))


