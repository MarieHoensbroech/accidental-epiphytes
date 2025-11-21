## 3) Predict suitability ---------------------
# ============================================================
# 1) ROBUST ALIGNMENT (the long version that worked) + CHECKS
#    - Uses terra::project explicitly (avoids projpred clash)
#    - Crops each source in native CRS, checks values,
#      projects to a fixed template, writes files, stacks them
#    - Includes a sanity-check block at the end
# ============================================================

library(terra)
library(sf)
library(stringr)
library(brms)
library(tidyterra)
library(tidyverse)

terra::terraOptions(progress = 1)

set.seed(1)
options(mc.cores = parallel::detectCores())
options(contrasts = c("contr.sum", "contr.poly"))  # sum-to-zero contrasts

out_dir <- "outputs/09_niche_model"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

env_joined <- fread("outputs/09_niche_model/model_data.csv")
dat <- fread("outputs/09_niche_model/model_data.csv")

mean_r <- rast(file.path(out_dir,"suitability_mean.tif"))
sd_r <- rast(file.path(out_dir,"suitability_sd.tif"))
# ---- 1) BBOX & INPUTS --------------------------------------
bbox_eu <- st_as_sfc(st_bbox(c(xmin = -10, ymin = 35, xmax = 15, ymax = 72),
                             crs = st_crs(4326)))
bb <- terra::vect(bbox_eu)   # SpatVector in EPSG:4326

# CHELSA files (adjust root if needed)
chelsa_files <- list.files("data/environmental/CHELSA/bio/",
                           pattern = glob2rx("CHELSA_*.tif"),
                           full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

chelsa_selected <- chelsa_files[str_detect(
  chelsa_files,
  paste0(c("CHELSA_bio7","CHELSA_bio3","CHELSA_bio6","CHELSA_bio15",
           "CHELSA_npp","CHELSA_vpd_max","CHELSA_sfcWind_max"), collapse = "|")
)]
stopifnot(length(chelsa_selected) > 0)

# Load rasters (as you did)
r_chelsa <- terra::rast(chelsa_selected)  # 7 layers, 0.008333° grid
r_access <- terra::rast("data/environmental/traveltime/access.min.NA.tif")
r_age    <- terra::rast("data/environmental/forest.age/2025417153414222_BGIForestAgeMPIBGC1.0.0.nc")[["ForestAge_TC000"]]
r_lai    <- terra::rast("data/environmental/lai/lai_july_filled.tif")

# Put in a named list for clear diagnostics
r_list <- list(CHELSA = r_chelsa, ACCESS = r_access, AGE = r_age, LAI = r_lai)

# ---- 2) TEMPLATE from CHELSA grid cropped to bbox -----------
chelsa_bb <- terra::project(bb, terra::crs(r_chelsa))                # IMPORTANT: terra::project
tmpl      <- terra::crop(r_chelsa[[1]], chelsa_bb, snap = "out")     # stable 0.008333° grid
stopifnot(terra::ncell(tmpl) > 0, inherits(tmpl, "SpatRaster"))

# ---- 3) ALIGN FUNCTION (with diagnostics) -------------------
pick_method <- function(r) if (any(terra::is.factor(r))) "near" else "bilinear"

out_dir <- "outputs/_tmp_align"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

align_one <- function(r, nm, i) {
  message(sprintf("\n---- [%s] source summary ----", nm))
  message("  - CRS: ", terra::crs(r))
  message("  - Ext: ", paste(format(terra::ext(r)), collapse = ", "))
  message("  - nlyr: ", nlyr(r), " | res: ", paste(terra::res(r), collapse = ", "))
  message("  - hasValues (quick): ", tryCatch(as.character(terra::hasValues(r)), error = function(e) "error"))
  
  # A) Crop to bbox in native CRS
  bb_r <- tryCatch(terra::project(bb, terra::crs(r)), error = function(e) NULL)
  if (is.null(bb_r)) { message("  ! Failed to project bbox to source CRS -> skip"); return(NULL) }
  
  r_crop <- tryCatch(terra::crop(r, bb_r, snap = "out"),
                     error = function(e) { message("  ! crop error: ", e$message); NULL })
  if (is.null(r_crop)) return(NULL)
  
  if (terra::ncell(r_crop) == 0 || terra::nrow(r_crop) == 0 || terra::ncol(r_crop) == 0) {
    message("  ! crop produced 0x0 -> skip"); return(NULL)
  }
  
  # B) Quick non-NA check (forces minimal read)
  has_any <- tryCatch(
    nrow(terra::spatSample(r_crop, 1, method = "regular", na.rm = TRUE)) > 0,
    error = function(e) FALSE
  )
  message("  - has at least one non-NA cell after crop? ", has_any)
  if (!has_any) { message("  ! all NA after crop -> skip"); return(NULL) }
  
  # C) Project directly to the template grid (writes to file)
  meth <- pick_method(r_crop)
  fn   <- file.path(out_dir, sprintf("aligned_%02d_%s.tif", i, nm))
  ok <- tryCatch({
    terra::project(
      r_crop, tmpl, method = meth,
      filename = fn, overwrite = TRUE,
      wopt = list(
        datatype = if (meth == "near") "INT2S" else "FLT4S",
        gdal     = c("COMPRESS=LZW", if (meth == "bilinear") "PREDICTOR=2" else NULL)
      )
    )
    TRUE
  }, error = function(e) {
    message("  ! project->file error: ", e$message); FALSE
  })
  
  if (!ok || !file.exists(fn)) { message("  ! no output file written -> skip"); return(NULL) }
  
  # D) Verify projected file has some finite values
  r_chk <- tryCatch(terra::rast(fn), error = function(e) NULL)
  if (is.null(r_chk)) { message("  ! cannot read output -> skip"); return(NULL) }
  
  any_finite <- tryCatch({
    s <- terra::spatSample(r_chk, 10, method = "regular", na.rm = TRUE)
    nrow(s) > 0
  }, error = function(e) FALSE)
  message("  - projected file has some finite values? ", any_finite)
  if (!any_finite) { message("  ! output appears all-NA -> delete & skip");
    try(unlink(fn), silent = TRUE); return(NULL) }
  
  message("  ✓ aligned -> ", fn)
  fn
}

# ---- 4) RUN ALIGNMENT + STACK --------------------------------
r_files_aligned <- mapply(function(r, nm, i) align_one(r, nm, i),
                          r_list, names(r_list), seq_along(r_list),
                          SIMPLIFY = FALSE, USE.NAMES = FALSE)

r_files_aligned <- unlist(Filter(Negate(is.null), r_files_aligned), use.names = FALSE)

if (length(r_files_aligned) == 0) {
  message("\n*** No aligned layers were produced. Quick intersections:")
  for (nm in names(r_list)) {
    r0 <- r_list[[nm]]
    bb_in_src <- tryCatch(terra::project(bb, terra::crs(r0)), error = function(e) NULL)
    inters <- if (is.null(bb_in_src)) FALSE else {
      tryCatch(terra::relate(terra::ext(r0), terra::ext(bb_in_src), "intersects"),
               error = function(e) FALSE)
    }
    message(sprintf(" - %-7s | CRS same as tmpl? %-5s | bbox intersects? %-5s",
                    nm, terra::same.crs(r0, tmpl), inters))
  }
  stop("No layers overlapped the template bbox with valid values.")
}

r_stack <- terra::rast(r_files_aligned)
print(r_stack)

# ============================================================
# 5) SANITY CHECKS -------------------------------------------
#    - all layers aligned to template grid?
#    - same CRS/ext/res?
#    - % non-NA cells per layer
#    - quick summary
# ============================================================

# A) CRS/geometry checks
stopifnot(terra::same.crs(r_stack, tmpl))
# extents and resolution of the stack should match the template
stopifnot(all.equal(terra::ext(r_stack), terra::ext(tmpl)))
stopifnot(all.equal(terra::res(r_stack), terra::res(tmpl)))

# B) Non-NA coverage per layer
nonNA_counts <- (terra::global(!is.na(r_stack), "sum", na.rm = TRUE))
nonNA_frac   <- as.numeric(nonNA_counts$sum) / terra::ncell(r_stack)
names(nonNA_frac) <- names(r_stack)
print(round(nonNA_frac, 4))

# C) Basic stats (min/max) to ensure numeric ranges look sane
print(terra::minmax(r_stack,compute=T))

# (Optional) quick visual check of the first layer
# terra::plot(r_stack[[1]], main = names(r_stack)[1])

# 5) RENAME LAYERS TO MATCH MODEL TERMS --------
# Make sure names(r_stack) exactly match your brms formula predictor names.
# Example (adjust to your formula and the order in r_stack):
names(r_stack)
names(r_stack) <- c("bio15","bio3","bio6","bio7","npp","sfcWind_max","vpd_max",
                    "access_min","sforest_age_fill","LAI")



r_stack[["access_min"]][ r_stack[["access_min"]] < 0 ] <- NA
# If your model has factors (e.g., is_forest, is_island), supply those layers too and ensure types/levels match.

# Assuming you already have a base environmental stack 'r_stack'
# names(r_stack) <- c("bio15","bio3","bio6","bio7","npp","sfcWind_max","vpd_max",
#                     "access_min","sforest_age_fill","LAI")

# --- (A) is_forest and forest_age_fill from a forest_age raster ---
# Suppose you have an unfilled (raw) forest_age raster (years), aligned to r_stack:
forest_age_raw <- r_stack[["sforest_age_fill"]]  # replace with your actual raw forest_age if different
# Derive is_forest: 1 where forest_age_raw is not NA; else 0
is_forest_r <- classify(is.na(forest_age_raw), rcl = matrix(c(0,1, 1,0), ncol=2, byrow=TRUE))
# Explanation: is.na==0 → not NA → forest → 1; is.na==1 → NA → non-forest → 0

# Create forest_age_fill: replace NA with 0 everywhere
forest_age_fill <- forest_age_raw
values(forest_age_fill)[is.na(values(forest_age_fill))] <- 0
names(forest_age_fill) <- "sforest_age_fill"

# --- (B) is_island from Natural Earth polygons (same logic as your points) ---
# 1) Download land polygons
land <- rnaturalearth::ne_download(scale=10, type="land", category="physical", returnclass="sf")

# 2) Equal-area projection and polygon-level area
land_eq <- st_transform(land, 6933) |>
  st_cast("POLYGON") |>
  mutate(poly_id = dplyr::row_number(),
         area_m2 = as.numeric(st_area(geometry)))

# 3) Choose the largest N as mainland (same N as training)
mainland_ids <- land_eq |>
  arrange(desc(area_m2)) |>
  filter(scalerank<3) %>% 
  pull(poly_id)

# 4) Make a mainland = 0, island = 1 raster aligned to r_stack
#    We first tag polygons as mainland/island:
land_eq$isl_flag <- ifelse(land_eq$poly_id %in% mainland_ids, "mainland", "island")

# 5) Bring to the r_stack CRS and rasterize
land_r <- st_transform(land_eq, crs(r_stack))
# Start with a dummy empty raster backbone
template <- r_stack[[1]]
#is_island_r <- rasterize(vect(land_r), template, field = "isl_flag") 
 #levels(is_island_r) <- data.frame(ID = 0:1, label = c("mainland", "island"))  
# Outside land (NA), you may want to set to NA (mask out ocean) or choose 0/1.
# If you want to mask ocean (common): keep NA.
names(is_island_r) <- "is_island"


###PLOT LAYERs-----------
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert to df with cell coords and layer columns
df <- as.data.frame(r_stack, xy = TRUE, na.rm = FALSE)

# Pivot to long format: one row per cell per layer
df_long <- df %>%
  pivot_longer(
    cols = -c(x, y),
    names_to = "layer",
    values_to = "value"
  )

p_rasters<-
ggplot(df_long, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  facet_wrap(~layer, ncol = 4, scales = "free") +
  scale_fill_viridis_c(na.value = "grey85", option = "C") +
  coord_equal() +
  labs(fill = "Value", title = "Environmental Raster Stack (NAs in light grey)") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))


##6) STANDARDISE_---------------
# Extract means/sds from the data you used to fit (dat)
# If you have them stored, use those. Otherwise:
z_params <- function(x) list(mu = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
par_bio3   <- z_params(dat$bio3)
par_bio6   <- z_params(dat$bio6)
par_bio7   <- z_params(dat$bio7)
par_npp    <- z_params(dat$npp)
par_vpd    <- z_params(dat$vpd_max)
par_wind   <- z_params(dat$sfcWind_max)
par_bio15  <- z_params(dat$bio15)
par_LAI    <- z_params(dat$LAI)
par_acc    <- z_params(dat$access_min)
par_age    <- z_params(dat$forest_age_fill)
par_lat    <- z_params(dat$lat_cdeg)  # if you scaled lat; if you only centered, set sd=1


z_std <- function(r, p) (r - p$mu) / p$sd

# Apply the SAME NA-to-0 rules BEFORE standardization where appropriate:
LAI_filled <- r_stack[["LAI"]]
values(LAI_filled)[is.na(values(LAI_filled))] <- 0

# Build a standardized stack
r_std <- c(
  z_std(r_stack[["bio7"]],        par_bio7),
  z_std(r_stack[["npp"]],         par_npp),
  z_std(r_stack[["bio3"]],        par_bio3),
  z_std(r_stack[["access_min"]],  par_acc),
  z_std(r_stack[["bio6"]],        par_bio6),
  z_std(r_stack[["vpd_max"]],     par_vpd),
  z_std(r_stack[["sfcWind_max"]], par_wind),
  z_std(r_stack[["bio15"]],       par_bio15),
  z_std(LAI_filled,               par_LAI),
  # lat_cdeg must be built from latitude and standardized with training params:
  # Create latitude raster:
  {
    lat_r <- init(r_stack[[1]], fun = "y") * (180/pi)  # if CRS is lat/long; otherwise compute lat in degrees appropriately
    z_std(lat_r, par_lat)
  },
  is_island_r,
  is_forest_r,
  z_std(forest_age_fill,          par_age)
)
names(r_std) <- c("bio7","npp","bio3","access_min","bio6","vpd_max","sfcWind_max",
                  "bio15","LAI","lat_cdeg","is_island","is_forest","forest_age_fill")


#7) PEDICT Suitability---------------
# Add offset


## --- INPUTS (assumes r_std, fit_s_re already in memory) ---
# Create a zero offset layer aligned to r_std[[1]]
log_effort <- r_std[[1]]; values(log_effort) <- 0; names(log_effort) <- "log_effort"

# Build predictor stack in the exact order the model expects (all standardized)
pred_stack <- c(
  r_std[["bio7"]], r_std[["npp"]], r_std[["bio3"]], r_std[["access_min"]],
  r_std[["bio6"]], r_std[["vpd_max"]], r_std[["sfcWind_max"]], r_std[["bio15"]],
  r_std[["LAI"]],  r_std[["lat_cdeg"]],
  r_std[["is_island"]], r_std[["is_forest"]], r_std[["forest_age_fill"]],
  log_effort
)

# Names required by the model (strict order)
need <- c("bio7","npp","bio3","access_min","bio6","vpd_max","sfcWind_max",
          "bio15","LAI","lat_cdeg","is_island","is_forest","forest_age_fill","log_effort")
stopifnot(all(need %in% names(pred_stack)))             # all present
pred_stack <- pred_stack[[need]]                        # enforce exact order
stopifnot(identical(names(pred_stack), need))           # and confirm

## --- SETTINGS (tune as needed) ---
NDRAWS_TOTAL <- 200      # total posterior draws
BATCH_DRAWS  <- 25       # draws per brms call
TARGET_CELLS <- 1e5      # target cells per block
OUT_MEAN     <- "outputs/09_niche_model/suitability_mean.tif"
OUT_SD       <- "outputs/09_niche_model/suitability_sd.tif"

terraOptions(todisk = TRUE, memfrac = 0.6)              # help with memory

## --- OUTPUT RASTERS (float32, compressed, same geometry as predictors) ---
tmpl   <- pred_stack[[1]]
mean_r <- rast(tmpl)
sd_r   <- rast(tmpl)
wopt   <- list(datatype="FLT4S", gdal=c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))

writeStart(mean_r, filename = OUT_MEAN, overwrite = TRUE, wopt = wopt)
writeStart(sd_r,   filename = OUT_SD,   overwrite = TRUE, wopt = wopt)
on.exit({ try(writeStop(mean_r), silent=TRUE); try(writeStop(sd_r), silent=TRUE) }, add = TRUE)

## --- GRID SHAPES ---
nrows_r  <- nrow(pred_stack)     # raster rows
ncols_r  <- ncol(pred_stack)     # raster columns (spatial)
nlayers  <- nlyr(pred_stack)     # number of predictor layers
block_rows <- max(1L, floor(TARGET_CELLS / ncols_r))  # rows per block

## --- PREDICTION HELPER: accumulate mean & sd across draw batches (Welford pooling) ---
pred_block <- function(df_block) {
  ok <- stats::complete.cases(df_block)
  n  <- nrow(df_block)
  m  <- numeric(n)        # running mean
  m2 <- numeric(n)        # running sum of squares
  k  <- 0L                # draws seen
  if (!any(ok)) return(list(mean = rep(NA_real_, n), sd = rep(NA_real_, n)))
  
  while (k < NDRAWS_TOTAL) {
    nd <- min(BATCH_DRAWS, NDRAWS_TOTAL - k)
    mu <- posterior_epred(                       # nd x sum(ok)
      fit_s_re,
      newdata    = df_block[ok, , drop = FALSE],
      re_formula = NA,
      ndraws     = nd
    )
    batch_mean <- colMeans(mu)
    batch_var  <- apply(mu, 2, var)
    if (k == 0L) {
      m[ok]  <- batch_mean
      m2[ok] <- batch_var * (nd - 1)
      k <- nd
    } else {
      m_old <- m[ok]; k_old <- k; k_new <- k_old + nd
      delta <- batch_mean - m_old
      m_new <- m_old + delta * (nd / k_new)
      M2_old <- m2[ok]; M2_batch <- batch_var * (nd - 1)
      M2_new <- M2_old + M2_batch + (delta^2) * (k_old * nd / k_new)
      m[ok]  <- m_new
      m2[ok] <- M2_new
      k <- k_new
    }
  }
  sd_vec <- rep(NA_real_, n)
  if (k > 1) sd_vec[ok] <- sqrt(m2[ok] / (k - 1))
  list(mean = m, sd = sd_vec)
}

## --- STREAMED IO OVER ROW BLOCKS ---
readStart(pred_stack)
on.exit(try(readStop(pred_stack), silent=TRUE), add=TRUE)

for (r1 in seq.int(1L, nrows_r, by = block_rows)) {
  r2    <- min(r1 + block_rows - 1L, nrows_r)
  nrows <- r2 - r1 + 1L
  
  # Read (cells_in_block x nlayers) as a matrix; fallback if terra returns a vector
  M <- readValues(pred_stack, row = r1, nrows = nrows, mat = TRUE)
  if (is.null(dim(M))) M <- matrix(M, ncol = nlayers, byrow = FALSE)
  stopifnot(ncol(M) == nlayers)
  
  colnames(M) <- names(pred_stack)
  df_block <- as.data.frame(M)
  
  # Force suitability (effort-invariant): offset must be 0
  if (!"log_effort" %in% names(df_block)) df_block$log_effort <- 0
  df_block$log_effort <- 0
  
  # Predict posterior expected probability for this block
  res <- pred_block(df_block)  # returns vectors of length nrows * ncols_r
  
  # Safety: ensure predicted lengths match cells in this row window
  expected_cells <- nrows * ncols_r
  stopifnot(length(res$mean) == expected_cells, length(res$sd) == expected_cells)
  
  # Write vectors; when v is a vector, supply nrows explicitly
  writeValues(mean_r, res$mean, r1, nrows = nrows)
  writeValues(sd_r,   res$sd,   r1, nrows = nrows)
}

# Finalize files (also handled by on.exit if an error occurs)
writeStop(mean_r)
writeStop(sd_r)



#8) MAP Suitability-------------

r_df <- as.data.frame(mean_r, xy = TRUE)

suitability_map<-
  ggplot(r_df, aes(x = x, y = y, fill = bio7)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Suitability") +
  coord_equal() +
  labs(title = "Predicted Suitability for Accidental Epiphytes") +
  my_theme14

suitability_map

ggsave("outputs/09_niche_model/suitability_map.png", width = 10, height = 8)








##9) PLOT SCENARIO MAPS---------
mk_const <- function(template, val) { r <- template; values(r) <- val; r }

pred_for_scenario <- function(is_forest_val, is_island_val) {
  pred_stack2 <- pred_stack
  pred_stack2[["is_forest"]] <- mk_const(pred_stack2[[1]], is_forest_val)
  pred_stack2[["is_island"]] <- mk_const(pred_stack2[[1]], is_island_val)
  
  df2 <- as.data.frame(pred_stack2, xy=FALSE, cells=FALSE, na.rm=FALSE)
  mu  <- posterior_epred(fit_s_re, newdata=df2, re_formula=NA, offset=df2$log_effort)
  out <- pred_stack2[[1]]; values(out) <- colMeans(mu)
  out
}

pred_forest_mainland <- pred_for_scenario(1, 0)
pred_forest_island   <- pred_for_scenario(1, 1)
pred_open_mainland   <- pred_for_scenario(0, 0)
pred_open_island     <- pred_for_scenario(0, 1)













############
# 6) PREDICT MEAN WITH STREAMING (terra::predict + brms) --------
# Optional: cache required predictors from the model (excludes response)
get_req_terms <- function(model) {
  # Works for common brms formulas; adjust if you use special terms
  f <- formula(model)                # response ~ predictors
  vars <- all.vars(f)
  resp <- all.vars(formula(model, resp = TRUE))[1]
  setdiff(vars, resp)
}
.req <- get_req_terms(fit_s_re)  # e.g., c("bio7","bio3",...,"LAI")

# 1) Coercions/transforms to mirror training (edit if you used scaling/factors)
coerce_predictors <- function(df) {
  # EXAMPLES (uncomment/adjust to match training):
  # df$bio7 <- (df$bio7 - mu_bio7) / sd_bio7
  # df$bio3 <- (df$bio3 - mu_bio3) / sd_bio3
  
  # If you truly have factors in the model and corresponding raster layers:
  # if ("is_forest" %in% names(df)) {
  #   df$is_forest <- factor(df$is_forest,
  #                          levels = c(0,1,2),
  #                          labels = c("baseline","level1","level2"))
  # }
  # if ("is_island" %in% names(df)) {
  #   df$is_island <- factor(df$is_island,
  #                          levels = c(0,1),
  #                          labels = c("mainland","island"))
  # }
  df
}

# 2) SAFE mean on response scale: always return vector length nrow(df)
brms_fun_mean <- function(df, model) {
  # Ensure all required columns exist; if missing, create as all-NA
  missing_cols <- setdiff(.req, names(df))
  if (length(missing_cols)) {
    for (nm in missing_cols) df[[nm]] <- NA_real_
  }
  
  df <- coerce_predictors(df)
  
  # Good rows = complete cases for required predictors
  good <- stats::complete.cases(df[, .req, drop = FALSE])
  out  <- rep(NA_real_, nrow(df))
  if (any(good)) {
    mu <- fitted(model, newdata = df[good, , drop = FALSE],
                 re_formula = NA, summary = TRUE)[, "Estimate"]
    out[good] <- as.numeric(mu)
  }
  out
}

# 3) SAFE uncertainty: posterior SD of expected mean (vector length nrow(df))
draws_n <- 300L
brms_fun_sd <- function(df, model) {
  missing_cols <- setdiff(.req, names(df))
  if (length(missing_cols)) {
    for (nm in missing_cols) df[[nm]] <- NA_real_
  }
  
  df <- coerce_predictors(df)
  
  good <- stats::complete.cases(df[, .req, drop = FALSE])
  out  <- rep(NA_real_, nrow(df))
  if (any(good)) {
    E <- posterior_epred(model, newdata = df[good, , drop = FALSE],
                         re_formula = NA, draws = draws_n)  # matrix: draws x n_good
    out[good] <- apply(E, 2, sd)
  }
  out
}

# 4) Run predictions (streamed)
pred_mean_r <- terra::predict(
  r_stack,
  fun       = brms_fun_mean,
  model     = fit_s_re,
  na.rm     = TRUE,
  filename  = "outputs/09_niche_model/suitability_mean.tif",
  overwrite = TRUE
)

pred_sd_r <- terra::predict(
  r_stack,
  fun       = brms_fun_sd,
  model     = fit_s_re,
  na.rm     = TRUE,
  filename  = "outputs/09_niche_model/suitability_uncertainty.tif",
  overwrite = TRUE
)

print(pred_mean_r)
print(pred_sd_r)








################

brms_fun_mean <- function(df, model) {
  # Coerce factor columns here if your model used them, e.g.:
  # df$is_forest <- factor(df$is_forest, levels = c(0,1,2), labels = c("baseline","level1","level2"))
  # df$is_island <- factor(df$is_island, levels = c(0,1))
  as.numeric(predict(model, newdata = df, re_formula = NA, summary = TRUE)[, "Estimate"])
}

pred_mean_r <- terra::predict(
  r_stack,
  fun    = brms_fun_mean,
  model  = fit_s_re,
  na.rm  = TRUE,
  filename = "outputs/09_niche_model/suitability_mean.tif",
  overwrite = TRUE
)

# 7) PREDICT POSTERIOR MEAN + SD (STREAMED; NO GIANT DATA.FRAME) --------
# Change 'draws' to taste (e.g., 500)
draws_n <- 500L
brms_fun_mean_sd <- function(df, model) {
  # Same coercions if factors are present as in brms_fun_mean()
  pred_mat <- posterior_predict(model, newdata = df, draws = draws_n)  # draws x rows
  mu <- colMeans(pred_mat)
  sd <- apply(pred_mat, 2, sd)
  cbind(mean = as.numeric(mu), sd = as.numeric(sd))
}

pred_mean_sd_r <- terra::predict(
  r_stack,
  fun    = brms_fun_mean_sd,
  model  = fit_s_re,
  na.rm  = TRUE,
  filename = "outputs/09_niche_model/suitability_mean_sd.tif",
  overwrite = TRUE
)
names(pred_mean_sd_r) <- c("suitability_mean", "suitability_sd")

# Also write a separate SD file if desired
writeRaster(pred_mean_sd_r[["suitability_sd"]],
            "outputs/09_niche_model/suitability_uncertainty.tif",
            overwrite = TRUE)

# 8) DONE — QUICK CHECK --------
print(rast("outputs/09_niche_model/suitability_mean.tif"))

## 3.3) Visualization ---------------------
# Plot suitability map
library(ggplot2)
r_df <- as.data.frame(r_pred_mean, xy = TRUE)
ggplot(r_df, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Suitability") +
  coord_equal() +
  labs(title = "Predicted Suitability for Accidental Epiphytes") +
  my_theme14
ggsave("outputs/09_niche_model/suitability_map.png", width = 10, height = 8)



## 4) Interpretation------------

# Desiccation sensitivity: VPD_max is the largest (negative) predictor 
# → high vapor pressure deficit strongly reduces presence, consistent with 
# moisture limitation for accidental epiphytes.

# Winter conditions: bio6 (min temp coldest month) strongly positive 
# → warmer winters favor presence (less freeze stress).

# Seasonality (bio7) positive conditional on bio6: After adjusting 
# for winter minimum temp, additional annual temperature range appears 
# associated with higher presence. Because bio6 and bio7 can be correlated, 
# interpret together; the combination is: warmer winters help, and, given 
# a winter minimum, more range is not necessarily harmful (but check collinearity below).

# Wind exposure: Higher max surface wind is negative, 
# consistent with dislodgement and desiccation risk.

# Canopy structure/moisture proxy: LAI positive, 
# consistent with more shelter and substrates in leafy canopies.

# Accessibility (access_min) strong negative: this likely captures 
# sampling bias/detectability (hard-to-reach sites less sampled/less detection), not ecology per se; including it both as offset(log_effort) (for field sites) and as a covariate can double-count survey effort effects (see “Modeling choices & pitfalls” below).

# Island/mainland: The coefficient is positive for “mainland” vs 
# the reference (“island”), suggesting higher odds on the mainland—not 
# negligible here (CrI excludes 0). Verify factor coding (see next section).

# Latitude: Slight negative trend (northward fewer), but CrI includes 0 slightly (borderline).

# Forest age: Small positive effect 

# NPP, bio15: Not meaningfully different from 0 after conditioning on other covariates.

#________________________----------------

