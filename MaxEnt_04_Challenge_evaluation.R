# ===============================
# MaxEnt Modeling with Presence and Absence Points
# Using terra and dismo
# ===============================

# Load libraries
library(dismo)
library(terra)
library(sp)
library(caret)
library(ENMeval)
library(raster)
library(pROC)

# -------------------------------
# 1. Load data
# -------------------------------

# Raster predictors (RasterStack or SpatRaster)
predictors_rc <- predictors_rc

# Presence points 
presence_data_c <- occurrences_in_chal

##Convert to lon, lat
presence_data_c <- presence_data_c %>%
  dplyr::select(geometry) %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(lon, lat)

# Ensure they're numeric
presence_data_c$lon <- as.numeric(presence_data_c$lon)
presence_data_c$lat <- as.numeric(presence_data_c$lat)

# Convert sf to SpatialPoints
presence_sp_c <- as_Spatial(occurrences_in_chal)

# Absence/background points (data.frame with lon, lat)
absence_data_c <- absent_sites_in_chal

##Convert to lon, lat
absence_data_c <- absence_data_c %>%
  dplyr::select(geometry) %>%
  mutate(
    lon = st_coordinates(geometry)[, 1],
    lat = st_coordinates(geometry)[, 2]
  ) %>%
  st_drop_geometry() %>%
  dplyr::select(lon, lat)

# Ensure they're numeric
absence_data_c$lon <- as.numeric(absence_data_c$lon)
absence_data_c$lat <- as.numeric(absence_data_c$lat)

# Convert to SpatialPoints
presence_sp_c <- SpatialPoints(presence_data_c[, c("lon","lat")],
                             proj4string = CRS("+proj=longlat +datum=WGS84"))

absence_sp_c <- SpatialPoints(
  as.matrix(absence_data_c[, c("lon", "lat")]),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

# ------------------------------
# 4. Align CRS
# ------------------------------
# Transform presence and absence points to raster CRS
presence_spc_utm <- spTransform(presence_sp_c, CRS(proj4string(predictors_rc)))
absence_spc_utm  <- spTransform(absence_sp_c, CRS(proj4string(predictors_rc)))

  
# Convert SpatialPoints to sf
presence_sfc <- st_as_sf(presence_spc_utm)

# Create 500 m buffer
presence_buf <- st_buffer(presence_sf, dist = 500)

# Convert raster to terra SpatRaster if not already
predictors_c_terra <- rast(predictors_rc)

# Crop raster to buffer extent
predictors_c_crop <- crop(predictors_c_terra, vect(presence_buf))
predictors_c_crop <- mask(predictors_c_crop, vect(presence_buf))

#Convert back to RasterStack
predictors_raster_c <- raster::stack(predictors_c_crop)

# ----------------------------
# 2. Fit MaxEnt with replicates
# ----------------------------
maxent_model_c <- maxent(
  x = predictors_raster_c,
  p = presence_spc_utm,
  a = absence_spc_utm,
  args = c("replicates=5", "replicatetype=crossvalidate")
)

# ----------------------------
# 3. Predict each replicate and average
# ----------------------------
pred_full_c <- raster::predict(predictors_raster_c, maxent_model@models[[1]])

pred_listc <- list()
for (i in seq_along(maxent_model@models)) {
  # Predict each replicate
  pred_i <- raster::predict(predictors_raster_c, maxent_model@models[[i]])
  
  # Force as RasterLayer
  if (!inherits(pred_i, "RasterLayer")) stop("Prediction is not RasterLayer")
  
  pred_list[[i]] <- pred_i
}

# Remove NULLs or invalid layers
pred_listc <- Filter(function(x) inherits(x, "RasterLayer"), pred_listc)

# Stack
pred_stack_c <- stack(pred_listc)


suitability_mean_c <- raster::stackApply(pred_stack_c,
                                       indices = rep(1, nlayers(pred_stack_c)),
                                       fun = base::mean)


threshold_val <- 0.5
suitability_binary_c <- suitability_mean_c > threshold_val

par(mfrow=c(1,2))
plot(suitability_mean_c, main="Mean MaxEnt Suitability")
plot(suitability_binary_c, main="Binary Presence/Absence Map")

# Convert raster to dataframe for ggplot ploting
suitability_dfm_c <- as.data.frame(rasterToPoints(suitability_mean_c))
colnames(suitability_dfm_c) <- c("x", "y", "suitability")

# Reproject shapefile and points to match raster
df_ib_c_utm <- st_transform(df_ib_c, crs(suitability_mean_c))
occurrences_c_utm <- st_transform(occurrences_sfc, crs(suitability_mean_c))
absent_sites_c_utm <- st_transform(absent_sites_in_chal, crs(suitability_mean_c))

# Plot again
ggplot() +
  geom_raster(data = suitability_dfm_c, aes(x = x, y = y, fill = suitability)) +
  scale_fill_viridis_c(name = "Suitability") +
  geom_sf(data = df_ib_c_utm, fill = NA, color = "white", size = 0.6) +
  geom_sf(data = occurrences_c_utm, color = "red", size = 2) +
  geom_sf(data = absent_sites_c_utm, color = "blue", size = 2) +
  labs(title = "Habitat Suitability") +
  theme_minimal() +
  coord_sf()


##----------------------------------------------------------------------------------------------------------------##
##Model Evaluation using ROC--------------------------------------------------------------------------------------##
##----------------------------------------------------------------------------------------------------------------##                     
# Use training/testing split for evaluation
presence_valuesc <- raster::extract(suitability_c, occurrences_in_chal)
background_pointsc <- randomPoints(predictors_rc, 500)  # Generate background points
background_valuesc <- raster::extract(suitability_c, background_pointsc)

# Combine predictions
labelsc <- c(rep(1, length(presence_valuesc)), rep(0, length(background_valuesc)))
predictionsc <- c(presence_valuesc, background_valuesc)

# Evaluate model performance using AUC
roc_curvec <- roc(labelsc, predictionsc)
plot(roc_curvec, main = "ROC Curve for MaxEnt Model")
auc_value <- auc(roc_curvec)  # Calculate AUC value
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), bty = "n")


# Convert ROC curve to a data frame
roc_dfc <- data.frame(
  specificity = rev(roc_curvec$specificities),
  sensitivity = rev(roc_curvec$sensitivities)
)

# Extract AUC
auc_valuec <- auc(roc_curvec)

# Plot with ggplot2
chalroc <- ggplot(roc_dfc, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1f78b4", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
  labs(title = "ROC Curve for MaxEnt Model (Challenge)",
       subtitle = paste("AUC =", round(auc_valuec, 3)),
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_manuscript()

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Chalenge ROC curve.pdf'), chalroc, width = 11, height = 10)


                     
# Perform k-fold cross-validation
library(sp)
presence_sp <- SpatialPoints(occurrences[, c("lon", "lat")], 
                             proj4string = CRS("+proj=longlat +datum=WGS84"))

maxent_model_cv <- dismo::maxent(x = predictors_r, p = presence_sp, args="replicates=5")

# Predict across study area
# Create a raster stack to store predictions
pred_stack <- stack()

for(i in 1:length(maxent_model_cv@models)) {
  pred_stack <- stack(pred_stack, predict(predictors_r, maxent_model_cv@models[[i]]))
}

pred_stack_terra <- rast(pred_stack) 

# Average across replicates
suitability_cv <- app(pred_stack_terra, fun = base::mean)

plot(suitability_cv)


##Evaluate with AIC

# Occurrences as matrix or data.frame
occ_matrix <- as.data.frame(occurrences[, c("lon", "lat")])

# Background points (optional, v2 can generate automatically)
bg_points <- NULL  # let ENMeval sample background if not provided

# Run ENMevaluate
eval <- ENMevaluate(
  occ = occ_matrix,          # presence points
  env = predictors_r,        # raster stack
  method = "block",          # spatial partitioning
  algorithm = "maxnet",      # ENMeval v2 uses 'maxnet' instead of maxent.jar
  parallel = TRUE
)

# View results including AICc
eval@results


# Check AIC results
eval@results  # Contains AICc, deltaAICc, and other metrics

# View AICc for each model
eval@results$AICc



##Evaluating with Uncertainties and confidence limits
# install.packages(c("ENMeval","terra","sf","maxnet"))  # maxnet backend is used by ENMeval v2+
library(ENMeval)
library(terra)
library(sf)

# Convert your RasterStack to SpatRaster if needed
predictors_rc <- terra::rast(predictors_rc)   # ensure it's a SpatRaster

# Convert land use to factor
predictors_rc$landuse_code <- as.factor(predictors_rc$landuse_code)

# Make sure occurrences are in matrix form
occs_matc <- sf::st_coordinates(occurrences_in_chal)

# Convert occurrence points to an sf object
occs_sfc <- st_as_sf(data.frame(occs_matc), coords = c("X","Y"), crs = 4326)  # WGS84

# Project to UTM zone 31N (same as your raster)
occs_utmc <- st_transform(occs_sfc, crs(predictors_rc))

# Convert back to matrix for ENMevaluate
occs_cleanc <- st_coordinates(occs_utmc)


# Create a background matrix (all raster cells)
bg_cellsc <- terra::as.data.frame(predictors_rc, xy = TRUE, cells = FALSE)
bg_coordsc <- bg_cellsc[, c("x","y")]

# Rename background columns to match occurrences
colnames(bg_coordsc) <- colnames(occs_matc)

# Extract predictor values for occurrences
# occ_vals <- terra::extract(predictors_r, occs_mat)
# keep_occ <- complete.cases(occ_vals)
# occs_clean <- occs_mat[keep_occ, ]

# Extract predictor values for background
bg_valsc <- terra::extract(predictors_rc, bg_coordsc)
keep_bgc <- complete.cases(bg_valsc)
bg_cleanc <- bg_coordsc[keep_bgc, ]


en_c <- ENMevaluate(
  occs = occs_cleanc,
  env = predictors_rc,
  bg = bg_cleanc,
  algorithm = "maxnet",
  partitions = "randomkfold",
  categoricals = "landuse_code",
  parallel = FALSE,
  tune.args = list(
    fc = c("L","LQ","H","LQH"),
    rm = c(0.5, 1, 2)
  )
)

#Explore the file(en_c)
en_c@results

# Check which model had the best AICc
en_c@results[which.min(en_c@results$AICc), ]

plot(en_c, "jackknife")


# ENMevaluate output object contains predictions for each run.
# Extract raster predictions: en@results or en@predictions depending on version
# For ENMeval v2, use en@predictions (list of rasters)
pred_listc <- en_c@predictions  # may be a list of SpatRaster objects or Raster*; check structure

# # If predictions are returned as a list of SpatRaster or RasterStack per run, convert to terra SpatRaster list
# library(terra)
# 
# # If you want a list with one raster
# rep_rastersc <- list(pred_listc[[1]])

# # Flatten into a list of SpatRasters (one per replicate)
# rep_rasters <- list()
# i <- 1
# for (p in pred_list) {
#   # p might be a RasterLayerStack with single layer or a list; convert safely:
#   pr <- rast(p)            # safe conversion from Raster* or list
#   # If pr has multiple layers (e.g. one per model variant), pick the layer you want.
#   # Here we assume a single layer per replicate.
#   rep_rasters[[i]] <- pr[[1]]
#   i <- i + 1
# }

# Stack them
#s <- rast(rep_rasters)     # SpatRaster with layers = replicates

# Define your own mean function
my_mean <- function(x) {
  mean(x, na.rm = TRUE)   # na.rm = TRUE ignores NA values
}

# Compute mean and SD per cell
mean_rc <- app(pred_listc, fun = my_mean, cores = 1)
sd_rc   <- app(pred_listc, fun = sd,   cores = 1)

# 95% CI by mean ± 1.96*sd (approximate)
ci_lower_approxc <- mean_rc - 1.96 * sd_rc
ci_upper_approxc <- mean_rc + 1.96 * sd_rc

# Define a function that returns the lower 2.5% quantile for a single cell across layers
ci_lower_fun <- function(x) {
  quantile(x, probs = 0.025, na.rm = TRUE)[[1]]  # select the first element
}

ci_upper_fun <- function(x) {
  quantile(x, probs = 0.975, na.rm = TRUE)[[1]]  # select the first element
}

# Better: percentile-based CI (empirical 2.5% and 97.5%)
ci_lowerc <- app(pred_listc, fun = ci_lower_fun)

ci_upperc <- app(pred_listc, fun = ci_upper_fun)

# Save outputs
writeRaster(mean_rc, "maxent_meanc.tif", overwrite = TRUE)
writeRaster(sd_rc,   "maxent_sdc.tif",   overwrite = TRUE)
writeRaster(ci_lowerc, "maxent_CI025c.tif", overwrite = TRUE)
writeRaster(ci_upperc, "maxent_CI975c.tif", overwrite = TRUE)

# Quick plots
plot(mean_r, main = "Mean suitability")
plot(sd_r, main = "SD of suitability (replicates)")
plot(ci_upper, main = "97.5% percentile (upper)")
plot(ci_lower, main = "25.0% percentile (lower)")

library(terra)

# Stack rasters
r_stackc <- c(mean_rc, sd_rc, ci_upperc, ci_lowerc)
names(r_stackc) <- c("Mean", "SD", "Upper 97.5%", "Lower 2.5%")

# Plot all in one figure
plot(r_stackc)


library(terra)
library(ggplot2)
library(gridExtra)

# Convert raster to data.frame
r_to_df <- function(r) as.data.frame(r, xy = TRUE)

df_meanc <- r_to_df(mean_rc)
df_sdc <- r_to_df(sd_rc)
df_upperc <- r_to_df(ci_upperc)
df_lowerc <- r_to_df(ci_lowerc)

# # Make ggplots
# p1 <- ggplot(df_meanc, aes(x=x, y=y, fill=lyr.1)) + geom_raster() + coord_fixed() + labs(title="Mean")
# p2 <- ggplot(df_sdc, aes(x=x, y=y, fill=sd)) + geom_raster() + coord_fixed() + labs(title="SD")
# p3 <- ggplot(df_upperc, aes(x=x, y=y, fill=lyr.1)) + geom_raster() + coord_fixed() + labs(title="97.5%")
# p4 <- ggplot(df_lowerc, aes(x=x, y=y, fill=lyr.1)) + geom_raster() + coord_fixed() + labs(title="2.5%")
# 
# # Arrange plots
# combined_plot <- grid.arrange(p1, p2, p3, p4, ncol=2)
# 
# # Save
# ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/Agugu_Suitability_plots.pdf"), combined_plot, width = 11, height = 10)

library(ggplot2)
library(gridExtra)

# Define a color palette
pal <- c("#f7fcf5","#c7e9c0","#41ab5d","#006d2c")  # light to dark green
pal_sd <- c("#f7f7f7","#cccccc","#969696","#252525") # grey scale for SD
# New color palette: purple to yellow
pal <- c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")  
pal <- c("#1b9e77", "#66c2a5", "#fbb4b9", "#d53e4f")


# Plot Mean suitability
p1c <- ggplot(df_meanc, aes(x=x, y=y, fill=lyr.1)) +
  geom_raster(color="black", size=0.1) +   # black outlines
  scale_fill_gradientn(colors = pal, na.value = "transparent") +
  coord_fixed() +
  labs(title="Mean suitability") +
  theme_minimal() +
  theme(axis.text=element_blank(), axis.title=element_blank())

# Plot SD of suitability
p2c <- ggplot(df_sdc, aes(x=x, y=y, fill= sd)) +
  geom_raster(color="black", size=0.1) +
  scale_fill_gradientn(colors = pal_sd, na.value = "transparent") +
  coord_fixed() +
  labs(title="SD of suitability") +
  theme_minimal() +
  theme(axis.text=element_blank(), axis.title=element_blank())

# Plot 97.5% percentile
p3c <- ggplot(df_upperc, aes(x=x, y=y, fill=lyr.1)) +
  geom_raster(color="black", size=0.1) +
  scale_fill_gradientn(colors = pal, na.value = "transparent") +
  coord_fixed() +
  labs(title="97.5% percentile") +
  theme_minimal() +
  theme(axis.text=element_blank(), axis.title=element_blank())

# Plot 2.5% percentile
p4c <- ggplot(df_lowerc, aes(x=x, y=y, fill=lyr.1)) +
  geom_raster(color="black", size=0.1) +
  scale_fill_gradientn(colors = pal, na.value = "transparent") +
  coord_fixed() +
  labs(title="2.5% percentile") +
  theme_minimal() +
  theme(axis.text=element_blank(), axis.title=element_blank())

# Arrange plots side by side (2x2)
combined_plotc <- grid.arrange(p1c, p2c, p3c, p4c, ncol=2)

# Save to PDF
ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", "Challenge_Suitability_plots.pdf"), combined_plotc, width = 11, height = 10)


