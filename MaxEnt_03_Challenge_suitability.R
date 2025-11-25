Sys.setenv(JAVA_HOME="C:/Program Files/Eclipse Adoptium/jdk-21.0.7.6-hotspot")
library(rJava)  
library(terra)
library(patchwork)
library(dismo)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


#1. Convert SpatRaster to RasterStack for MaxEnt
predictors_rc <- raster::stack(predictors_c_subset)

# 2. Make occurrence points a matrix of coordinates
occ_coordsc <- st_coordinates(occurrences_in_chal)  # gives X = lon, Y = lat


# Convert sf to matrix in raster CRS
occ_coordsc <- st_coordinates(st_transform(occurrences_in_chal, crs(predictors_rc)))
# Get raster extent
terra::ext(predictors_rc)

# Extract values for all layers
occ_valsc <- terra::extract(predictors_c, vect(st_transform(occurrences_in_chal, crs(predictors_c))))

occ_valsc

# Remove rows where any predictor is NA
valid_idxc <- complete.cases(occ_valsc[, -1])  # exclude ID column
occ_coords_validc <- occ_coordsc[valid_idxc, ]

#3. Run MaxEnt
maxent_model_c <- maxent(
  x = predictors_rc, 
  p = occ_coords_validc,
  factors = "landuse_code")

# 4. View summary
print(maxent_model_c)

# ---------------------------
# 4. Jackknife test
# ---------------------------

# Re-run MaxEnt with jackknife enabled
jkc <- maxent(predictors_rc, occ_coords_validc, args = c("jackknife=true"))

# Base R jackknife plot
plot(jkc, type = "jackknife")

# Extract numeric results
jkc_results <- jkc@results
head(jkc_results)

##Explore feature importance direction
response(jkc, var = "avg_rad")

# Convert results into a dataframe
resc <- as.data.frame(jkc_results)

# Extract variable names as they appear in the result names
varsc <- gsub("\\.contribution.*", "", grep("\\.contribution", rownames(resc), value = TRUE))

# Build a tidy dataframe of jackknife results
jkc_df <- lapply(varsc, function(v) {
  data.frame(
    variable = v,
    with_only = as.numeric(resc[paste0("Training.gain.with.only.", v), 1]),
    without   = as.numeric(resc[paste0("Training.gain.without.", v), 1])
  )
}) %>%
  bind_rows()

# Add the "with all variables" value (same for all)
jk_dfc <- jkc_df %>%
  mutate(all_vars = as.numeric(resc["Regularized.training.gain", 1]))

##Rename variable names
jk_dfc <- jk_dfc %>%
  mutate(variable = recode(variable,
                           "EVI.2"  = "EVI(June 2024)",
                           "EVI.3"  = "EVI(July 2024)",
                           "NDMI.1" = "NDMI (May 2024)",
                           "NDMI.2" = "NDMI (June 2024)",
                           "NDMI.3" = "NDMI (July 2024)",
                           "NDWI.1" = "NDWI (May 2024)",
                           "NDWI.2" = "NDWI (June 2024)",
                           "NDWI.3" = "NDWI (July 2024)",
                           "angle_mean" = "Mean Angle",
                           "avg_rad" = "Avg. Radiance from NTL",
                           "distance2water_30arcsec" = "Distance to Water",
                           "gpw_v4_population_density_rev11_2020_1_deg" = "Population Density",
                           "landuse_code" = "Land Use",
                           "log_area" = "Building Area",
                           "lyr.1" = "Land Surface Temperature",
                           "nndist_mean" = "Mean Nearest Neighbor Distance",
                           "shape_mean" = "Mean Shape"))

# Reshape to long format for ggplot
jkc_long <- jk_dfc %>%
  pivot_longer(cols = c("with_only", "without", "all_vars"),
               names_to = "condition",
               values_to = "gain")

# Plot jackknife bar chart
ggplot(jkc_long, aes(x = variable, y = gain, fill = condition)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("with_only" = "#1f78b4",   # blue
                               "without"  = "#e31a1c",   # red
                               "all_vars" = "#33a02c")) + # green
  labs(title = "MaxEnt Jackknife Test of Variable Importance",
       x = "Environmental Variable",
       y = "Training Gain",
       fill = "Condition") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


##Another Plot
# Separate out "all_vars"
jk_sepc <- jkc_long %>% 
  filter(condition != "all_vars")

jkc_all <- jkc_long %>% 
  filter(condition == "all_vars")

# Plot
jk_chal <- ggplot(jk_sepc, aes(x = gain, y = variable, fill = condition)) +
  geom_col(position = "dodge") +  # side-by-side bars
  geom_col(data = jk_all, aes(x = gain, y = variable), 
           fill = "grey", width = 0.5, inherit.aes = FALSE) +
  labs(title = "MaxEnt Jackknife Test of Variable Importance",
       x = "Training Gain", y = "Environmental Variable") +
  scale_fill_manual(values = c("with_only" = "lightgreen",
                               "without" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Jack n Knife Test for Challenge.pdf'), jk_chal , width = 11, height = 10)



##Replot for Manuscript
# Reorder variables by 'with_only' gain
jk_sepc <- jk_sepc %>% filter(condition !="all_vars") %>%  
  group_by(variable) %>%
  mutate(max_with_only = ifelse(condition == "with_only", gain, NA)) %>%
  ungroup() %>%
  group_by(variable) %>%
  fill(max_with_only, .direction = "downup") %>%
  ungroup() %>%
  mutate(variable = fct_reorder(variable, max_with_only, .desc = TRUE))


# Then plot
jk_chal <- ggplot(jk_sepc, aes(x = gain, y = variable, fill = condition)) +
  geom_col(position = "dodge") +
  #geom_col(data = jk, aes(x = gain, y = variable), 
           #fill = "aliceblue", width = 0.5, inherit.aes = FALSE) +
  labs(title = "MaxEnt Jackknife Test of Variable Importance (Challenge)",
       x = "Training Gain", y = "Environmental Variable") +
  scale_fill_manual(values = c("with_only" = "lightgreen",
                               "without" = "plum",
                               "all_vars")) +
  theme_manuscript() +
  theme(legend.position = "bottom")

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Jack n Knife Test for challenge.pdf'), jk_chal , width = 8, height = 11)


##Evaluate Variable Importance direction 
# Open PDF device
pdf("Challenge jackknife_responses.pdf", width = 11, height = 16)

response(jkc, var = "avg_rad")

response(jkc, var = "NDWI.1")

response(jkc, var = "angle_mean")

# Close PDF device
dev.off()

#Open PDF automatically
shell.exec("Agugu jackknife_responses.pdf")



###-----------------------------------------------------------------------------
# Predict across study area
##------------------------------------------------------------------------------

suitability_c <- predict(predictors_rc, maxent_model_c)

crs(suitability_c) <- "EPSG:32631"

# Plot results
plot(suitability_c, main="Habitat Suitability")
points(occurrences_in_chal, col="red", pch=20)
plot(ward_vectc_mi, border="blue", add=TRUE)


##Plot using ggplot2

# Convert raster to dataframe for ggplot
suitability_dfc <- as.data.frame(rasterToPoints(suitability_c))
colnames(suitability_dfc) <- c("x", "y", "suitability")

# Reproject shapefile and points to match raster
df_ib_c_utm <- st_transform(df_ib_c, crs(suitability_c))
occurrences_utmc <- st_transform(occurrences_in_chal, crs(suitability_c))

# Plot again
ggplot() +
  geom_raster(data = suitability_dfc, aes(x = x, y = y, fill = suitability)) +
  scale_fill_viridis_c(name = "Suitability") +
  geom_sf(data = df_ib_c_utm, fill = NA, color = "white", size = 0.6) +
  geom_sf(data = occurrences_utmc, color = "red", size = 2) +
  #geom_sf(data = absent_sites_utm, color = "blue", size = 2) +
  labs(title = "Habitat Suitability") +
  theme_minimal() +
  coord_sf()


##Plot only within ward shape file
# Make sure raster and ward have the same CRS
raster_crsc <- st_crs(suitability_c)

ward_vect_utmc <- st_transform(df_ib_c, crs = raster_crsc)
library(terra)

# Convert RasterLayer to SpatRaster
suitabilityc_terra <- rast(suitability_c)

# Mask with SpatVector
suitability_maskc <- mask(suitabilityc_terra, ward_vect_utmc)

# Convert to dataframe for ggplot
suit_dfc <- as.data.frame(suitability_maskc, xy = TRUE)
names(suit_dfc)[3] <- "Suitability"

# Plot clipped map
habsuitc <- ggplot() +
  geom_raster(data = suit_dfc, aes(x = x, y = y, fill = Suitability)) +
  scale_fill_viridis_c(option = "plasma") +
  geom_sf(data = st_as_sf(ward_vect_utmc), fill = NA, color = "black", size = 0.5) +
  #geom_sf(data = occurrences_utmc, color = "green", size = 3.5) +
  #geom_sf(data = hh_pos_a_df, color = "black", size = 1) +
  theme_manuscript() +
  labs(title = "Habitat Suitability Challenge", fill = "Suitability")

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Habitat Suitability for Challenge.pdf'), habsuitc, width = 11, height = 10)


#Convert suitability plot to categories

# Load raster
suitability_rastc <- suitability_c

# Create categories
library(classInt)
fisher_breaks <- classIntervals(suit_dfc$Suitability, n = 4, style = "fisher")$brks
fisher_breaks

mc <- matrix(c(-Inf, fisher_breaks[2], 0,
              fisher_breaks[2], fisher_breaks[3], 1,
              fisher_breaks[3], fisher_breaks[4], 2,
              fisher_breaks[4], Inf, 3),
            ncol = 3, byrow = TRUE)


suitability_catc <- reclassify(suitability_rastc, mc)


# # Define colors for categories
# cat_colors <- c("lightgreen","lightblue", "yellow", "red")  # very low, low, medium, high
# 
# # Plot raster with specified colors
# plot(suitability_catc, 
#      col = cat_colors,
#      legend = TRUE,
#      main = "Suitability Categories (0 = Very Low,1 = Low, 2 = Medium, 3 = High)")

#Replot raster with ggplot
# 1. Mask raster to shapefile extent
suitability_maskedc <- mask(suitability_catc, df_ib_c_utm)

# 2. Convert raster to dataframe
suit_dfcatc <- as.data.frame(suitability_maskedc, xy = TRUE)
colnames(suit_dfcatc)[3] <- "class"

# 3. Attach class labels (0–3 → categorical names)
suit_dfcatc$class <- factor(suit_dfcatc$class,
                        levels = 0:3,
                        labels = c("Very Low", "Low", "Medium", "High"))

suit_dfcatc <- suit_dfcatc %>%
  filter(!is.na(class))

# 4. Plot with ggplot2
habsuitcatc <- ggplot() +
  geom_raster(data = suit_dfcatc, aes(x = x, y = y, fill = class)) +
  geom_sf(data = st_as_sf(ward_vect_utmc), fill = NA, color = "black", size = 0.5)+
  coord_sf() +
  scale_fill_manual(values = c(
    "Very Low" = "#d9f0a3",  # light green
    "Low" = "#addd8e",       # medium green
    "Medium" = "#fee08b",    # yellow
    "High" = "#fc8d59"       # orange-red
  )) +
  labs(title = "Suitability Categories",
       fill = "Suitability") +
  theme_manuscript()

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Habitat Suitability Categories for Challenge.pdf'), habsuitcatc, width = 11, height = 10)



#Boxplot showing suitability
library(raster)
library(dplyr)

# Extract values from the raster
suit_valuesc <- raster::extract(suitability_rastc, 1:ncell(suitability_rastc))
suit_valuesc <- data.frame(value = suit_valuesc)

# # Reclassify numeric suitability into categories using your matrix 'mc'
# suit_valuesc$class <- cut(
#   suit_valuesc$value,
#   breaks = c(-Inf, fisher_breaks[2], fisher_breaks[3], fisher_breaks[4], Inf),
#   labels = c("Very Low", "Low", "Moderate", "High"),
#   include.lowest = TRUE
# )


# Remove NAs first
suit_valuesc_clean <- suit_valuesc %>%
  filter(!is.na(value))

# Assuming your dataframe is called suit_values
suit_values_clnc <- suit_valuesc %>%
  filter(!is.na(value) & !is.na(class))


# Quantile values
q25 <- 0.1403995
q50 <- 0.3497512
q75 <- 0.6006984
q80 <- 0.6557162
q90 <- 0.7909244

# Reclassify using these exact quantile values
suit_values_clnc$class <- cut(
  suit_values_clnc$value,
  breaks = c(-Inf, q25, q50, q75, q80, q90, Inf),
  labels = c("Very Low", "Low", "Moderate", "High", "Very High", "Extreme"),
  include.lowest = TRUE
)

# Check class distribution
table(suit_values_clnc$class)



ggplot(suit_values_clnc, aes(x = class, y = value, fill = class)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # alpha for transparency, hide default outliers
  #geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "black") +  # add jittered points
    labs(
    x = "Suitability Class",
    y = "Suitability Value",
    title = "Distribution of Suitability Values by Class"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")



ggsave(paste0(LuDir, '/plots/', Sys.Date(), 'Box plot of larval densitities by settlement n season.pdf'), ld_bxp, width = 10, height = 11)



##Estimate Confidence Intervals
#============================#
#      Required Libraries    #
#============================#
library(raster)
library(dplyr)
library(ggplot2)
library(boot)
library(tidyr)
library(RColorBrewer)

#============================#
#   1. Load Raster Data      #
#============================#
suitability_rast <- raster("suitability_raster.tif")

# Extract raster values into a dataframe
suit_values <- data.frame(
  value = getValues(suitability_rast)
)

# Remove NA values
suit_values <- suit_values %>%
  filter(!is.na(value))

#============================#
#  2. Classify by Percentile #
#============================#
# 75th percentile cutoff for "High" class
#q25 <- quantile(suit_values_clnc$value, 0.25)
#q50 <- quantile(suit_values_clnc$value, 0.50)
#q75 <- quantile(suit_values_clnc$value, 0.75)

q25 <- 0.1403995
q50 <- 0.3497512
q75 <- 0.6006984
q80 <- 0.6557162
q90 <- 0.7909244

suit_values_clnc$class <- cut(
  suit_values_clnc$value,
  breaks = c(-Inf, q25, q50, q75, Inf),
  labels = c("Very Low", "Low", "Moderate", "High"),
  include.lowest = TRUE
)

#============================#
#  3. Bootstrap for CI       #
#============================#
# Function to compute mean
boot_mean <- function(data, indices) {
  d <- data[indices]
  return(mean(d, na.rm = TRUE))
}

# Initialize dataframe for storing CI results
ci_summary <- data.frame(
  class = character(),
  mean_value = numeric(),
  lower_ci = numeric(),
  upper_ci = numeric(),
  stringsAsFactors = FALSE
)

# Compute bootstrap CI per class
set.seed(123)  # for reproducibility
classes <- unique(suit_values_clnc$class)

for (cls in classes) {
  class_data <- suit_values_clnc$value[suit_values_clnc$class == cls]
  
  # Bootstrap 1000 resamples
  boot_res <- boot(class_data, boot_mean, R = 1000)
  ci <- boot.ci(boot_res, type = "perc")
  
  ci_summary <- rbind(ci_summary, data.frame(
    class = cls,
    mean_value = mean(class_data),
    lower_ci = ci$percent[4],  # 2.5% percentile
    upper_ci = ci$percent[5]   # 97.5% percentile
  ))
}

#============================#
#  4. Visualization          #
#============================#
# Ensure the class column in both datasets is a factor with the same levels and order
levels_order <- c("Very Low", "Low", "Moderate", "High", "Very High", "Extreme")  # include all your levels

suit_values_clnc$class <- factor(suit_values_clnc$class, levels = levels_order)
ci_summary$class <- factor(ci_summary$class, levels = levels_order)

# Make plot
ggplot(suit_values_clnc, aes(x = class, y = value, fill = class)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  #geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "black") +
  geom_errorbar(
    data = ci_summary,
    aes(x = class, ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    color = "red",
    size = 1,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Suitability Class",
    y = "Suitability Value",
    title = "Habitat Suitability Distribution with 95% Confidence Intervals"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


##Check again
# Categorize suitability
suit_summary <- suit_values_clnc %>%
  mutate(
    suitability_category = ifelse(value > 0.79, "Highly suitable", "Not highly suitable")
  ) %>%
  summarise(
    total = n(),
    highly_suitable = sum(suitability_category == "Highly suitable"),
    proportion_highly_suitable = highly_suitable / total
  )

suit_summary

