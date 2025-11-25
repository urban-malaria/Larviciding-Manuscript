library(dismo)  # For MaxEnt modeling
library(raster) # For working with spatial raster data
library(maps)   # For map visualization
library(terra)
library(sf)
library(dismo)
library(randomForest)
library(gbm)
library(caret)
library(usdm)
source("functions.R")

# --- Load rasters by type ---
##EVI
evi_dir <- "C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/HLS30m/EVI/Ibadan"

# EVI files
# List all EVI files
evi_files <- list.files(evi_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract month and year from filename
evi_info <- data.frame(
  file = evi_files,
  year = as.numeric(str_extract(evi_files, "(?<=month\\d{2}_)\\d{4}")),
  month = as.numeric(str_extract(evi_files, "(?<=month)\\d{2}(?=_)"))
)

# Filter to May–Aug 2024
evi_study <- evi_info %>% filter(year == 2024, month %in% 5:7)

# Stack the filtered rasters
evi_stack <- stack(evi_study$file)

# Check
print(evi_stack)
plot(evi_stack)


# NDWI files
ndWi_dir <- "C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/field_study_ndwi_30m"

# List all NDWI files in that directory
ndWi_files <- list.files(ndWi_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract month and year from filename 
ndWi_info <- data.frame(
  file = ndWi_files,
  year  = as.numeric(sapply(strsplit(basename(ndWi_files), "_"), `[`, 5)),  # 5th element is year
  month = as.numeric(sub("\\.tif$", "", sapply(strsplit(basename(ndWi_files), "_"), `[`, 6)))  # 6th element is month
)

# Filter to the months of interest (e.g., May–Aug 2024)
ndWi_study <- ndWi_info %>% filter(year == 2024, month %in% 5:7)

# Stack only the filtered rasters
ndWi_stack <- stack(ndWi_study$file)

# Check
print(ndWi_stack)
plot(ndWi_stack)


# Directory where NDMI rasters are stored
ndmi_dir <- "C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/field_study_ndmi_30m"

# List all NDMI files in that directory
ndmi_files <- list.files(ndmi_dir, pattern = "\\.tif$", full.names = TRUE)

# Extract month and year from filename (adjust regex based on your naming convention)
ndmi_info <- data.frame(
  file = ndmi_files,
  year  = as.numeric(sapply(strsplit(basename(ndmi_files), "_"), `[`, 5)),  # 5th element is year
  month = as.numeric(sub("\\.tif$", "", sapply(strsplit(basename(ndmi_files), "_"), `[`, 6)))  # 6th element is month
)


# Filter to the months of interest (e.g., May–Aug 2024)
ndmi_study <- ndmi_info %>% filter(year == 2024, month %in% 5:7)

# Stack only the filtered rasters
ndmi_stack <- stack(ndmi_study$file)

# Check
print(ndmi_stack)
plot(ndmi_stack)

# --- Combine all into one stack ---
#------------------------------------------------------------------------------
# Convert RasterStack to SpatRaster
##EVI
evi_stack_terrac <- rast(evi_stack)
crs(evi_stack_terrac) <- "EPSG:32631"

# Make sure polygon is in the same CRS
#df_ib_c_proj <- st_transform(df_ib_c, crs(evi_stack_terra))
df_ib_c_geom <- df_ib_c_proj["geometry"]
ward_vectc <- vect(df_ib_c_geom)

# Reproject ward polygon to UTM 31N
ward_vectc_utm <- project(ward_vectc, evi_stack_terrac)

# Crop raster to polygon extent
evi_stack_cropc <- crop(evi_stack_terrac, ward_vectc_utm)

# Mask raster to polygon boundary
evi_stack_maskc <- mask(evi_stack_cropc, ward_vectc_utm)

# Plot to check
plot(evi_stack_maskc, main="Cropped + Masked EVI (May-July 2024)")
lines(ward_vectc_utm, col="blue")


# Crop & mask
#evi_stack_cropc <- crop(evi_stack_terrac, ward_vectc)
#evi_stack_maskc <- mask(evi_stack_cropc, ward_vectc)

# Plot
plot(evi_stack_maskc[[1]])
plot(ward_vectc, add=TRUE)

#--------------------------------------------------------------------------------
#NDWI
# Convert RasterStack to SpatRaster
ndWi_stack_terra <- rast(ndWi_stack)

# Make sure polygon is in the same CRS
df_ib_c_projWi <- st_transform(df_ib_c, crs(ndWi_stack_terra))
df_ib_c_geomWi <- df_ib_c_projWi["geometry"]
ward_vectc_Wi <- vect(df_ib_c_geomWi)

# Crop & mask
ndWi_stack_cropc <- crop(ndWi_stack_terra, ward_vectc_Wi)
ndWi_stack_maskc <- mask(ndWi_stack_terra, ward_vectc_Wi)

# Plot
plot(ndWi_stack_cropc[[1]])
plot(ward_vectc_Wi, add=TRUE)

#-----------------------------------------------------------
#NDMI
# Convert RasterStack to SpatRaster
ndmi_stack_terra <- rast(ndmi_stack)

# Make sure polygon is in the same CRS
df_ib_c_projmi <- st_transform(df_ib_c, crs(ndmi_stack_terra))
df_ib_c_geommi <- df_ib_c_projmi["geometry"]
ward_vectc_mi <- vect(df_ib_c_geommi)

# Now crop & mask
ndmi_stack_cropc <- crop(ndmi_stack_terra, ward_vectc_mi)
ndmi_stack_maskc <- mask(ndmi_stack_terra, ward_vectc_mi)

# Plot
plot(ndmi_stack_cropc[[1]])
plot(ward_vectc_mi, add=TRUE)

#-------------------------------------------------------------------------------
##NTL---------------------------------------------------------------------------
#-------------------------------------------------------------------------------

ntl_rast <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/night_timel_lights/VIIRS_NTL_2024_Nigeria.tif")

plot(ntl_rast)

# Crop raster to polygon extent
ntl_cropc <- crop(ntl_rast, df_ib_c)

# Convert RasterStack to SpatRaster
ntl_stack_terra <- rast(ntl_rast)

# Make sure polygon is in the same CRS
df_ib_c_projntl <- st_transform(df_ib_c, crs(ntl_stack_terra))
df_ib_c_geomntl <- df_ib_c_projntl["geometry"]
ward_vectc_ntl <- vect(df_ib_c_geomntl)

# Crop & mask
ntl_stack_cropc <- crop(ntl_stack_terra, ward_vectc_ntl)
#ntl_stack_maskc <- mask(ntl_stack_terra, ward_vectc_ntl)

# Plot
plot(ntl_cropc[[1]])
plot(ward_vectc_ntl, add=TRUE)

#-------------------------------------------------------------------------------
##Population Density
#-------------------------------------------------------------------------------
popn_den <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/Population/NGA_pop_density/gpw_v4_population_density_rev11_2020_1_deg.tif") 

# Transform CRS to match raster
df_ib_c_proj <- st_transform(df_ib_c, crs(popn_den))

# Convert to SpatVector
ward_vectc <- vect(df_ib_c_proj["geometry"])
pop_vals <- terra::extract(popn_den, ward_vectc)

# Check extracted values
head(pop_vals)
summary(pop_vals)

#Crop & mask
popnd_stack_cropc <- crop(popn_den, ward_vectc)

popnd_stack_maskc <- mask(popn_den, ward_vectc)

# Plot
ext(ward_vectc)
plot(popnd_stack_maskc)
plot(ward_vectc, add=TRUE)

##--------------------------------------------------------------
# pop_crop <- crop(popn_den, ext(ward_vectc_popnd))
# pop_mask <- mask(pop_crop, ward_vectc_popnd)
# 
# plot(pop_mask, main="Population raster for ward")
# plot(ward_vectc_popnd, add=TRUE)
# 
# plot(pop_mask)
# plot(ward_vectc_popnd, add=TRUE)
# zoom(ext(ward_vectc_popnd))

##-------------------------------------------------------------------------
##Distance water bodies
dwb <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/distance_to_water_bodies/distance2water_30arcsec.tif")

dwb_cropc <- crop(dwb, df_ib_c)
plot(dwb_cropc)
chal_vect <- vect(df_ib_c)  # convert sf → SpatVector
plot(chal_vect, add = TRUE, border = "red", lwd = 2)



##------------------------------------------------------------------------------
#Land Use
##------------------------------------------------------------------------------
bfp <- st_read("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/building_footprints/nigeria_footprints/nigeria blocks 2/Nigeria_Blocks_V1.shp")

bfp_oyo <- bfp %>% 
  dplyr::filter(state == "Oyo")

bfp_oyo <- sf::st_zm(bfp_oyo, drop = TRUE, what = "ZM")

# Check validity
sum(!st_is_valid(bfp_oyo))

# Fix invalid geometries
bfp_oyo <- st_make_valid(bfp_oyo)

#Plot to see extent
ggplot(bfp_oyo) +
  geom_sf(aes(fill = landuse), color = NA) +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal() +
  labs(title = "Landuse in Oyo State", fill = "Land Use")

##Filter to Challenge
st_crs(bfp_oyo) <- st_crs(df_ib_c)

chal_bfp <- st_intersection(bfp_oyo, st_union(df_ib_c))

#Plot to see extent
ggplot(chal_bfp) +
  geom_sf(aes(fill = landuse), color = NA) +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal() +
  labs(title = "Landuse in Challenge", fill = "Land Use")

#Convert to raster to ensure ease of stacking
# Convert sf → SpatVector
chal_vect <- vect(chal_bfp)

# Encode landuse categories as numeric codes
chal_vect$landuse_code <- as.numeric(as.factor(chal_vect$landuse))

# Ensure rasters are in the same extent/res) using EVI
templatec <- evi_stack_cropc[[1]]  # first layer of your EVI stack

# Make sure chal_vect is in the same CRS as template
chal_vect_utm <- project(chal_vect, crs(templatec))

# Crop template to chal_vect extent
template_cropc <- crop(templatec, ext(chal_vect_utm))

# Rasterize landuse using numeric codes
landusec_raster <- rasterize(chal_vect_utm, templatec, field = "landuse_code")

plot(landusec_raster)

landuse_key <- data.frame(
  code = as.numeric(as.factor(chal_bfp$landuse)),
  category = levels(as.factor(chal_bfp$landuse))
)
print(landuse_key)

##--------------------------------------------------------------------------
#LST
##--------------------------------------------------------------------------
lst_cropc <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/kano_ibadan/kano_ibadan_ento/challenge_lstraster.tif")


##------------------------------------------------------------------------------
##Building morphology
##------------------------------------------------------------------------------
chal_bfp <- read.csv("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/building_morphology/challenge_agugu/cbound_vectors/Challenge/Challenge_cbound_with_coords.csv")

chal_bf_stack


--------------------------------------------------------------------------------
# 1. Make sure all rasters have the same CRS (use EVI as template)
ndWi_stack_cropc <- project(ndWi_stack_cropc, crs(evi_stack_cropc))
ndmi_stack_cropc <- project(ndmi_stack_cropc, crs(evi_stack_cropc))
ntl_stack_cropc <- project(ntl_cropc, crs(evi_stack_cropc))
popdn_stack_cropc <- project(popnd_stack_cropc, crs(evi_stack_cropc))
dwb_stack_cropc <- project(dwb_cropc, crs(evi_stack_cropc))
lst_stack_cropc <- project(lst_cropc, crs(evi_stack_cropc))
landuse_stack_cropc <-project(landusec_raster, crs(evi_stack_cropc))
chal_bf_stack_cropc <- project(chal_bf_stack, crs(evi_stack_cropc))

# 2. Align extents and resolution
ndWi_stack_resc <- resample(ndWi_stack_cropc, evi_stack_cropc, method="bilinear")
ndmi_stack_resc <- resample(ndmi_stack_cropc, evi_stack_cropc, method="bilinear")
ntl_stack_resc <- resample(ntl_stack_cropc, evi_stack_cropc, method="bilinear")
popdn_stack_resc <- resample(popdn_stack_cropc, evi_stack_cropc, method="bilinear")
dwb_stack_resc <- resample(dwb_stack_cropc, evi_stack_cropc, method="bilinear")
lst_stack_resc <- resample(lst_stack_cropc, evi_stack_cropc, method="bilinear")
landuse_stack_resc <- resample(landuse_stack_cropc, evi_stack_cropc, method="bilinear")
chal_bf_stack_resc <- resample(chal_bf_stack_cropc, evi_stack_cropc, method="bilinear")

##Stack
predictors_c <- c(evi_stack_cropc, ndWi_stack_resc, ndmi_stack_resc, ntl_stack_resc, popdn_stack_resc,
                  dwb_stack_resc, lst_stack_resc, landuse_stack_resc, chal_bf_stack_resc)

# Check results
print(predictors_c)
plot(predictors_c)



##Explore Correlation
# Sample values across your study area
valsc <- values(predictors_c)

# Remove rows with NA
valsc <- na.omit(valsc)

# Compute correlation matrix
cor_matc <- cor(valsc)

# Find highly correlated pairs (r > 0.7 or r < -0.7)
high_corc <- which(abs(cor_matc) > 0.7 & abs(cor_matc) < 1, arr.ind = TRUE)

high_corc

library(ggcorrplot)

chalcorr <- ggcorrplot(cor_matc, 
           hc.order = TRUE,       # hierarchical clustering
           type = "upper",
           lab = TRUE,            # show correlation values
           lab_size = 3,
           colors = c("blue", "white", "red"), 
           title = "Correlation Matrix of Environmental Variables")

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Correlation Matrix for Challenge.pdf'), chalcorr, width = 11, height = 10)


##VIF as correlation check

vif_resultc <- vifstep(predictors_c, th = 5)
vif_resultc

# List selected layers (after dropping collinear ones)
vif_resultc@selected

# names to keep
vars_keep <- c("EVI.2", "EVI.3",
               "NDWI.1", "NDWI.2", "NDWI.3",
               "NDMI.1", "NDMI.2", "NDMI.3",
               "avg_rad", "lyr.1", "landuse_code", 
               "gpw_v4_population_density_rev11_2020_1_deg",
               "distance2water_30arcsec", "nndist_mean", "log_area", "angle_mean", "shape_mean")

# subset the SpatRaster
predictors_c_subset <- predictors_c[[vars_keep]]

# check result
predictors_c_subset

# Save as multi-layer GeoTIFF
writeRaster(predictors_c_subset, 
            filename = "predictors_c_subset.tif", 
            overwrite = TRUE)

predictors_c_subset <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/kano_ibadan/kano_ibadan_ento/predictors_c_subset.tif")



# # Drop layers by name
# predictors_selectedc <- predictors_c[[ !names(predictors_c) %in% c("EVI.1", "EVI.2", "NDMI.1", 
#                                                                   "gpw_v4_population_density_rev11_2020_1_deg") ]]
# 
# plot(env_stack_selected)


