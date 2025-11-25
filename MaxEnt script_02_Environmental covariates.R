library(raster)
library(terra)
library(dismo)
library(terra)
library(sp)
library(caret)
        
source("functions.R")


# --- Load rasters by type ---

# EVI files
evi_dir <- "C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/HLS30m/EVI/Ibadan"


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


# # NDVI files
# ndvi_files <- list.files(pattern = "^NDVI_.*\\.tif$", full.names = TRUE)
# ndvi_stack <- stack(ndvi_files)

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

# Extract month and year from filename 
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


#------------------------------------------------------------------------------
#Night time Lights
#------------------------------------------------------------------------------

ntl_rast <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/night_timel_lights/VIIRS_NTL_2024_Nigeria.tif")

plot(ntl_rast)

##convert raster
# ntl_stack_terra <- rast(evi_stack)
# crs(evi_stack_terra) <- "EPSG:32631"

# Crop raster to polygon extent
ntl_crop <- crop(ntl_rast, df_ib_a)

plot(ntl_crop)
ag_vect <- vect(df_ib_a)  # convert sf → SpatVector
plot(ag_vect, add = TRUE, border = "red", lwd = 2)


#--------------------------------------------------------------------------
##Soil Wetness
#--------------------------------------------------------------------------
# soil_moist <- rast("C:/Users/ebamgboye/Downloads/GIOVANNI-g4.timeAvgOverlayMap.GLDAS_NOAH10_M_2_1_SoilMoi100_200cm_inst.20240501-20250930.2E_3N_15E_15N.tif")
# soil_moist_crop <- crop(soil_moist, df_ib_a)
# plot(soil_moist_crop)
# plot(ag_vect, add = TRUE, border = "red", lwd = 2)

##Population density
popn_den <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/NGA_pop_density/gpw_v4_population_density_rev11_2020_1_deg.tif") 

# Crop raster to polygon extent
popn_den_crop <- crop(popn_den, df_ib_a)
plot(popn_den_crop)
ag_vect <- vect(df_ib_a)  # convert sf → SpatVector
plot(ag_vect, add = TRUE, border = "red", lwd = 2)

#------------------------------------------------------------------------------
##Distance to Water Bodies
#------------------------------------------------------------------------------
dwb <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/Raster_files/distance_to_water_bodies/distance2water_30arcsec.tif")

# Crop raster to polygon extent
dwb_crop <- crop(dwb, df_ib_a)
plot(dwb_crop)
ag_vect <- vect(df_ib_a)  # convert sf → SpatVector
plot(ag_vect, add = TRUE, border = "red", lwd = 2)


#LST
agugu_lst <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/kano_ibadan/kano_ibadan_ento/agugu_lstraster.tif")
plot(agugu_lst)
plot(ag_vect, add = TRUE, border = "red", lwd = 2)


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

ggplot(bfp_oyo) +
  geom_sf(aes(fill = landuse), color = NA) +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal() +
  labs(title = "Landuse in Oyo State", fill = "Land Use")

##Filter to Agugu
st_crs(bfp_oyo) <- st_crs(df_ib_a)

agu_bfp <- st_intersection(bfp_oyo, st_union(df_ib_a))

ggplot(agu_bfp) +
  geom_sf(aes(fill = landuse), color = NA) +
  scale_fill_viridis_d(option = "plasma") +
  theme_minimal() +
  labs(title = "Landuse in Agugu", fill = "Land Use")

#Convert to raster to ensure ease of stacking
# Convert sf → SpatVector
agu_vect <- vect(agu_bfp)

# Encode landuse categories as numeric codes
agu_vect$landuse_code <- as.numeric(as.factor(agu_vect$landuse))

# Use one of your predictor rasters as template (ensures same extent/res)
templatea <- evi_stack_crop[[1]]  # first layer of your EVI stack

# Make sure agu_vect is in the same CRS as your template
agu_vect_utm <- project(agu_vect, crs(templatea))

# Crop template to chal_vect extent
template_cropa <- crop(templatea, ext(agu_vect_utm))

# Rasterize landuse using numeric codes
landusea_raster <- rasterize(agu_vect_utm, templatea, field = "landuse_code")

plot(landusea_raster)

# landuse_key <- data.frame(
#   code = as.numeric(as.factor(chal_bfp$landuse)),
#   category = levels(as.factor(chal_bfp$landuse))
# )
# print(landuse_key)

#-------------------------------------------------------------------------------
##Make sure all rasters are aligned and stack
##------------------------------------------------------------------------------
#EVI
# Convert RasterStack to SpatRaster
evi_stack_terra <- rast(evi_stack)

crs(evi_stack_terra) <- "EPSG:32631"

#Make sure polygon is in the same CRS
df_ib_a_proj <- st_transform(df_ib_a, crs(evi_stack_terra))
df_ib_a_geom <- df_ib_a_proj["geometry"]
ward_vect <- vect(df_ib_a_geom)

ward_vect_proj <- terra::project(ward_vect, evi_stack_terra)
evi_stack_crop <- terra::crop(evi_stack_terra, ward_vect_proj)
evi_stack_mask <- terra::mask(evi_stack_crop, ward_vect_proj)

# Now crop & mask
evi_stack_crop <- crop(evi_stack_terra, ward_vect_proj)
evi_stack_mask <- mask(evi_stack_crop, ward_vect_proj)

# Plot
plot(evi_stack_mask[[1]])
plot(ward_vect_proj, add=TRUE)

#--------------------------------------------------------------------------------
#NDWI
# Convert RasterStack to SpatRaster
ndWi_stack_terra <- rast(ndWi_stack)

# Make sure polygon is in the same CRS
df_ib_a_projWi <- st_transform(df_ib_a, crs(ndWi_stack_terra))
df_ib_a_geomWi <- df_ib_a_projWi["geometry"]
ward_vect_Wi <- vect(df_ib_a_geomWi)

# Now crop & mask
ndWi_stack_crop <- crop(ndWi_stack_terra, ward_vect_Wi)
ndWi_stack_mask <- mask(ndWi_stack_terra, ward_vect_Wi)

# Plot
plot(ndWi_stack_crop[[1]])
plot(ward_vect_Wi, add=TRUE)

#-----------------------------------------------------------
#NDMI
# Convert RasterStack to SpatRaster
ndmi_stack_terra <- rast(ndmi_stack)

# Make sure polygon is in the same CRS
df_ib_a_projmi <- st_transform(df_ib_a, crs(ndmi_stack_terra))
df_ib_a_geommi <- df_ib_a_projmi["geometry"]
ward_vect_mi <- vect(df_ib_a_geommi)

# Now crop & mask
ndmi_stack_crop <- crop(ndmi_stack_terra, ward_vect_mi)
ndmi_stack_mask <- mask(ndmi_stack_terra, ward_vect_mi)

# Plot
plot(ndmi_stack_crop[[1]])
plot(ward_vect_mi, add=TRUE)




##Stack Environmental Covariates

# 1. Make sure all rasters have the same CRS
ndWi_stack_crop <- project(ndWi_stack_crop, crs(evi_stack_crop))
ndmi_stack_crop <- project(ndmi_stack_crop, crs(evi_stack_crop))
lst_stack_crop  <- project(agugu_lst,  crs(evi_stack_crop))
ntl_stack_crop  <- project(ntl_crop,  crs(evi_stack_crop))
popn_den_crop  <-  project(popn_den_crop,  crs(evi_stack_crop))
dwb_crop       <-  project(dwb_crop,  crs(evi_stack_crop))
landusea_stack_crop  <- project(landusea_raster,  crs(evi_stack_crop))


# 2. Align extents and resolution
ndWi_stack_res <- resample(ndWi_stack_crop, evi_stack_crop, method="bilinear")
ndmi_stack_res <- resample(ndmi_stack_crop, evi_stack_crop, method="bilinear")
lst_stack_res  <- resample(lst_stack_crop,  evi_stack_crop, method="bilinear")
ntl_stack_res  <- resample(ntl_stack_crop,  evi_stack_crop, method="bilinear")
popn_den_res  <- resample(popn_den_crop,  evi_stack_crop, method="bilinear")
dwb_res  <- resample(dwb_crop,  evi_stack_crop, method="bilinear")
landusea_stack_res  <- resample(landusea_stack_crop,  evi_stack_crop, method="bilinear")


# 3. Combine into one multilayer stack for analysis
all_env_vars <- c(evi_stack_crop, ndWi_stack_res, ndmi_stack_res,
                  ntl_stack_res, popn_den_res, dwb_res, lst_stack_res, landusea_stack_res)

# Convert to stack
env_stack <- stack(all_env_vars)

plot(all_env_vars)


##Explore Correlation
# Sample values across your study area
vals <- getValues(env_stack)

# Remove rows with NA
vals <- na.omit(vals)

# Compute correlation matrix
cor_mat <- cor(vals)

# Find highly correlated pairs (r > 0.7 or r < -0.7)
high_cor <- which(abs(cor_mat) > 0.7 & abs(cor_mat) < 1, arr.ind = TRUE)

high_cor


aggcorr <- ggcorrplot(cor_mat, 
           hc.order = TRUE,       # hierarchical clustering
           type = "upper",
           lab = TRUE,            # show correlation values
           lab_size = 3,
           colors = c("blue", "white", "red"), 
           title = "Correlation Matrix of Environmental Variables")

ggsave(paste0(LuDir, '/plots/', Sys.Date(), "/", 'Correlation Matrix for Agugu.pdf'), aggcorr, width = 11, height = 10)


##Use VIF to cross check
vif_result <- vifstep(all_env_vars, th = 5)
vif_result

# List selected layers (after dropping collinear ones)
vif_result@selected

# names to keep
vars_keep <- c("EVI.1","EVI.2", "EVI.3",
               "NDWI.1", "NDWI.2", "NDWI.3",
                "NDMI.2", "NDMI.3",
               "avg_rad",
               "distance2water_30arcsec")

##
# subset the SpatRaster
predictors_a_subset <- all_env_vars[[vars_keep]]

# check result
predictors_a_subset

# Save as multi-layer GeoTIFF
writeRaster(predictors_a_subset, 
            filename = "predictors_a_subset.tif", 
            overwrite = TRUE)

# Example: keep only selected layers
env_stack_selected <- dropLayer(env_stack, c("NDMI.1", "NDWI.3"))
plot(env_stack_selected)

# Save as multi-layer GeoTIFF
writeRaster(env_stack_selected, 
            filename = "env_stack_selected.tif", 
            overwrite = TRUE)


env_stack_selected <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/kano_ibadan/kano_ibadan_ento/env_stack_selected.tif")

predictors_a_subset <- rast("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/kano_ibadan/kano_ibadan_ento/predictors_a_subset.tif")













