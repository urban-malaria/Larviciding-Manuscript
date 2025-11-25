# Load necessary libraries
library(dismo)  
library(raster) 
library(maps)   
library(rJava)
library(terra)
library(sf)
library(dismo)
library(randomForest)
library(gbm)
library(caret)

##Run Java locally before loading library
Sys.setenv(JAVA_HOME = "C:/Program Files/Eclipse Adoptium/jdk-21.0.7.6-hotspot")
library(rJava)

#Load occurence data
lav_data <- read.csv("C:/Users/ebamgboye/Urban Malaria Proj Dropbox/urban_malaria/data/nigeria/kano_ibadan/kano_ibadan_ento/Wet Season Data_Ibadan/lav_coords_bsw1.csv")

anopheles_sites <- lav_data %>% 
  dplyr::filter(anophw == "Yes")

# Keep only latitude and longitude columns
occurrences_in_agugu <- anopheles_sites[, 2:3]
colnames(occurrences_in_agugu) <- c("lon", "lat")


# Ensure they're numeric
occurrences_in_agugu$lon <- as.numeric(occurrences_in_agugu$lon)
occurrences_in_agugu$lat <- as.numeric(occurrences_in_agugu$lat)

##Ensure points fall within agugu ward
# Step 1: Convert to sf object 
occurrences_sf <- st_as_sf(occurrences_in_agugu, coords = c("lon", "lat"), crs = 4326)

# Step 2: Reproject to match the CRS of df_ib_a
occurrences_proj <- st_transform(occurrences_sf, crs = st_crs(df_ib_a))

# Plot to verify visually
plot(st_geometry(df_ib_a), col = "lightblue", main = "Occurrence Points within Agugu")
plot(st_geometry(occurrences_sf), col = "red", pch = 20, add = TRUE)


# Plot using ggplot with geom_sf for both layers
ggplot() +
  geom_sf(data = df_ib_a, fill = NA, color = "black") +
  geom_sf(data = occurrences_proj, color = "red", size = 2, alpha = 0.8) +
  labs(
    title = "Positive Breeding Sites in Agugu",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  coord_sf()
  