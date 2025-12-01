# Load required libraries
library(dplyr)
library(geosphere)
library(sf)
library(purrr)

##Load breeding site data(from data preparation scripts)
lav_bs_dry <- lav_df-dry

##Extract for Slums (Agugu) Dry Season
lav_bs_slum_dry <- lav_bs_dry %>% 
  dplyr::filter(`Settlement Type` == "Slum") %>% 
  filter(!`site_label` %in% c(17, 27, 40))

breeding_sites <- lav_bs_slum_dry

# Function to generate sampled paths
set.seed(123)  # For reproducibility
generate_sample_paths <- function(data, n_paths, fixed_start, fixed_end, sample_fraction) {
  map_dfr(1:n_paths, ~ {
    # Sample points excluding the fixed start and end
    sampled_points <- data %>%
      filter(!site_label %in% c(fixed_start$site_label, fixed_end$site_label)) %>%
      sample_n(size = floor(sample_fraction * (nrow(data) - 2)), replace = FALSE)
    
    # Ensure the pathway starts at fixed_start and ends at fixed_end
    path <- bind_rows(fixed_start, sampled_points, fixed_end) %>%
      mutate(Path_ID = .x)
    
    return(path)
  })
}

# Function to compute total distance for each path
compute_path_distance <- function(df) {
  df %>%
    arrange(Path_ID, latitude, longitude) %>%
    group_by(Path_ID) %>%
    summarise(
      Total_Distance_km = sum(geosphere::distVincentySphere(
        cbind(longitude, latitude), 
        cbind(dplyr::lead(longitude), dplyr::lead(latitude))
      ), na.rm = TRUE) / 1000  # Convert meters to km
    )
}

# Function to compute Sample Area (km²)
compute_sample_area <- function(total_distance, effective_strip_width) {
  total_distance * (effective_strip_width)
}

# Function to compute Breeding Site Density
compute_breeding_site_density <- function(num_positive_sites, sampled_area) {
  ifelse(sampled_area > 0, num_positive_sites / sampled_area, NA)  # Avoid division by zero
}

##Using Field data
# Define start and end points
fixed_start <- breeding_sites %>% filter(site_label == 1)
fixed_end <- breeding_sites %>% filter(site_label == 50)

sampled_paths <- generate_sample_paths(
  breeding_sites, 
  n_paths = 100, 
  fixed_start = fixed_start, 
  fixed_end = fixed_end, 
  sample_fraction = 0.8
)

##Convert to sf points for visualization
sampled_paths_sf <- st_as_sf(sampled_paths, coords = c("longitude", "latitude"), crs = 4326)
breeding_sites_sf <- st_as_sf(breeding_sites, coords = c("longitude", "latitude"), crs = 4326)


##Visualize paths
# Ensure sampled_paths is an sf object
if (!inherits(sampled_paths, "sf")) {
  sampled_paths_sf <- st_as_sf(sampled_paths, coords = c("longitude", "latitude"), crs = 4326)
} else {
  sampled_paths_sf <- sampled_paths
}

# Convert points into LINESTRING for each Path_ID
sampled_paths_sf <- sampled_paths_sf %>%
  group_by(Path_ID) %>%
  summarise(geometry = st_combine(geometry)) %>%  # Merge points into a multipoint geometry
  st_cast("LINESTRING")  # Convert multipoint geometry into a LINESTRING

# Plot the paths
path_plot <- ggplot() +
  # Add ward boundaries
  geom_sf(data = df_ib_a, fill = NA, color = "gray60", size = 0.5) +
  
  # Add sampled paths
  geom_sf(data = sampled_paths_sf, aes(color = as.factor(Path_ID)), size = 0.8, alpha = 0.5, show.legend = FALSE) +
  
  # Add breeding site points
  geom_sf(data = breeding_sites_sf, color = "black", size = 2, shape = 21, fill = "yellow") +
  
  # Labels and theme
  labs(title = "Sampled Pathways Between Fixed Breeding Sites",
       subtitle = "Each line represents a different sampled path",
       x = "Longitude", y = "Latitude") +
  theme_manuscript()

ggsave(paste0(LuDir, '/plots/', Sys.Date(), 'Sampled Pathways Between Fixed Breeding Sites(slum dryseason).pdf'), path_plot, width = 8, height = 6)

# Compute path distances
path_distances <- compute_path_distance(sampled_paths)

# Assuming effective strip width is 0.05 km (50m)
effective_strip_width <- 0.05  # km

# Compute sample area for each path
path_distances <- path_distances %>%
  mutate(Sample_Area_km2 = compute_sample_area(Total_Distance_km, effective_strip_width))

# Count the number of positive breeding sites per path
positive_breeding_sites <- sampled_paths %>%
  filter(anophw == "Yes") %>%  # Filter positive breeding sites
  group_by(Path_ID) %>%
  summarise(Number_Positive_Sites = n(), .groups = "drop")

# Merge positive breeding sites count with path distances
path_summary <- path_distances %>%
  left_join(positive_breeding_sites, by = "Path_ID") %>%
  mutate(Number_Positive_Sites = ifelse(is.na(Number_Positive_Sites), 0, Number_Positive_Sites))  # Replace NA with 0

# Apply function to compute breeding site density
path_summary <- path_summary %>%
  mutate(Breeding_Site_Density = compute_breeding_site_density(Number_Positive_Sites, Sample_Area_km2))

# View the results
head(path_summary)

# Calculate summary statistics
mean_density <- mean(path_summary$Breeding_Site_Density)
sd_density <- sd(path_summary$Breeding_Site_Density)
min_density <- min(path_summary$Breeding_Site_Density)
max_density <- max(path_summary$Breeding_Site_Density)

# Determine annotation position
x_pos <- quantile(path_summary$Breeding_Site_Density, 0.75) # 75th percentile for better placement
y_pos <- max(table(cut(path_summary$Breeding_Site_Density, breaks = 30))) * 0.8 # Just below the highest bar

ggplot(path_summary, aes(x = Breeding_Site_Density)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Breeding Site Densities with ESW of 50m(slumdry)",
       x = "Density (sites/km²)",
       y = "Frequency") +
  theme_minimal()+
  annotate("text", x = x_pos, y = y_pos, 
           label = paste("Mean: ", round(mean_density, 2), "\n",
                         "SD: ", round(sd_density, 2), "\n",
                         "Min: ", round(min_density, 2), "\n",
                         "Max: ", round(max_density, 2)),
           vjust = -0.9, hjust = -0.9, size = 3, color = "black")


##Convert summaries into a dataframe
Agudry_summary <- data.frame(
  Ward         = "Agugu",
  mean_density = mean(path_summary$Breeding_Site_Density, na.rm = TRUE),
  sd_density   = sd(path_summary$Breeding_Site_Density, na.rm = TRUE),
  min_density  = min(path_summary$Breeding_Site_Density, na.rm = TRUE),
  max_density  = max(path_summary$Breeding_Site_Density, na.rm = TRUE)
)


Agudry_summary$season <- "Dry"

Agudry_summary$settlment <- "Slum"




