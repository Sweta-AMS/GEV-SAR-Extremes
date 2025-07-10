rm(list=ls())
source("~/Desktop/GEV-SAR/RScripts/CaseStudy/RCM-dataLoad.R")

## packages required
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(fields)
library(maps)
# library(maptools)
library(rworldmap) 

# Spatial Domain Plot:
set.panel(1,1)

### North-America-with water
image.plot(x=lon_data,
           y=lat_data,
           matrix(1,
                  nrow=nrow(lon_data),
                  ncol=ncol(lon_data)),
           col=alpha("grey", 0.5),
           xlab="Longitude",
           ylab="Latitude",
           main="Spatial Domain NA-CORDEX")
world(add=TRUE)

### ---  Just extract the USA spatial fields --- ###
# Convert longitude & latitude to spatial points
coords_df <- data.frame(lon = as.vector(lon_data),
                        lat = as.vector(lat_data))
coords_sf <- st_as_sf(coords_df,
                      coords = c("lon", "lat"),
                      crs = 4326)  # Convert to spatial object

# Load USA boundaries from Natural Earth
usa <- ne_countries(country = "United States of America",
                    returnclass = "sf")

# Perform spatial filtering: Keep only points inside the USA
usa_points <- coords_sf[st_within(coords_sf, 
                                  usa,
                                  sparse = FALSE), ]

# Extract subsetted longitude & latitude
usa_lon <- st_coordinates(usa_points)[, 1]
usa_lon <- usa_lon[(usa_lon > -124)]
usa_lat <- st_coordinates(usa_points)[, 2]
usa_lat <- usa_lat[(usa_lat < 52)]

# Filter precipitation data for the USA
library(FNN)  # For nearest neighbor search

find_nearest_index <- function(values, reference) {
  nn <- get.knnx(reference, values, k = 1)  # Find nearest neighbors
  return(nn$nn.index[, 1])  # Return the index of the closest match
}

# Find indices for USA coordinates
lon_indices <- find_nearest_index(usa_lon, as.vector(lon_data))
lat_indices <- find_nearest_index(usa_lat, as.vector(lat_data))

# Extract precipitation data for the USA
## Convert 3D to 2D precipitation data:
precip_data2D <- matrix(precip_data, 
                        nrow=(dim(precip_data)[1]*dim(precip_data)[2]), 
                        ncol=dim(precip_data)[3])
dim(precip_data2D)


valid_indices <- intersect(lon_indices, lat_indices)

# Extract precipitation data for the USA only for valid indices
usa_precip_data <- precip_data2D[valid_indices, , drop = FALSE]
dim(usa_precip_data ) # 7324 11323

lon_vec <- as.vector(lon_data)
lat_vec <- as.vector(lat_data)

lon_usa <- lon_vec[valid_indices]
lat_usa <- lat_vec[valid_indices]

plot(lon_usa, lat_usa, cex=0.5)
world(add=T, col='grey', lwd=3)
