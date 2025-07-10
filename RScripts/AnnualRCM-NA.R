rm(list=ls())
library(tictoc)
library(fields)
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/CaseStudy/RCM-dataLoad.R")


# Spatial Domain Plot:



## -- Annual-Maximum -- :
# Extract the year from the time variable
library(lubridate)
years <- year(time_dates)

# Initialize an array to store annual max precipitation (lat, lon, years)
annual_maxima <- array(NA, 
                       dim=c(dim(precip_data)[1],
                             dim(precip_data)[2],
                             length(unique(years))))

# Loop through each year and compute annual max at each grid point
for (i in seq_along(unique(years))) {
  print(i)
  year_mask <- (years == unique(years)[i])  # Mask for the current year
  annual_maxima[,,i] <- apply(precip_data[, ,year_mask],
                              c(1,2),
                              max,
                              na.rm = TRUE)
}
# print dimensions of the computed annual maxima
print(dim(annual_maxima))  # Should be (lat, lon, years)

# # Save the data
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/CaseStudy")
# save(annual_maxima, file='AnnualMaxima-acrossNA.npy')



# ### North-America- with water 
# image.plot(x=lon_data,
#            y=lat_data,
#            matrix(1,
#                   nrow=nrow(lon_data),
#                   ncol=ncol(lon_data)),
#            col=alpha("grey", 0.5),
#            xlab="Longitude",
#            ylab="Latitude",
#            main="Spatial Domain NA-CORDEX")
# world(add=TRUE)


## ---  Focusing just on the land area  ---:
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(fields)
library(maps)
# library(maptools)
library(rworldmap)   # getMap

land <- getMap(resolution = "low")
plot(land,
     col = "grey80",
     border = "black",
     xlim = c(-160, -50), 
     ylim = c(10, 80),
     xlab = "Longitude", 
     ylab = "Latitude",
     main = "Spatial Domain NA-CORDEX (Land Only)")

### --- America ---:
land <- getMap(resolution = "low")

# # Extract only the USA
# usa <- land[land$ADMIN == "United States of America", ]
# 
# # Plot only the USA region
# plot(usa, col = "grey80", border = "black",
#      xlim = c(-130, -68), ylim = c(20, 47),  # Focus on CONUS (Continental US)
#      xlab = "Longitude", ylab = "Latitude",
#      main = "Spatial Domain: USA Only")
# 
# 
# 
# ###  --- Just extracting locations for the USA ---:
# coords_df <- data.frame(lon = as.vector(lon_data), 
#                         lat = as.vector(lat_data))
# coords_sf <- st_as_sf(coords_df,
#                       coords = c("lon", "lat"),
#                       crs = 4326)  # Convert to spatial format
# 
# # Perform spatial intersection: Keep only points within the USA
# usa_points <- coords_sf[st_within(coords_sf,
#                                   usa,
#                                   sparse = FALSE), ]
# 
# # Extract subsetted coordinates
# usa_lon <- st_coordinates(usa_points)[, 1]
# usa_lat <- st_coordinates(usa_points)[, 2]
# 
# # Load USA boundaries
# usa <- ne_countries(country = "United States of America",
#                     returnclass = "sf")
# 
# # Use nearest-neighbor search instead of exact matching
# usa_lon_idx <- sapply(usa_lon, 
#                       function(x) which.min(abs(lon_data - x)))
# usa_lat_idx <- sapply(usa_lat, function(x) which.min(abs(lat_data - x)))
# 
# # Ensure indices are valid
# valid_idx <- !is.na(usa_lon_idx) & !is.na(usa_lat_idx)
# 
# # Extract USA precipitation data
# usa_precip_data <- precip_data[usa_lon_idx[valid_idx], usa_lat_idx[valid_idx], ]
# 
# # Print dimensions to confirm correct extraction
# print(dim(usa_precip_data))
# 
# 
# 
# # 
# # # Filter precipitation data for the USA
# # usa_precip_data <- precip_data[match(usa_lon, as.vector(lon_data)),
# #                                match(usa_lat, as.vector(lat_data)),]
# 
# # Plot USA subset
# plot(st_geometry(usa), col = "grey80", border = "black",
#      xlab = "Longitude", ylab = "Latitude", main = "NA-CORDEX Subset Over USA")
# points(usa_lon, usa_lat, col = "blue", pch = 16, cex = 0.5)
# 
# # Save the subset as a new NetCDF file (Optional)
# usa_nc <- rast(nrows = length(usa_lat),
#                ncols = length(usa_lon),
#                crs = "EPSG:4326",
#                xmin = min(usa_lon),
#                xmax = max(usa_lon),
#                ymin = min(usa_lat),
#                ymax = max(usa_lat))
# 
# values(usa_nc) <- usa_precip_data  # Assign precipitation data
# 
# writeRaster(usa_nc,
#             filename = "USA_NA-CORDEX_subset.nc",
#             format = "CDF", 
#             overwrite = TRUE)
# 
# print("USA subset of NA-CORDEX saved successfully.")


