## Loading the data
rm(list=ls())
library(tictoc)

## ERA-Int- WRF data timeline: 1980-01-01to 2010-12-31
## ~79 km resolution
# Load necessary libraries
library(ncdf4)    
library(fields)  
library(ggplot2)
library(sf)
library(raster)
library(dplyr)
library(terra)
library(scales)  # For transparency


# Set file path
file_path <- "/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Data/prec.eval.ERA-Int.WRF.day.NAM-22.raw.nc"

# Open the NetCDF file
nc_file <- nc_open(file_path)

# List all available variables in the NetCDF file
variables <- names(nc_file$var)
print("Available Variables:")
print(variables)

# Extract Precipitation, Latitude, Longitude, and Time
precip_data <- ncvar_get(nc_file,
                         "prec")  # Check if variable name is "pr" instead of "prec"
lat_data <- ncvar_get(nc_file, 
                      "lat")   

lon_data <- ncvar_get(nc_file,
                      "lon")   
time_data <- ncvar_get(nc_file,
                       "time")  # Extract actual time dimension
time_bnds <- ncvar_get(nc_file,
                       "time_bnds")  # Check time boundaries

proj_info <- ncatt_get(nc_file, "Lambert_Conformal")  # Extract projection details
print(proj_info)  # Display Lambert projection parameters

# Print dimensions to verify structure
print(paste("Precipitation Dimensions:",
            paste(dim(precip_data),
                  collapse=" x ")))
print(paste("Latitude Dimensions:",
            paste(dim(lat_data),
                  collapse=" x ")))
print(paste("Longitude Dimensions:",
            paste(dim(lon_data),
                  collapse=" x ")))
print(paste("Time Dimensions:",
            length(time_data)))



## -- Extract time variable --
time_data <- ncvar_get(nc_file,
                       "time")

# Get time units from metadata
time_units <- ncatt_get(nc_file,
                        "time", 
                        "units")$value
print(time_units)  # Example: "days since 1979-01-01 00:00:00"

# Extract the reference date from the units
time_origin <- gsub("days since ",
                    "", 
                    time_units)  # Extract only the date part

# Convert numeric time values to actual dates
time_dates <- as.Date(time_data,
                      origin=time_origin)

# Get start and end time
start_time <- min(time_dates,
                  na.rm = TRUE)
end_time <- max(time_dates,
                na.rm = TRUE)

# Print results
print(paste("Start Time:", start_time))
print(paste("End Time:", end_time))


# ## Data visvualization
# image.plot(lon_data,
#            lat_data,
#            (precip_data[,,1]-median((precip_data[,,1])))/IQR(precip_data[,,1]))

