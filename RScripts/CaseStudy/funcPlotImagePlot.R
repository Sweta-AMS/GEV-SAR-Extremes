# # Create image plot
# imagePlot <- function(loc, y)
# {
#   ug <- sort(unique(loc[,1]))
#   vg <- sort(unique(loc[,2]))
#   
#   ind_ug <- match(loc[,1], ug)
#   ind_vg <- match(loc[,2], vg)
#   
#   img <- matrix(NA,
#                 nrow = length(ug),
#                 ncol = length(vg))
#   data_index <- cbind(ind_ug, ind_vg)
#   
#   img[data_index] <- y
#   return(list(ug, vg, img))
# }

imagePlot <- function(loc, y) {
  # Get unique x and y locations (grid points)
  ug <- sort(unique(loc[,1]))  # Unique x (longitude)
  vg <- sort(unique(loc[,2]))  # Unique y (latitude)
  
  # Find the index positions in the grid
  ind_ug <- match(loc[,1], ug)  # Row indices
  ind_vg <- match(loc[,2], vg)  # Column indices
  
  # Create an empty matrix with NA values
  img <- matrix(NA, 
                nrow = length(vg), 
                ncol = length(ug))  # Reverse dimensions
  
  # ✅ Correct indexing method
  for (i in 1:length(y)) {
    img[ind_vg[i], ind_ug[i]] <- y[i]  # Use row and column indexing separately
  }
  
  return(list(ug = ug, vg = vg, img = t(img)))  # Transpose for correct plotting
}

# # Example Usage
# set.seed(123)
# loc <- matrix(runif(100, min = 0, max = 10), ncol = 2)  # Random locations
# y <- rnorm(50)  # Example values
# 
# # Generate the structured image data
# image_data <- imagePlot(loc, y)
# 
# # Plot the image
# image.plot(image_data$ug, image_data$vg, image_data$img, col = terrain.colors(100), 
#       xlab = "Longitude", ylab = "Latitude", main = "Fixed Spatial Image Plot")
