######### Color Palette ########################################################
color_palette <- colorRampPalette(c("darkblue", 
                                    "lightblue",
                                    "white",
                                    "orange",
                                    "red",
                                    "darkred"))

library(RColorBrewer)
jet.colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

# Generate the desired number of colors from this palette
nbcol <- 50
color <- jet.colors(nbcol)
ramp <- colorRamp(c("blue", "white", "red"))
################################################################################