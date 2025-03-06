### SPATA2 Spatial Trajectory Functions ###

library(SPATA2)
library(ggplot2)
library(stringr)
library(dplyr)
library(patchwork)

### Plots a box defining the full size of the trajectory based on a width provided
plotTrajectoryLayout <- function(data, start, end, width, ids = "horizontal_mid", color_by = "communities", clrp_adjust = NULL){
    
    traj_df <- tibble::tibble(
        x = c(start[1], end[1]),
        y = c(start[2], end[2])
    )

    starting_pos <- base::as.numeric(traj_df[1, c("x", "y")])
    final_pos <- base::as.numeric(traj_df[2, c("x", "y")])
    drvc <- final_pos - starting_pos  # Directional vector

    line_width <- width  # Use manually defined width

    # Compute midpoint of the main line
    midpoint <- (starting_pos + final_pos) / 2

    # Compute perpendicular vector
    perp_vec <- c(-drvc[2], drvc[1])  # Rotate by 90 degrees

    # Normalize and scale
    perp_vec <- (perp_vec / sqrt(sum(perp_vec^2))) * (line_width / 2)

    # Compute endpoints of the perpendicular line
    perp_start <- midpoint - perp_vec
    perp_end <- midpoint + perp_vec

    # Compute rectangle corners correctly
    corner1 <- starting_pos - perp_vec
    corner2 <- starting_pos + perp_vec
    corner3 <- final_pos + perp_vec
    corner4 <- final_pos - perp_vec

    # Create rectangle data frame
    box_df <- data.frame(
        x = c(corner1[1], corner2[1], corner3[1], corner4[1], corner1[1]),
        y = c(corner1[2], corner2[2], corner3[2], corner4[2], corner1[2])
    )
    
    # Generate plot
    p <- plotSpatialTrajectories(
        object = data, 
        ids = ids,
        color_by = color_by,
        clrp_adjust = clrp_adjust
    ) + 
        theme_classic() +
        # Add perpendicular line
        geom_segment(data = data.frame(
            x = perp_start[1], y = perp_start[2],
            xend = perp_end[1], yend = perp_end[2]
        ), 
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "red", size = 1.2, linetype = "dashed") +
        # Add midpoint
        geom_point(data = data.frame(
            x = midpoint[1], y = midpoint[2]
        ), aes(x = x, y = y), color = "blue", size = 3) +
        # Add rectangle (box)
        geom_polygon(data = box_df, aes(x = x, y = y), 
                     fill = NA, color = "green", size = 1.5)

    return(p)
}



## Plots the relative bins of the trajectory
PlotTrajectoryBins <- function(object, width, variables, id = "horizontal_mid", distance = "dte", unit = "px", resolution = SPATA2::recSgsRes(object),core = FALSE, angle_span = c(0,360), verbose = FALSE){
    
    coords_df <- getCoordsDf(object)
    
    distance <- SPATA2::getTrajectoryLength(object, id = id, unit = unit)
    resolution <- SPATA2::as_unit(resolution, unit = unit, object = object)
    coords_df_st <-
        SPATA2::getCoordsDfST(
          object = object,
          id = id,
          width = width,
          variables = variables,
          dist_unit = "px", # ensure that distance is computed in correct unit
          verbose = verbose
        ) %>%
        dplyr::filter(rel_loc == "inside")
    p <- ggplot(coords_df_st, aes(x = x, y = y, color = bins_dist)) +
      geom_point(size = 2, alpha = 0.7) +  # Scatter points
      scale_color_viridis_d() +  # Better color scale for bins
      labs(title = "Scatter Plot of Bins", x = "X Coordinate", y = "Y Coordinate", color = "Bins") +
      theme_classic() + theme(legend.position = "none")  # Removes the legend
    return(p)

}