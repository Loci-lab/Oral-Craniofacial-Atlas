# Load necessary libraries
library(sf)
library(dplyr)
library(tidyr)

# Define helper functions
create_windows <- function(data, window_size, step_size) {
  bbox <- st_bbox(data)
  x_breaks <- seq(bbox["xmin"], bbox["xmax"] - window_size, by = step_size)
  y_breaks <- seq(bbox["ymin"], bbox["ymax"] - window_size, by = step_size)
  expand.grid(x = x_breaks, y = y_breaks)
}

unique_neighbors <- function(results) {
  results %>%
    distinct() %>%
    rowwise() %>%
    mutate(neighbors = list(unique(neighbors)))
}

find_neighbors_in_window <- function(window, data, window_size) {
  xmin <- window$x
  xmax <- window$x + window_size
  ymin <- window$y
  ymax <- window$y + window_size
  
  window_data <- data %>%
    filter(st_coordinates(.)[, 1] >= xmin & st_coordinates(.)[, 1] < xmax &
             st_coordinates(.)[, 2] >= ymin & st_coordinates(.)[, 2] < ymax)
  
  if (nrow(window_data) < 2) {
    return(data.frame(cell_id = character(0), neighbors = I(list())))
  }
  
  distances <- st_distance(window_data)
  diag(distances) <- NA  # Exclude self-distance
  
  neighbors <- lapply(seq_len(nrow(window_data)), function(i) {
    near <- order(distances[i, ], na.last = NA)[1:5]  # Get indices of the 5 closest neighbors
    if (length(near) > 0) as.character(window_data$cell_id[near]) else NA
  })
  
  data.frame(cell_id = window_data$cell_id, neighbors = I(neighbors))
}

# Parameters
window_size <- 300  # Size of the window
step_size <- 250  # Step size for sliding window

# Initialize an empty list to store group results
group_results <- list()
data2$CellID=1:nrow(data2)
# Iterate over each group
unique_groups <- unique(data2$Group)
for (group_name in unique_groups) {
  cat("Processing Group:", group_name, "\n")
  
  # Subset data for the current group
  group_data <- data2 %>% filter(Group == group_name)
  
  # Convert to sf object
  data_sub_sf <- st_as_sf(group_data, coords = c("X", "Y"), crs = NA)
  data_sub_sf$cell_id <- group_data$CellID
  
  # Create sliding windows
  windows <- create_windows(data_sub_sf, window_size, step_size)
  
  # Find neighbors within each window
  results <- do.call(rbind, lapply(seq_len(nrow(windows)), function(i) {
    find_neighbors_in_window(windows[i, ], data_sub_sf, window_size)
  }))
  
  # Remove duplicates
  results <- unique_neighbors(results)
  
  # Combine results by group
  results_df <- results %>%
    group_by(cell_id) %>%
    summarise(neighbors = list(unique(unlist(neighbors)))) %>%
    filter(!is.na(neighbors)) %>%
    filter(!is.na(cell_id))
  
  # Expand results for pairwise connections
  results_df_expanded <- results_df %>%
    unnest(neighbors)
  
  colnames(results_df_expanded)[2] <- "Neighbor_ID"
  
  # Store in the group results list
  group_results[[group_name]] <- results_df_expanded
  
  cat("Completed Group:", group_name, "\n")
}

# Combine all group results into a single data frame
final_results <- bind_rows(group_results, .id = "Group")

final_results=final_results[which(final_results$Neighbor_ID!="NA"),]

final_results$cell_id=as.numeric(final_results$cell_id)
final_results$Neighbor_ID=as.numeric(final_results$Neighbor_ID)

final_results$TACIT_cellID=data2$TACIT[final_results$cell_id]
final_results$TACIT_Neighbor=data2$TACIT[final_results$Neighbor_ID]

final_results$X_cellID=data2$X[final_results$cell_id]
final_results$X_Neighbor=data2$X[final_results$Neighbor_ID]

final_results$Y_cellID=data2$Y[final_results$cell_id]
final_results$Y_Neighbor=data2$Y[final_results$Neighbor_ID]



# Prepare data for heatmap
# Step 1: Create the connection matrix
connection_matrix <- table(final_results$TACIT_cellID, final_results$TACIT_Neighbor)

# Step 2: Convert the table to a matrix for manipulation
connection_matrix <- as.matrix(connection_matrix)

# Step 3: Add the matrix to its transpose to make it symmetric
symmetric_matrix <- connection_matrix + t(connection_matrix)

# Optional: Convert back to a table if needed
symmetric_table <- as.table(symmetric_matrix)

library(pheatmap)
pheatmap(
  connection_matrix,
  color = colorRampPalette(c("white", "red"))(50),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  main = "Connection Frequency Heatmap"
)










library(ggplot2)
final_results_sub=final_results[which(final_results$Group=="Group 1"),]
# Plot connections

ggplot(final_results_sub) +
  geom_segment(aes(x = X_cellID, y = Y_cellID, xend = X_Neighbor, yend = Y_Neighbor, color = TACIT_cellID), size = 0.5)+ 
  theme_classic(base_size = 15) +
  labs(x = "X", y = "Y", title = "Neighbor Connections", color = "Cell Type") +
  guides(color = guide_legend(override.aes = list(size = 3)))

