# Load required libraries
library(readr)
library(readxl)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)  # For kde2d function
library(igraph)
library(RColorBrewer)

# Set working directory
setwd("C:/Users/huynhk4/Documents/")


ligands_list <- readxl::read_xlsx("LR methods/human_lr_pair_full.xlsx")


data2 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_Tongue_ID62_TACIT.csv") #Import datasets
data2$Group="Tongue62"
data2$Group_MG="Mucosal"

data3 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_Tongue12_TACIT.csv") #Import datasets
data3$Group="Tongue12"
data3$Group_MG="Mucosal"

data4 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_Gingiva_2_TACIT.csv") #Import datasets
data4$Group="Gingiva2"
data4$Group_MG="Mucosal"

data5 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_Gingiva1_TACIT.csv") #Import datasets
data5$Group="Gingiva1"
data5$Group_MG="Mucosal"

data6 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_BuccalMucosa_ID_5_TACIT.csv") #Import datasets
data6$Group="BuccalMucosa"
data6$Group_MG="Mucosal"

data7 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_geneBuccalMucosaMSG1_TACIT.csv") #Import datasets
data7$Group="MSG1"
data7$Group_MG="Glands"

data8 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_geneBuccalMucosaMSG8_TACIT.csv") #Import datasets
data8$Group="MSG8"
data8$Group_MG="Glands"

data9 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_geneBuccalMucosaMSG9_TACIT.csv") #Import datasets
data9$Group="MSG9"
data9$Group_MG="Glands"

data10 <- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_Parotid_39_TACIT.csv") #Import datasets
data10$Group="Parotid39"
data10$Group_MG="Glands"

data11<- read_csv("TACIT_annotate_MERSCOPE/cell_by_gene_Parotid_ID_35_TACIT.csv") #Import datasets
data11$Group="Parotid35"
data11$Group_MG="Glands"

data12<- read_csv("TACIT_annotate_MERSCOPE/cell_by_geneParotid_20_TACIT.csv") #Import datasets
data12$Group="Parotid20"
data12$Group_MG="Glands"

data13<- read_csv("TACIT_annotate_MERSCOPE/cell_by_geneParotid_40_TACIT.csv") #Import datasets
data13$Group="Parotid40"
data13$Group_MG="Glands"

data14<- read_csv("TACIT_annotate_MERSCOPE/cell_by_geneSubmandibular_53_56.csv") #Import datasets
data14$Group="Submandibular53_56"
data14$Group_MG="Glands"



standardize_and_offset <- function(data, x_offset = 0, y_offset = 0) {
  # Standardize column names
  common_columns <- intersect(colnames(data9), intersect(colnames(data2), colnames(data14)))
  data <- data[common_columns]  # Reorder columns
  
  # Apply offset
  data <- data %>%
    mutate(X = X + x_offset, Y = Y + y_offset)
  
  return(data)
}

# List of datasets
datasets <- list(data10,data11,data12,data13)

# Initialize offsets
x_offset <- 0
y_offset <- 0
increment <- 10000  # Smaller increment to make the datasets closer

# Standardize and combine datasets with offsets
data1 <- datasets[[1]]
for (i in 2:length(datasets)) {
  x_offset <- x_offset + increment
  y_offset <- y_offset  + increment
  datasets[[i]] <- standardize_and_offset(datasets[[i]], x_offset, y_offset)
  data1 <- bind_rows(data1, datasets[[i]])
}

# View the combined data
head(data1)


# Find common columns between ligands/receptors and data1
common <- intersect(c(ligands_list$Ligante, ligands_list$Receptores), colnames(data1))





#---------------------------------------------

# Do not change

data_sub=data1[,c("X","Y","TACIT","Group",common)]
data_sub=data.frame(cellID=1:nrow(data_sub),data_sub)
selected_rows <- rowSums(data1[, common] > 0) > 2
data_sub=data_sub[selected_rows,]



library(sf)
library(dplyr)

# Assuming you have data in data_sub with columns X, Y, cellID, and Group

# Convert data to an sf object without a CRS (non-geographic)
data_sub_sf <- st_as_sf(data_sub, coords = c("X", "Y"), crs = NA)
data_sub_sf$cell_id <- data_sub$cellID

# Function to create a grid of windows
create_windows <- function(data, window_size, step_size) {
  bbox <- st_bbox(data)
  x_breaks <- seq(bbox["xmin"], bbox["xmax"] - window_size, by = step_size)
  y_breaks <- seq(bbox["ymin"], bbox["ymax"] - window_size, by = step_size)
  return(expand.grid(x = x_breaks, y = y_breaks))
}


unique_neighbors <- function(results) {
  results %>%
    distinct() %>%
    rowwise() %>%
    mutate(neighbors = list(unique(neighbors)))
}

find_neighbors_in_window <- function(window, data, window_size, distance_cutoff) {
  xmin <- window$x
  xmax <- window$x + window_size
  ymin <- window$y
  ymax <- window$y + window_size
  
  window_data <- data %>%
    filter(st_coordinates(.)[,1] >= xmin & st_coordinates(.)[,1] < xmax &
             st_coordinates(.)[,2] >= ymin & st_coordinates(.)[,2] < ymax)
  
  if (nrow(window_data) < 2) {
    return(data.frame(cell_id = character(0), neighbors = I(list())))
  }
  
  distances <- st_distance(window_data)
  diag(distances) <- NA
  distances[distances > distance_cutoff] <- NA
  
  neighbors <- lapply(seq_len(nrow(window_data)), function(i) {
    near <- which(!is.na(distances[i, ]))
    if (length(near) > 0) as.character(window_data$cell_id[near]) else NA
  })
  
  return(data.frame(cell_id = window_data$cell_id, neighbors = I(neighbors)))
}
# Parameters
window_size <- 500  # Size of the window
step_size <- 400  # Step size for sliding window
distance_cutoff <- 50  # Distance cutoff for neighbors

# Create windows with overlap
windows <- create_windows(data_sub_sf, window_size, step_size)

# Find neighbors within each window and combine results
results <- do.call(rbind, lapply(seq_len(nrow(windows)), function(i) {
  find_neighbors_in_window(windows[i, ], data_sub_sf, window_size, distance_cutoff)
}))

# Remove duplicates
results <- unique_neighbors(results)


# Combine results by group
results_df <- results %>%
  group_by(cell_id) %>%
  summarise(neighbors = list(unique(unlist(neighbors))))

# Display results
print(results_df)



# Assuming 'results_df' is already loaded and contains the neighbors column as lists
results_df_expanded <- results_df %>%
  mutate(neighbors = ifelse(is.na(neighbors), list(NA), neighbors)) %>%  # Ensure NA is handled as a list for consistency
  unnest(neighbors) %>%
  rename(Neighbor_ID = neighbors)  # Rename the column for clarity

results_df_expanded1=results_df_expanded


data2_sub=data1[,c("X","Y","TACIT","Group",common)]
results_df_sum=results_df_expanded1[which(results_df_expanded1$Neighbor_ID!="NA"),]

data2_sub=data.frame(cell_id=1:nrow(data2_sub),data2_sub)

# Function to check interactions using vectorized operations and pre-indexing
check_interaction <- function(data, ligands, results) {
  # Initialize the interaction matrix
  interaction_matrix <- results
  interaction_names <- paste(ligands$Ligante, ligands$Receptores, sep = "_")
  interaction_matrix[interaction_names] <- 0
  
  # Create a data frame with only relevant columns for easier manipulation
  data_subset <- data[, c("cell_id", ligands$Ligante, ligands$Receptores)]
  
  # Pre-calculate which cells express each ligand and receptor
  ligand_expressed <- sapply(ligands$Ligante, function(l) data_subset[data_subset[[l]] > 0, "cell_id"])
  receptor_expressed <- sapply(ligands$Receptores, function(r) data_subset[data_subset[[r]] > 0, "cell_id"])
  
  # Vectorized comparison for each ligand-receptor pair
  for (i in seq_along(interaction_names)) {
    if(length(seq_along(interaction_names))==1){
      ligand_cells <- ligand_expressed
      receptor_cells <- receptor_expressed
    }else{
      ligand_cells <- ligand_expressed[[i]]
      receptor_cells <- receptor_expressed[[i]]
    }
    # Find matching cell_id and Neighbor_ID pairs
    interactions <- results$cell_id %in% ligand_cells & results$Neighbor_ID %in% receptor_cells
    interaction_matrix[interactions, interaction_names[i]] <- 1
  }
  
  return(interaction_matrix)
}



ligands_list=ligands_list[which(ligands_list$Ligante%in%common),]
ligands_list=ligands_list[which(ligands_list$Receptores%in%common),]


# Apply the function
interaction_matrix <- check_interaction(data2_sub, ligands_list, results_df_sum)

table(interaction_matrix$CXCL12_CXCR4)

unique_columns <- !duplicated(colnames(interaction_matrix))
interaction_matrix <- interaction_matrix[, unique_columns]


interaction_matrix <- interaction_matrix %>%
  mutate(across(where(is.list), ~ map_dbl(.x, ~ as.numeric(.x[[1]])), .names = "first_{col}"))


# Add coordinates for cell_id
interaction_matrix <- interaction_matrix %>%
  left_join(data2_sub %>% dplyr::select(cell_id, X, Y), by = "cell_id") %>%
  dplyr::rename(X_cell = X, Y_cell = Y)

interaction_matrix$Neighbor_ID=as.integer(interaction_matrix$Neighbor_ID)

# Add coordinates for Neighbor_ID
interaction_matrix <- interaction_matrix %>%
  left_join(data2_sub %>% dplyr::select(cell_id, X, Y), by = c("Neighbor_ID" = "cell_id")) %>%
  dplyr::rename(X_neighbor = X, Y_neighbor = Y)

TACIT_Ligands=data1$TACIT[interaction_matrix$cell_id]
TACIT_Receptors=data1$TACIT[interaction_matrix$Neighbor_ID]


count=interaction_matrix[,3:165]
interaction_matrix=interaction_matrix[which(rowSums(count)>0),]



count=interaction_matrix[,3:165]

data1_sub_plot <- data.frame(
  X = interaction_matrix$X_cell,
  Y = interaction_matrix$Y_cell,
  ID = 1:nrow(interaction_matrix),
  count
)
library(dplyr)
library(tidyr)

# Assume data1_sub_plot is already defined and includes interaction counts
# Define grid size and calculate dimensions of each grid cell
grid_size <- 10000
x_breaks <- seq(min(data1_sub_plot$X), max(data1_sub_plot$X), length.out = grid_size + 1)
y_breaks <- seq(min(data1_sub_plot$Y), max(data1_sub_plot$Y), length.out = grid_size + 1)
cell_width <- x_breaks[2] - x_breaks[1]
cell_height <- y_breaks[2] - y_breaks[1]

data1_sub_plot <- data1_sub_plot %>%
  mutate(
    grid_x = cut(X, breaks = x_breaks, include.lowest = TRUE, labels = FALSE),
    grid_y = cut(Y, breaks = y_breaks, include.lowest = TRUE, labels = FALSE),
    grid_id = interaction(grid_x, grid_y, drop = TRUE)  # drop=TRUE to eliminate unused levels
  )



data_wide <- data1_sub_plot


# Summarize the data by grid_id
data_wide <- data1_sub_plot %>%
  group_by(grid_id) %>%
  summarize(
    avg_grid_x = mean(grid_x, na.rm = TRUE),
    avg_grid_y = mean(grid_y, na.rm = TRUE),
    across(all_of(colnames(count)), sum, na.rm = TRUE)
  )

# View the data structure
head(data_wide)

library(MASS)
library(dplyr)
library(tidyr)
library(fields)

# Function to compute density for each ligand-receptor pair with adjustments
compute_density_exact <- function(data, col, x_col, y_col) {
  data_sub <- data %>%
    filter(!!sym(col) == 1)
  
  if (nrow(data_sub) > 1) {
    density <- kde2d(data_sub[[x_col]], data_sub[[y_col]], h = c(600, 600), n = 100)
    
    # Interpolate the density at the given x and y coordinates
    density_values <- interp.surface(list(x = density$x, y = density$y, z = density$z),
                                     cbind(data[[x_col]], data[[y_col]]))
    density_df <- data.frame(X = data[[x_col]], Y = data[[y_col]])
    density_df[[col]] <- density_values
    density_df[[col]][is.na(density_df[[col]])] <- 0  # Handle NA values
    return(density_df)
  } else {
    density_df <- data.frame(X = data[[x_col]], Y = data[[y_col]])
    density_df[[col]] <- 0
    return(density_df)
  }
}


density_list <- data_wide
i=1

#interaction_matrix=interaction_matrix[,-1]
lig_rec_cols=colnames(data_wide)[4:166]

#data_wide=as.data.frame(data_wide)
data_wide[is.na(data_wide)==T]=0

# Iterate over each ligand-receptor pair and compute the density matrix
for (i in 1:length(lig_rec_cols)) {
  if(sum(data_wide[,lig_rec_cols[i]])>4){
    density_df <- compute_density_exact(data_wide, lig_rec_cols[i], "avg_grid_x", "avg_grid_y")
    density_list[,i+3]=density_df[,3]
    i=i+1
  }else{
    density_list[,i+3]=0
    i=i+1
  }
  print(i)
}



lig_rec_cols=colnames(density_list)[4:166]

density_list_final=density_list#[order(density_list$ID),]
#colnames(density_list_final)[1:51]=lig_rec_cols

density_list_final$Adjusted_X=density_list$avg_grid_x
density_list_final$Adjusted_Y=density_list$avg_grid_y
#density_list_final$Group=interaction_matrix$Group



data_wide

ggplot(density_list_final, aes(x = Adjusted_X, y = Adjusted_Y, color = as.factor(ifelse(density_list_final$IGF2_INSR>median(density_list_final$IGF2_INSR),1,0)))) +
  geom_point(size = 0.1)  +
  theme_minimal()

ggplot(data_wide, aes(x = avg_grid_x, y = avg_grid_y, color = as.factor(ifelse(data_wide$IGF2_INSR>0,1,0)))) +
  geom_point(size = 0.1)  +
  theme_minimal()


# 3. Perform k-means clustering with an appropriate number of clusters (e.g., k = 3)
set.seed(123)
k <- 20  # Adjust based on the elbow method plot
library(irlba)
density_list_cluster=density_list[,4:166]
pca_result <- prcomp(as.matrix(density_list_cluster[,which(colSums(density_list_cluster)>0)]), scale. = TRUE)
#pca_result <- irlba(as.matrix(density_list_cluster[,which(colSums(density_list_cluster)>0)]), nv = 10, scale. = T)
pca_df <- pca_result$u[, 1:160] * pca_result$d[1:160]


# Create a data frame with PCA results and cluster assignments
pca_df <- as.data.frame(pca_result$x)


kmeans_result <- kmeans(pca_df[,c(1:41)], centers = k, nstart = 25)

pca_df=as.data.frame(pca_df)


pca_df$Cluster <- as.factor(kmeans_result$cluster)
colnames(pca_df)[c(1:2)]=c("PC1","PC2")
# Plot the PCA results with cluster assignments
ggplot(pca_df, aes(x = PC1, y = PC2, color = pca_df$Cluster)) +
  geom_point(size=0.1) +
  labs(title = "PCA Plot of Clusters",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()

density_list_sub=density_list
ggplot(density_list_sub, aes(x = avg_grid_x, y = avg_grid_y, color = as.factor(pca_df$Cluster))) +
  geom_point(size = 1)  +
  theme_minimal()


data_cluster=data.frame(grid_id=density_list_final$grid_id,Cluster=kmeans_result$cluster)
data_cluster=data_cluster[unique(data_cluster$grid_id),]

data_final_sub=merge(data_cluster,data1_sub_plot,"grid_id")
data_final_sub=data_final_sub[unique(data_final_sub$ID),]

data_final_sub=data_final_sub[,c("ID","Cluster")]

interaction_matrix$ID=1:nrow(interaction_matrix)
data_final=merge(data_final_sub,interaction_matrix,"ID")


ggplot() +
  geom_segment(data = data_final, aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor,color = as.factor(t(data_final$Cluster))), na.rm = TRUE,size=0.1) +
  theme_classic(base_size = 15)+
  scale_color_manual(values = user_colors)  +
  labs(x = "X", y = "Y", title = "") + guides(color = guide_legend(override.aes = list(size = 3),title = "Cluster"))



###Mean LR for each MCIMs

library(dplyr)
# Group the data by the predicted label and calculate the mean for each column
data_plot=data.frame(TACIT=kmeans_result$cluster,density_list[,c(4:166)])
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.99))
#    mean_values_TACIT[ , -1][mean_values_TACIT[ , -1] > 0] <- 1
mean_values_TACIT=as.data.frame(mean_values_TACIT)
rownames(mean_values_TACIT)=mean_values_TACIT$TACIT
mean_values_TACIT <- as.data.frame((mean_values_TACIT[,-1]))


my.breaks <- c(seq(-2, 0, by=0.1),seq(0.1, 2, by=0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), colorRampPalette(colors = c("white", "red"))(length(my.breaks)/2))

aa=scale(mean_values_TACIT)
aa <- as.data.frame(aa)

aa <- aa %>%
  select_if(~ !all(is.na(.)))

library(pheatmap)
pheatmap(aa,
         cluster_cols = T,
         cluster_rows = T,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15,clustering_method="ward.D2")







library(Seurat)
data = data_wide[,colnames(count)]
orig_values = as.matrix(data)  # Consider only marker data
rownames(orig_values) = 1:nrow(data)
orig_values = t(orig_values)
orig_values_metadata = data.frame("CellID" = 1:nrow(data))
rownames(orig_values_metadata) = orig_values_metadata$CellID

scfp = CreateSeuratObject(counts = orig_values, meta.data = orig_values_metadata)

scfp = NormalizeData(scfp, normalization.method = "CLR", margin = 2)

scfp = FindVariableFeatures(scfp, selection.method = "vst", nfeatures = 100)

scfp = ScaleData(scfp, features = rownames(scfp))

#scfp = RunPCA(scfp, features = VariableFeatures(object = scfp))
#scfp = RunUMAP(scfp, reduction = "pca", dims = 1:10,metric = "correlation")

#scfp = FindNeighbors(scfp, dims = 1:10)

#scfp = FindClusters(scfp, resolution = 1.2)

Idents(scfp)=kmeans_result$cluster


scfp.markers = FindAllMarkers(scfp, only.pos = TRUE)

top5 = scfp.markers %>%
  group_by(cluster) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  as.data.frame()

DoHeatmap(scfp, features = top5$gene) + NoLegend()



top5$gene <- gsub("-", "_", top5$gene)


data_plot=data.frame(Cluster=data_final$Cluster,X_cell=data_final$X_cell,Y_cell=data_final$Y_cell,
                     X_neighbor=data_final$X_neighbor,Y_neighbor=data_final$Y_neighbor,data_final[,top5$gene])



# Define colors
user_colors2 <- c("0" = "grey90", "1" = "red")



library(dplyr)
# Loop through each cluster and generate the plot
for (cluster in unique(data_plot$Cluster)) {
  # Subset data for the specific cluster
  data_final_sub <- data_plot %>% filter(Cluster == cluster)
  
  # Calculate the sums of the columns and select the top 6
  final_results_sub=top5[which(top5$cluster==cluster),]
  top_columns <- final_results_sub$gene
  
  data_final_sub=as.data.frame(data_final_sub)
  
  
  melted_data=data_final_sub[,c("Cluster", "X_cell", "Y_cell", "X_neighbor", "Y_neighbor",top_columns)]
  # Melt the data for ggplot
  melted_data <- melted_data%>%
    pivot_longer(cols = top_columns, names_to = "Gene", values_to = "Value")
  
  # Generate the plot
  p <- ggplot(melted_data) +
    geom_segment(aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor, color = as.factor(Value)), 
                 na.rm = TRUE, size = 0.1) +
    facet_wrap(~ Gene) +
    theme_classic(base_size = 15) +
    scale_color_manual(values = user_colors2) +
    labs(x = "X", y = "Y", title = paste("Cluster", cluster)) +
    guides(color = guide_legend(override.aes = list(size = 3), title = "Value"))
  
  # Save the plot to a file
  ggsave(filename = paste0("Pariotid_Fibroblast/","cluster_all_top_", cluster, ".png"), plot = p, width = 14, height = 8)
}













































