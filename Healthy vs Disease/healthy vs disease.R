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

#---------------------------------------------

#cell_type=c("Fibroblasts")

# Make change here

# Read data

ligands_list <- read_csv("C:/Users/huynhk4/Downloads/merscope_LR.csv")


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

data15<- read_csv("TACIT_annotate_MERSCOPE/data_final.csv") #Import datasets
#data15$Group="Submandibular53_56"
data15$Group_MG="Disease"

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
datasets <- list(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14)

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








# Filter data by cell type
#data1 <- data1[data1$TACIT %in% "Fibroblasts",]

# Find common columns between ligands/receptors and data1
common <- intersect(c(ligands_list$Ligante, ligands_list$Receptores), colnames(data1))



user_colors <- c(
  "Arterioles"="darkgrey",
  "B Cells"="#FFA500",
  "Capillaries"="#aa6e28",
  "CD4+ T Cells"="#FF0000",
  "T Helper"="#FF0000",
  "CD8+ Effector T Cells"="#CC0000",
  "CD8+ Exhausted T Cells"="#FF6347",
  "CD8+ T Cells exhausted"="#FF6347",
  "CD8+ T Cells"="#FF6347",
  "Dendritic Cells"="#FFD700",
  "Ductal Cells"="#00FF00",
  "Ductal Progenitors"="#008000",
  "Ductal Proliferating"="#008000",
  "Fibroblasts"="#f032e6",
  "Fibroblast Progen"="#f032e6",
  "IgA Plasma Cells"="#7B68EE",
  "IgG Plasma Cells"="purple",
  "Intermediate Epithelium"="powderblue",
  "Ionocytes"="lightgreen",
  "M1 Macrophages"="yellow3",
  "M2 Macrophages"="gold3",
  "Macrophage"="gold3",
  "Mucous Acinar Cells"="cyan",
  "Acinar Cells"="cyan",
  "Myoepithelium"="blue",
  "Myoepitelial Cells"="blue",
  "NK Cells"="darkred",
  "Adipocyte"="#FFC0CB",
  "Memory T cell"="brown3",
  "Regulatory T Cells"="brown2",
  "Treg"="brown2",
  "Seromucous Acinar Cells"="royalblue2",
  "Smooth Muscle"="azure3",
  "T Cell Progenitors"="brown3",
  "Venules"="orange3",
  "VEC"="orange3",
  "VEC Progen"="orange1",
  "LECs"="orange",
  "Others"="grey90",
  "Mast Cell"="gold",
  "Lymphoid"="darkcyan",
  "1"="green",
  "2"="red",
  "3"="blue",
  "4"="gold",
  "5"="darkcyan",
  "6"="brown3",
  "7"="royalblue2",
  "8"="cyan",
  "9"="lightgreen",
  "10"="yellow3",
  "11"="purple",
  "12"="powderblue",
  "13"="azure3",
  "14"="royalblue1",
  "15"="cyan3",
  "16"="#FFC0CB",
  "17"="yellow1",
  "18"="purple4",
  "19"="darkred",
  "20"="#f032e6"
)

data1=data1[!(data1$TACIT%in%c("Basal Keratincytes","Suprabasal Keratinocytes","Ductal Epithelial Cells","Acinar Cells")),]

ggplot(data1, aes(x = X, y = Y, color = (TACIT))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = user_colors) +
  theme_classic(base_size = 15)+ guides(color = guide_legend(override.aes = list(size = 3),title = "Cell type"))



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


# ligands_list=ligands_list[which(ligands_list$Ligante%in%c("CXCL14",
#                                                           "IL18",
#                                                           "CXCL17",
#                                                           "CXCL5",
#
#                                                           "CXCL6",
#                                                           "CXCL13")),]



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



interaction_matrix=interaction_matrix[which(rowSums(interaction_matrix[,3:53])>0),]

count=interaction_matrix[,3:53]

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
grid_size <- 20000
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




process_chunk <- function(data_chunk, count_columns) {
  pivot_longer(
    data_chunk,
    cols = all_of(count_columns),
    names_to = "interaction_type",
    values_to = "count"
  )
}

split_data_into_chunks <- function(data, chunk_size) {
  split(data, (seq(nrow(data)) - 1) %/% chunk_size)
}

chunk_size <- 20000  # Define an appropriate chunk size
data_chunks <- split_data_into_chunks(data1_sub_plot, chunk_size)

# Assuming `count` is a vector of column names you want to pivot
count_columns <- colnames(count)

# Apply the process_chunk function to each chunk and combine results
data_long_list <- lapply(data_chunks, process_chunk, count_columns = count_columns)
data_long <- bind_rows(data_long_list)













data_long <- pivot_longer(
  data1_sub_plot,
  cols = colnames(count),  # Assuming you want to calculate density for columns starting with CCL2
  names_to = "interaction_type",
  values_to = "count"
)

# Aggregate data and calculate density
grid_data <- data_long %>%
  group_by(grid_id, grid_x, grid_y, interaction_type) %>%
  summarise(
    Total = sum(count),
    Density = Total,  # Adjust units as necessary
    .groups = 'drop'
  )
# Pivot back to wide format if needed
data_wide <- pivot_wider(
  grid_data,
  names_from = interaction_type,
  values_from = Density,
  names_prefix = ""
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
    density <- kde2d(data_sub[[x_col]], data_sub[[y_col]], h = c(10, 10), n = 200)
    
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
lig_rec_cols=colnames(data_wide)[4:54]

#data_wide=as.data.frame(data_wide)
data_wide[is.na(data_wide)==T]=0
#density_list=as.data.frame(density_list)

#952,1143

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



lig_rec_cols=colnames(density_list)[4:54]

density_list_final=density_list#[order(density_list$ID),]
#colnames(density_list_final)[1:51]=lig_rec_cols

density_list_final$Adjusted_X=density_list$avg_grid_x
density_list_final$Adjusted_Y=density_list$avg_grid_y
#density_list_final$Group=interaction_matrix$Group



data_wide



# 3. Perform k-means clustering with an appropriate number of clusters (e.g., k = 3)
set.seed(123)
k <- 20  # Adjust based on the elbow method plot
library(irlba)
density_list_cluster=density_list[,4:54]
library(Seurat)
library(reticulate)

reticulate::virtualenv_create("r-reticulate")
reticulate::use_virtualenv("r-reticulate", required = TRUE)
# Install leidenalg within the virtual environment
reticulate::py_install("leidenalg", envname = "r-reticulate")
reticulate::py_install("pandas", envname = "r-reticulate")
# Test if leidenalg can be loaded
py <- import("leidenalg")
print(py)

# Create a Seurat object from your data frame
seurat_object <- CreateSeuratObject(counts = t(as.matrix(density_list_cluster)), project = "ClusterAnalysis")

# Normalize the data
seurat_object <- NormalizeData(seurat_object)

# Find variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_object <- ScaleData(seurat_object)

# Run PCA for dimensionality reduction
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Optionally, run UMAP or t-SNE
seurat_object <- RunUMAP(seurat_object, dims = 1:10)
# Alternatively, for t-SNE:
# seurat_object <- RunTSNE(seurat_object, dims = 1:10)
# Find neighbors
seurat_object <- FindNeighbors(seurat_object, dims = 1:5)


# Use Leiden algorithm for clustering
seurat_object <- FindClusters(seurat_object, resolution = 0.005)  # Algorithm 4 is Leiden






pca_result <- prcomp(as.matrix(density_list_cluster[,which(colSums(density_list_cluster)>0)]), scale. = TRUE)
#pca_result <- irlba(as.matrix(density_list_cluster[,which(colSums(density_list_cluster)>0)]), nv = 10, scale. = T)
#pca_df <- pca_result$u[, 1:40] * pca_result$d[1:40]


# Create a data frame with PCA results and cluster assignments
pca_df <- as.data.frame(pca_result$x)


#kmeans_result <- kmeans(pca_df[,c(1:41)], centers = k, nstart = 25)

pca_df=as.data.frame(pca_df)


pca_df$Cluster <- as.factor(seurat_object$seurat_clusters)
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


data_cluster=data.frame(grid_id=density_list_final$grid_id,Cluster=seurat_object$seurat_clusters)
data_cluster=data_cluster[unique(data_cluster$grid_id),]

data_final_sub=merge(data_cluster,data1_sub_plot,"grid_id")
data_final_sub=data_final_sub[unique(data_final_sub$ID),]

data_final_sub=data_final_sub[,c("ID","Cluster")]

interaction_matrix$ID=1:nrow(interaction_matrix)
data_final=merge(data_final_sub,interaction_matrix,"ID")


#data_final_sub=data_final[which(data_final$Cluster%in%c("1","3","20")),]


# Assuming 'user_colors' is a named vector with colors corresponding to each cluster
library(ggplot2)

# Correct ggplot with geom_segment
ggplot() +
  geom_segment(data = data_final, 
               aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor, 
                   color = as.factor(Cluster)),  # Correct reference to Cluster
               na.rm = TRUE, size = 0.1) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = user_colors) +  # Ensure user_colors matches the factors in Cluster
  labs(x = "X", y = "Y", title = "") +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cluster"))



library(dplyr)
# Group the data by the predicted label and calculate the mean for each column
data_plot=data.frame(TACIT=as.numeric(seurat_object$seurat_clusters),density_list[,c(4:54)])
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5))
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
# Custom color palette from blue to white to red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the color scale
breaks <- seq(-3, 3, length.out = 101)  # Adjust to match your data range
library(pheatmap)

pheatmap(aa,
         cluster_cols = T,
         cluster_rows = T,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15,scale = "none",clustering_method="ward.D2",color = my.colors,breaks = my.breaks)

pheatmap(log10(aa+0.000000001),
         cluster_cols = T,
         cluster_rows = T,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15,scale = "none",clustering_method="ward.D2")



library(umap)   # or use library(uwot) if you installed it
library(ggplot2)
library(reticulate)
use_virtualenv("r-reticulate", required = TRUE)  # Adjust the name if your virtualenv is named differently

# Prepare the data: Exclude non-numeric columns
umap_data <-density_list[,c(4:54)]

# Using the umap package

py_install("umap-learn", envname = "r-reticulate")
# If using the uwot package
umap_result <- umap(scale(umap_data), n_neighbors = 15, min_dist = 0.1, metric = "correlation",method = "umap-learn")
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Group <- data_plot$TACIT  # Assuming you want to color points by group
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = as.factor(Group))) +
  geom_point(alpha = 0.7, size = 1.5) +
  theme_minimal() +
  labs(title = "UMAP Projection", x = "UMAP1", y = "UMAP2")  +
  scale_color_manual(values = user_colors)



# Load necessary libraries
library(dplyr)
library(ggplot2)
library(reshape2)

data_final$Group=data1$Group[data_final$cell_id]
data_final$Cluster=as.numeric(data_final$Cluster)

data_final$GroupMG=data1$Group_MG[data_final$cell_id]
data_final$GroupNICHES=ifelse(data_final$Group%in%c("Gingiva1","Gingiva1"),"Gingiva",data_final$Group)
data_final$GroupNICHES=ifelse(data_final$GroupNICHES%in%c("MSG1","MSG8","MSG9"),"MSG",data_final$GroupNICHES)
data_final$GroupNICHES=ifelse(data_final$GroupNICHES%in%c("Parotid20","Parotid35","Parotid39","Parotid40"),"Parotid",data_final$GroupNICHES)
data_final$GroupNICHES=ifelse(data_final$GroupNICHES%in%c("Tongue12","Tongue62"),"Tongue",data_final$GroupNICHES)



# List of neighbors to iterate over
neighbors <- unique(data_final$Group)

for (neighbor in neighbors) {
  
  data_final_sub = data_final[which(data_final$Group == neighbor),]
  
  library(ggplot2)
  
  # Ensure that `Cluster` is a factor in your data
  data_final_sub$Cluster <- as.factor(data_final_sub$Cluster)
  
  # Correct ggplot with geom_segment
  p=ggplot() +
    geom_segment(data = data_final_sub, 
                 aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor, 
                     color = Cluster),  # Use Cluster directly as a factor
                 na.rm = TRUE, size = 0.1) +
    theme_classic(base_size = 15) +
    scale_color_manual(values = user_colors) +  # Ensure `user_colors` aligns with Cluster levels
    labs(x = "X", y = "Y", title = "") + 
    guides(color = guide_legend(override.aes = list(size = 3), title = "Cluster"))
  
  
  
  
  # Save the plot as an SVG file
  ggsave(paste0("C:/Users/huynhk4/Downloads/OCF_MERSCOPE_0409/Healthy/n=15/slide/", neighbor, "_LR_spatial.svg"), plot = p, width = 10, height = 8, dpi = 300)
}



data_plot=data.frame(Group=data_final$Group,Cluster=data_final$Cluster)


# Load necessary library
library(ggplot2)
library(dplyr)

# Calculate proportions of each Group within each Cluster
data_plot_prop <- data_plot %>%
  group_by(Cluster, Group) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

# Define custom colors for groups based on your group labels
group_colors <- c(
  "BuccalMucosa" = "#FF9999", "Gingiva1" = "blue", "Gingiva2" = "#99FF99",
  "MSG1" = "#FFCC99", "MSG8" = "#FFB266", "MSG9" = "#FF66B2",
  "Parotid20" = "#66FFB2", "Parotid35" = "#B266FF", "Parotid39" = "#66FF66",
  "Parotid40" = "#FF66FF", "Submandibular53_56" = "#6699FF", 
  "Tongue12" = "red", "Tongue62" = "purple"
)

# Create the bar plot with custom colors
ggplot(data_plot_prop, aes(x = factor(Cluster), y = proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Cluster", y = "Proportion", title = "Proportion of Group in Each Slides in LR Cluster") +
  theme_classic(base_size = 15)+
  coord_flip()  # This flips the coordinates





data_plot=data.frame(Group=data_final$Group,Cluster=data_final$Cluster)
data_plot$Group2=ifelse(data_plot$Group%in%c("Gingiva1","Gingiva2"),"Gingiva",data_plot$Group)
data_plot$Group2=ifelse(data_plot$Group%in%c("MSG1","MSG8","MSG9"),"MSG",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group%in%c("Tongue12","Tongue62"),"Tongue",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group%in%c("Parotid20","Parotid35","Parotid39","Parotid40"),"Parotid",data_plot$Group2)

# Load necessary library
library(ggplot2)
library(dplyr)

# Calculate proportions of each Group within each Cluster
data_plot_prop <- data_plot %>%
  group_by(Cluster, Group2) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

# Define custom colors for groups based on your group labels
group_colors <- c(
  "BuccalMucosa" = "grey90", "Gingiva1" = "blue", "Gingiva" = "#99FF99",
  "MSG" = "#FFCC99", "MSG8" = "#FFB266", "MSG9" = "#FF66B2",
  "Parotid20" = "#66FFB2", "Parotid" = "#B266FF", "Parotid39" = "#66FF66",
  "Parotid40" = "#FF66FF", "Submandibular53_56" = "#6699FF", 
  "Tongue" = "lightcoral", "Tongue62" = "purple"
)

# Create the bar plot with custom colors
ggplot(data_plot_prop, aes(x = factor(Cluster), y = proportion, fill = Group2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Cluster", y = "Proportion", title = "Proportion of Group in Each Niches in LR Cluster") +
  theme_classic(base_size = 15)+
  coord_flip()  # This flips the coordinates


library(dplyr)
library(ggplot2)

# Define the desired group order for plotting
mucosal_groups <- c("Tongue", "BuccalMucosa", "Gingiva")
gland_groups <- c("MSG", "Parotid", "Submandibular53_56")

# Create a new column 'GroupOrder' to define the order
data_plot_prop <- data_plot_prop %>%
  mutate(GroupOrder = case_when(
    Group2 %in% mucosal_groups ~ 1,
    Group2 %in% gland_groups ~ 2,
    TRUE ~ 3  # If any other groups are present, they will come last
  )) %>%
  arrange(GroupOrder, Group2, desc(proportion)) %>%
  mutate(Cluster = factor(Cluster, levels = unique(Cluster)))

# Plot the data with the clusters ordered by the Group and proportion
ggplot(data_plot_prop, aes(x = Cluster, y = proportion, fill = Group2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Cluster", y = "Proportion", title = "Proportion of Group in Each Cluster") +
  theme_classic(base_size = 15) +
  coord_flip()  # This flips the coordinates




data_plot=data.frame(Group=data_final$Group,Cluster=data_final$Cluster)
data_plot$Group2=ifelse(data_plot$Group%in%c("MSG1","MSG8","MSG9","Parotid20","Parotid35","Parotid39","Parotid40","Submandibular53_56"),"Glands","Mucosal")

# Load necessary library
library(ggplot2)
library(dplyr)

# Calculate proportions of each Group within each Cluster
data_plot_prop <- data_plot %>%
  group_by(Cluster, Group2) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

# Define custom colors for groups based on your group labels
group_colors <- c(
  "BuccalMucosa" = "#FF9999", "Glands" = "blue", "Gingiva" = "#99FF99",
  "MSG" = "#FFCC99", "MSG8" = "#FFB266", "MSG9" = "#FF66B2",
  "Parotid20" = "#66FFB2", "Parotid" = "#B266FF", "Parotid39" = "#66FF66",
  "Parotid40" = "#FF66FF", "Submandibular53_56" = "#6699FF", 
  "Mucosal" = "red", "Tongue62" = "purple"
)


library(dplyr)
library(ggplot2)

# Calculate the Mucosal proportion for each Cluster and create a custom order
data_plot_prop <- data_plot_prop %>%
  group_by(Cluster) %>%
  mutate(mucosal_proportion = sum(proportion[Group2 == "Mucosal"])) %>%
  ungroup() %>%
  mutate(Cluster = factor(Cluster, levels = unique(Cluster[order(mucosal_proportion, decreasing = F)])))

# Plot the data with the clusters ordered by Mucosal group proportion
ggplot(data_plot_prop, aes(x = Cluster, y = proportion, fill = Group2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = group_colors) +
  labs(x = "Cluster", y = "Proportion", title = "Proportion of Group in Each Group in LR Cluster") +
  theme_classic(base_size = 15) +
  coord_flip()  # This flips the coordinates



data_final


# Convert data to an sf object without a CRS (non-geographic)
data_sub_sf <- st_as_sf(data_final, coords = c("X_cell", "Y_cell"), crs = NA)
data_sub_sf$cell_id <- data_final$ID

# Function to create a grid of windows
create_windows <- function(data, window_size, step_size) {
  bbox <- st_bbox(data)
  x_breaks <- seq(bbox["xmin"], bbox["xmax"] - window_size, by = step_size)
  y_breaks <- seq(bbox["ymin"], bbox["ymax"] - window_size, by = step_size)
  return(expand.grid(x = x_breaks, y = y_breaks))
}

# Function to find unique neighbors
unique_neighbors <- function(results) {
  results %>%
    distinct() %>%
    rowwise() %>%
    mutate(neighbors = list(unique(neighbors)))
}

# Function to find neighbors in a sliding window
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
distance_cutoff <- 20  # Distance cutoff for neighbors

# Function to process each group
process_group <- function(group_data) {
  # Create windows for the group
  windows <- create_windows(group_data, window_size, step_size)
  
  # Find neighbors within each window and combine results
  results <- do.call(rbind, lapply(seq_len(nrow(windows)), function(i) {
    find_neighbors_in_window(windows[i, ], group_data, window_size, distance_cutoff)
  }))
  
  # Remove duplicates and group neighbors
  results <- unique_neighbors(results)
  
  # Combine results by cell_id
  results_df <- results %>%
    group_by(cell_id) %>%
    summarise(neighbors = list(unique(unlist(neighbors))))
  
  return(results_df)
}

# Apply the process to each group and combine results
all_results <- data_sub_sf %>%
  group_by(Group) %>%
  group_map(~ process_group(.x)) %>%
  bind_rows()

# Display combined results
print(all_results)




all_results_glands_v2=all_results

results_df=all_results_glands_v2
# Filter out rows where the 'neighbors' column contains NA values
results_df <- results_df %>%
  filter(!is.na(neighbors))
table(is.na(results_df$cell_id))
results_df <- results_df %>%
  filter(!is.na(cell_id))
table(is.na(results_df$neighbors))


results_df_expanded <- results_df %>%
  unnest(neighbors)

# Check column names after unnesting
colnames(results_df_expanded)[2]="Neighbor_ID"


results_df_expanded$Neighbor_cell_id=data_final$Cluster[match(results_df_expanded$cell_id, data_final$ID)]
results_df_expanded$Neighbor_Neigbor_ID=data_final$Cluster[match(results_df_expanded$Neighbor_ID, data_final$ID)]
results_df_expanded$Group=data_final$Group[match(results_df_expanded$cell_id, data_final$ID)]


# Load the required libraries
library(dplyr)
library(igraph)

# Prepare the data
# Filter out rows where Neighbor_cell_id is equal to Neighbor_Neigbor_ID
filtered_data <- results_df_expanded %>%
  filter(Neighbor_cell_id != Neighbor_Neigbor_ID)


table(filtered_data$Group)

filtered_data=filtered_data[!(filtered_data$Group%in%c("Gingiva2","Tongue12","Tongue62")),]
#filtered_data=filtered_data[which(filtered_data$Neighbor_cell_id!=11),]
#filtered_data=filtered_data[which(filtered_data$Neighbor_Neigbor_ID!=11),]


# Count the frequency of each Neighbor_Neigbor_ID for each Neighbor_cell_id
connection_counts <- filtered_data %>%
  group_by(Neighbor_cell_id, Neighbor_Neigbor_ID) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate the total connections for each Neighbor_cell_id
total_connections <- connection_counts %>%
  group_by(Neighbor_cell_id) %>%
  summarise(total = sum(count))

# Merge the total connections back to calculate the proportion
connection_counts <- connection_counts %>%
  left_join(total_connections, by = "Neighbor_cell_id") %>%
  mutate(proportion = count / total)

# For each Neighbor_cell_id, select the top connection by proportion
top_connections <- connection_counts %>%
  group_by(Neighbor_cell_id) %>%
  slice_max(order_by = proportion, n = 3) %>%
  ungroup()

# Create the edge list for the graph
edges <- data.frame(from = top_connections$Neighbor_cell_id,
                    to = top_connections$Neighbor_Neigbor_ID,
                    weight = top_connections$proportion)

# Create the graph object
graph <- graph_from_data_frame(edges, directed = TRUE)
vertex_labels <- V(graph)$name

V(graph)$color <- user_colors[vertex_labels]
plot(graph,
     edge.width = E(graph)$weight * 10,  # Thickness of the edges proportional to the connection strength
     vertex.size = 15,                   # Size of the vertices
     vertex.label = V(graph)$name,       # Label with the cell names
     vertex.color = V(graph)$color,      # Color of the vertices
     edge.arrow.size = 0.5,              # Arrow size for directed edges
     main = "Motif Glands")




results_df_expanded$GroupNICHES=data_final$GroupNICHES[match(results_df_expanded$cell_id, data_final$ID)]


# Load the required libraries
library(dplyr)
library(igraph)

# Get unique groups
groups <- unique(results_df_expanded$GroupNICHES)

# Loop through each group and generate the graph
for (group_name in groups) {
  
  # Filter the data for the current group
  group_data <- results_df_expanded %>%
    filter(GroupNICHES == group_name)
  
  # Prepare the data: Filter out rows where Neighbor_cell_id equals Neighbor_Neigbor_ID
  filtered_data <- group_data %>%
    filter(Neighbor_cell_id != Neighbor_Neigbor_ID)
  
  #  filtered_data=filtered_data[which(filtered_data$Neighbor_cell_id!="11"),]
  #  filtered_data=filtered_data[which(filtered_data$Neighbor_Neigbor_ID!="11"),]
  
  
  # Count the frequency of each Neighbor_Neigbor_ID for each Neighbor_cell_id
  connection_counts <- filtered_data %>%
    group_by(Neighbor_cell_id, Neighbor_Neigbor_ID) %>%
    summarise(count = n()) %>%
    ungroup()
  
  # Calculate the total connections for each Neighbor_cell_id
  total_connections <- connection_counts %>%
    group_by(Neighbor_cell_id) %>%
    summarise(total = sum(count))
  
  # Merge the total connections back to calculate the proportion
  connection_counts <- connection_counts %>%
    left_join(total_connections, by = "Neighbor_cell_id") %>%
    mutate(proportion = count / total)
  
  # For each Neighbor_cell_id, select the top connection by proportion
  top_connections <- connection_counts %>%
    group_by(Neighbor_cell_id) %>%
    slice_max(order_by = proportion, n =3) %>%
    ungroup()
  
  # Create the edge list for the graph
  edges <- data.frame(from = top_connections$Neighbor_cell_id,
                      to = top_connections$Neighbor_Neigbor_ID,
                      weight = top_connections$proportion)
  
  # Create the graph object
  graph <- graph_from_data_frame(edges, directed = TRUE)
  
  
  # Ensure that the vertex labels (V(graph)$name) are present in user_colors
  vertex_labels <- V(graph)$name
  
  # Assign colors to vertices based on the user_colors vector
  V(graph)$color <- user_colors[vertex_labels]
  
  # Handle any unmatched labels with a default color (e.g., "grey" if a label is not in user_colors)
  V(graph)$color[is.na(V(graph)$color)] <- "grey"
  
  
  # Plot the graph using igraph
  plot(graph,
       edge.width = E(graph)$weight * 10,  # Thickness of the edges proportional to the connection strength
       vertex.size = 15,                   # Size of the vertices
       vertex.label = V(graph)$name,       # Label with the cell names
       vertex.color = V(graph)$color,      # Color of the vertices
       edge.arrow.size = 0.5,              # Arrow size for directed edges
       main = paste("Motif for Group:", group_name))
  
  # Optionally save the plot as an image
  filename <- paste0("C:/Users/huynhk4/Downloads/OCF_MERSCOPE_0409/Healthy/n=15/niches/OCF_motif_", group_name, ".svg")
  svg(filename)
  plot(graph,
       edge.width = E(graph)$weight * 10,  # Thickness of the edges proportional to the connection strength
       vertex.size = 15,                   # Size of the vertices
       vertex.label = V(graph)$name,       # Label with the cell names
       vertex.color = V(graph)$color,      # Color of the vertices
       edge.arrow.size = 0.5,              # Arrow size for directed edges
       main = paste("Motif for Group:", group_name))
  dev.off()
}




data_cluster=data.frame(grid_id=density_list_final$grid_id,Cluster=seurat_object$seurat_clusters)
data_cluster=data_cluster[unique(data_cluster$grid_id),]

data_final_sub=merge(data_cluster,data1_sub_plot,"grid_id")
data_final_sub=data_final_sub[unique(data_final_sub$ID),]

data_final_sub=data_final_sub[,c("ID","Cluster","grid_id")]

interaction_matrix$ID=1:nrow(interaction_matrix)
data_final_heatmap=merge(data_final_sub,interaction_matrix,"ID")

data_final_heatmap$Group=data1$Group[data_final_heatmap$cell_id]
data_final_heatmap=data_final_heatmap[,c("grid_id","Group")]




library(dplyr)
# Group the data by the predicted label and calculate the mean for each column
data_plot=data.frame(TACIT=as.numeric(seurat_object$seurat_clusters),grid_id=density_list$grid_id,density_list[,c(4:54)])
mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.5))
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
# Custom color palette from blue to white to red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the color scale
breaks <- seq(-3, 3, length.out = 101)  # Adjust to match your data range
library(pheatmap)

pheatmap(aa,
         cluster_cols = T,
         cluster_rows = T,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15,scale = "row",clustering_method="ward.D2",color = my.colors,breaks = my.breaks)




# Load necessary libraries
library(dplyr)
library(pheatmap)
library(ggplot2)

# Define custom colors and breaks for pheatmap
my.breaks <- c(seq(-2, 0, by=0.1), seq(0.1, 2, by=0.1))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("white", "red"))(length(my.breaks)/2))

# Loop through each group in data_final_heatmap
unique_groups <- unique(data_final_heatmap$Group)

# Loop through each group in data_final_heatmap
unique_groups <- unique(data_final_heatmap$Group)

library(svglite)
for (i in 13:length(unique_groups)) {
  # Filter grid_id for the current group
  grid_ids <- data_final_heatmap %>% 
    filter(Group == unique_groups[i]) %>%
    pull(grid_id)
  grid_ids=unique(grid_ids)
  # Filter data_plot using the selected grid_ids
  filtered_data <- data_plot[data_plot$grid_id %in% grid_ids, ]
  
  # Remove non-numeric columns and calculate the mean for each numeric column
  numeric_columns <- sapply(filtered_data, is.numeric)  # Identify numeric columns
  numeric_columns[c("TACIT", "grid_id")] <- FALSE      # Exclude TACIT and grid_id specifically
  
  filtered_data=filtered_data[,-2]
  mean_values_TACIT <- filtered_data %>%
    group_by(TACIT) %>%
    summarise_all(~quantile(., 0.5))
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
  # Custom color palette from blue to white to red
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # Define breaks for the color scale
  breaks <- seq(-3, 3, length.out = 101)  # Adjust to match your data range
  library(pheatmap)
  
  p=pheatmap(aa,
             cluster_cols = T,
             cluster_rows = T,
             show_colnames = TRUE,
             fontsize_col = 15,
             fontsize_row = 15,scale = "row",clustering_method="ward.D2",color = my.colors,breaks = my.breaks,
             main = paste("Heatmap for Group:", unique_groups[i]))
  
  
  file_name <- paste0("C:/Users/huynhk4/Downloads/OCF_MERSCOPE_0409/Healthy/n=20/heatmap_", unique_groups[i], ".svg")  
  
  # Use svglite to save the heatmap as an SVG
  svglite(file_name, width = 14, height = 8)
  grid::grid.draw(p$gtable)  # Draw the pheatmap gtable to the SVG device
  dev.off()  #
  
}










data1$Cell_ID=1:nrow(data1)

LR=colnames(count)

for (i in 1:length(LR)) {
  gene=strsplit(LR[i], "_")[[1]]
  # Match Cell_IDs from data1 to data_final_lung$cell_id
  gene1_col <- gene[1]  # First gene column
  gene2_col <- gene[2]  # Second gene column
  
  # Match Cell_IDs from data1 to data_final_lung$cell_id and extract gene1 column
  gene1_values <- data1[match(data_final$cell_id, data1$Cell_ID), gene1_col]
  
  # Match Neighbor_IDs from data1 to data_final_lung$Neighbor_ID and extract gene2 column
  gene2_values <- data1[match(data_final$Neighbor_ID, data1$Cell_ID), gene2_col]
  
  value=gene1_values*gene2_values
  
  data_final[,LR[i]]=value
}







library(Seurat)
data = data_final[,colnames(count)]
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

Idents(scfp)=data_final$Cluster


scfp.markers = FindAllMarkers(scfp, only.pos = TRUE)

top5 = scfp.markers %>%
  group_by(cluster)  %>%
  ungroup() %>%
  as.data.frame()







gene_pairs <- top5$gene[which(top5$cluster%in%c(15,13,12,11,8,6,7,4))]

# Split the pairs and unlist into a vector of individual genes
individual_genes <- unlist(strsplit(gene_pairs, "-"))


library(clusterProfiler)
library(org.Hs.eg.db)

genes <-unique(individual_genes)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check the conversion
print(entrez_ids)
# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa')

# View the enrichment results
head(kegg_enrichment)
# Visualize top pathways in a barplot
barplot(kegg_enrichment, showCategory = 10, title = "KEGG Pathway Enrichment for Ligand-Receptor Pairs")

# Or a dotplot
dotplot(kegg_enrichment, showCategory = 10, title = "KEGG Pathway Enrichment for Ligand-Receptor Pairs")







gene_pairs <- top5$gene[!(top5$cluster%in%c(15,13,12,11,8,6,7,4))]

# Split the pairs and unlist into a vector of individual genes
individual_genes <- unlist(strsplit(gene_pairs, "-"))


library(clusterProfiler)
library(org.Hs.eg.db)

genes <-unique(individual_genes)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Check the conversion
print(entrez_ids)
# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = 'hsa')

# View the enrichment results
head(kegg_enrichment)
# Visualize top pathways in a barplot
barplot(kegg_enrichment, showCategory = 10, title = "KEGG Pathway Enrichment for Ligand-Receptor Pairs")

# Or a dotplot
dotplot(kegg_enrichment, showCategory = 10, title = "KEGG Pathway Enrichment for Ligand-Receptor Pairs")








glands_gene_pairs <- top5$gene[which(top5$cluster%in%c(3,4,5,6,7,16,20))]
mucosal_gene_pairs <- top5$gene[which(top5$cluster%in%c(1,2,8,13,15,14,18))]



# Unique gene pairs for mucosal
unique_mucosal_gene_pairs <- setdiff(mucosal_gene_pairs, glands_gene_pairs)

# Unique gene pairs for glands
unique_glands_gene_pairs <- setdiff(glands_gene_pairs, mucosal_gene_pairs)

# Print the results
print("Unique gene pairs for mucosal:")
print(unique_mucosal_gene_pairs)

print("Unique gene pairs for glands:")
print(unique_glands_gene_pairs)


unique_gene=c(unique_glands_gene_pairs,unique_mucosal_gene_pairs)

data_final_sub=data_final[,c("Group",unique_gene)]
data_final_sub=data_final_sub[which(rowSums(data_final_sub[,-1])>0),]
data_final_sub$Group=ifelse(data_final_sub$Group%in%c("MSG1","MSG8","MSG9","Parotid20","Parotid35","Parotid39","Parotid40","Submandibular53_56"),"Glands","Mucosal")

# Load necessary libraries
library(dplyr)
library(pheatmap)

# Calculate the proportion of each column greater than 0 for each Group
data_proportions <- data_final_sub %>%
  group_by(Group) %>%
  summarise(across(everything(), ~mean(. > 0)))

# Convert the data to a matrix for heatmap plotting
data_matrix <- as.matrix(data_proportions[,-1])  # Remove the Group column for heatmap
rownames(data_matrix) <- data_proportions$Group  # Set Group as rownames

# Plot the heatmap
pheatmap(data_matrix, 
         scale = "column", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         main = "Proportion of Each Gene Pair > 0 Across Groups")






library(Seurat)
data = data_final[,colnames(count)]
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

Idents(scfp)=ifelse(data_final$Group%in%c("MSG1","MSG8","MSG9","Parotid20","Parotid35","Parotid39","Parotid40","Submandibular53_56"),"Glands","Mucosal")


# Run differential expression analysis
markers <- FindMarkers(scfp, ident.1 = "Mucosal", ident.2 = "Glands", test.use = "wilcox")

# View the top markers
head(markers)

# Add a column for significance based on adjusted p-value and logFC threshold
markers$significance <- with(markers, ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, 
                                             ifelse(avg_log2FC > 0.25, "Up in Mucosal", "Up in Glands"), 
                                             "Not Significant"))

# View the updated markers table
head(markers)

# Set a minimum threshold for the p-values to avoid extreme values in the plot
markers$p_val_adj_plot <- pmax(markers$p_val_adj, 1e-300)  # Set a lower bound for p-values

# Add the significance column again for classification if necessary
markers$significance <- with(markers, ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, 
                                             ifelse(avg_log2FC > 0.25, "Up in Mucosal", "Up in Glands"), 
                                             "Not Significant"))

# Load necessary libraries
library(ggplot2)

# Create the volcano plot with adjusted p-values
# Load necessary library
library(ggrepel)

# Add text labels for significant genes
ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj_plot), color = significance)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_text_repel(data = subset(markers, significance != "Not Significant"),
                  aes(label = rownames(subset(markers, significance != "Not Significant"))),
                  size = 4, max.overlaps = 10) +  # Adjust size and overlaps as needed
  scale_color_manual(values = c("Up in Mucosal" = "red", "Up in Glands" = "blue", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot of Differentially LR pairs", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value") +
  theme_classic(base_size = 15)

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-2, 2, length.out = 101)  # Adjust breaks for the color scale

data_final$TACIT_Ligands=data1$TACIT[data_final$cell_id]
data_final$TACIT_Receptors=data1$TACIT[data_final$Neighbor_ID]

TACIT_Ligands=data1$TACIT[interaction_matrix$cell_id]
TACIT_Receptors=data1$TACIT[interaction_matrix$Neighbor_ID]


data_plot=data.frame(TACIT_Ligands,TACIT_Receptors,count,Group=data_final$Group)
data_plot$Group=ifelse(data_plot$Group%in%c("MSG1","MSG8","MSG9","Parotid20","Parotid35","Parotid39","Parotid40","Submandibular53_56"),"Glands","Mucosal")

data_plot=data_plot[which(data_plot$Group=="Mucosal"),]
data_plot=data_plot[,-54]
# Function to calculate top interactions for each cell type
get_top_interactions <- function(cell_type, data_plot, top_n = 40) {
  filtered_data <- data_plot %>%
    filter(TACIT_Ligands == cell_type)
  
  interaction_sums <- colSums(filtered_data[,-c(1,2)])
  top_interactions <- sort(interaction_sums, decreasing = TRUE)[1:top_n]
  
  top_df <- data.frame(
    Interaction = names(top_interactions),
    Count = as.numeric(top_interactions),
    CellType = cell_type
  )
  
  return(top_df)
}

# Calculate top interactions for each cell type
cell_types <- unique(data_plot$TACIT_Ligands)
top_interactions_list <- lapply(cell_types, get_top_interactions, data_plot = data_plot)
top_interactions_df <- do.call(rbind, top_interactions_list)



# Convert to wide format for heatmap
heatmap_data <- dcast(top_interactions_df, Interaction ~ CellType, value.var = "Count", fill = 0)
rownames(heatmap_data) <- heatmap_data$Interaction
heatmap_matrix <- as.matrix(heatmap_data[,-1])

# Plot heatmap using pheatmap
library(pheatmap)

# Specify the desired order of the cell types
desired_order <- c(
  "Ductal Epithelial Cells", 
  "Basal Keratincytes", 
  "Suprabasal Keratinocytes", 
  "Ionocytes", 
  "Merkel Cells", 
  "Acinar Cells", 
  "Fibroblasts", 
  "LECs", 
  "Mural Cells", 
  "VECs", 
  "Glial/Neuron",
  "Skeletal Myocytes",
  "B Cells", 
  "T Cells",
  "CD4 T Cells", 
  "CD8 T Cells", 
  "gd T Cells", 
  "NK Cells", 
  "Dendritic Cells", 
  "Monocyte-Macrophage", 
  "Plasma Cells", 
  "Mast Cells",
  "Langerhans Cells",
  "Others"
)

# Identify the columns that are present in both `heatmap_matrix` and `desired_order`
common_columns <- intersect(desired_order, colnames(heatmap_matrix))

# Reorder `heatmap_matrix` based on `common_columns`
heatmap_matrix <- heatmap_matrix[, common_columns]


# Create the heatmap with a red-to-blue color palette
pheatmap(heatmap_matrix[which(rowSums(heatmap_matrix) > 0), ], 
         main = "Heatmap of Top Interactions for Each Cell Type in Mucosal",
         cluster_rows = F, 
         cluster_cols = F, 
         scale = "row", 
         color = color_palette,
         breaks = breaks)# Red-to-blue color gradient










data_plot=data.frame(TACIT_Ligands,TACIT_Receptors,count,Group=data_final$Group)
data_plot$Group=ifelse(data_plot$Group%in%c("MSG1","MSG8","MSG9","Parotid20","Parotid35","Parotid39","Parotid40","Submandibular53_56"),"Glands","Mucosal")

data_plot=data_plot[which(data_plot$Group=="Glands"),]
data_plot=data_plot[,-54]
# Function to calculate top interactions for each cell type
get_top_interactions <- function(cell_type, data_plot, top_n = 40) {
  filtered_data <- data_plot %>%
    filter(TACIT_Ligands == cell_type)
  
  interaction_sums <- colSums(filtered_data[,-c(1,2)])
  top_interactions <- sort(interaction_sums, decreasing = TRUE)[1:top_n]
  
  top_df <- data.frame(
    Interaction = names(top_interactions),
    Count = as.numeric(top_interactions),
    CellType = cell_type
  )
  
  return(top_df)
}

# Calculate top interactions for each cell type
cell_types <- unique(data_plot$TACIT_Ligands)
top_interactions_list <- lapply(cell_types, get_top_interactions, data_plot = data_plot)
top_interactions_df <- do.call(rbind, top_interactions_list)



# Convert to wide format for heatmap
heatmap_data <- dcast(top_interactions_df, Interaction ~ CellType, value.var = "Count", fill = 0)
rownames(heatmap_data) <- heatmap_data$Interaction
heatmap_matrix <- as.matrix(heatmap_data[,-1])

# Plot heatmap using pheatmap
library(pheatmap)
# Identify the columns that are present in both `heatmap_matrix` and `desired_order`
common_columns <- intersect(desired_order, colnames(heatmap_matrix))

# Reorder `heatmap_matrix` based on `common_columns`
heatmap_matrix <- heatmap_matrix[, common_columns]


# Create the heatmap with a red-to-blue color palette
pheatmap(heatmap_matrix[which(rowSums(heatmap_matrix) > 0), ], 
         main = "Heatmap of Top Interactions for Each Cell Type in Glands",
         cluster_rows = F, 
         cluster_cols = F, 
         scale = "row", 
         color = color_palette,
         breaks = breaks)  # Red-to-blue color gradient






data_plot=data.frame(TACIT_Ligands,TACIT_Receptors,count,Group=data_final$Group)
data_plot$Group2=ifelse(data_plot$Group%in%c("Gingiva1","Gingiva2"),"Gingiva",data_plot$Group)
data_plot$Group2=ifelse(data_plot$Group%in%c("MSG1","MSG8","MSG9"),"MSG",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group%in%c("Tongue12","Tongue62"),"Tongue",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group%in%c("Parotid20","Parotid35","Parotid39","Parotid40"),"Parotid",data_plot$Group2)

data_plot=data_plot[which(data_plot$Group2%in%c("BuccalMucosa")),]
data_plot=data_plot[,-c(55,54)]
# Function to calculate top interactions for each cell type
get_top_interactions <- function(cell_type, data_plot, top_n = 40) {
  filtered_data <- data_plot %>%
    filter(TACIT_Ligands == cell_type)
  
  interaction_sums <- colSums(filtered_data[,-c(1,2)])
  top_interactions <- sort(interaction_sums, decreasing = TRUE)[1:top_n]
  
  top_df <- data.frame(
    Interaction = names(top_interactions),
    Count = as.numeric(top_interactions),
    CellType = cell_type
  )
  
  return(top_df)
}

# Calculate top interactions for each cell type
cell_types <- unique(data_plot$TACIT_Ligands)
top_interactions_list <- lapply(cell_types, get_top_interactions, data_plot = data_plot)
top_interactions_df <- do.call(rbind, top_interactions_list)



# Convert to wide format for heatmap
heatmap_data <- dcast(top_interactions_df, Interaction ~ CellType, value.var = "Count", fill = 0)
rownames(heatmap_data) <- heatmap_data$Interaction
heatmap_matrix <- as.matrix(heatmap_data[,-1])

# Plot heatmap using pheatmap
library(pheatmap)

# Identify the columns that are present in both `heatmap_matrix` and `desired_order`
common_columns <- intersect(desired_order, colnames(heatmap_matrix))

# Reorder `heatmap_matrix` based on `common_columns`
heatmap_matrix <- heatmap_matrix[, common_columns]
# Create the heatmap with a red-to-blue color palette
pheatmap(heatmap_matrix[which(rowSums(heatmap_matrix) > 0), ], 
         main = "Heatmap of Top Interactions for Each Cell Type in BuccalMucosa",
         cluster_rows = F, 
         cluster_cols = F, 
         scale = "row", 
         color = color_palette,
         breaks = breaks)  # Specify breaks for the legend)  # Red-to-blue color gradient





data_plot=data_final#[,c(21:33,56:60)]
data_plot$Group2=ifelse(data_plot$Group%in%c("Gingiva1","Gingiva2"),"Gingiva",data_plot$Group)
data_plot$Group2=ifelse(data_plot$Group%in%c("MSG1","MSG8","MSG9"),"MSG",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group%in%c("Tongue12","Tongue62"),"Tongue",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group%in%c("Parotid20","Parotid35","Parotid39","Parotid40"),"Parotid",data_plot$Group2)
data_plot$Group2=ifelse(data_plot$Group2%in%c("BuccalMucosa","Gingiva","Tongue"),"Mucosal","Glands")
# Load necessary libraries
library(dplyr)

# Calculate the proportion of non-zero values for each ligand-receptor pair by group
ligand_receptor_signature <- data_plot %>%
  group_by(Group2) %>%
  summarise(across(colnames(count), ~ mean(. > 0)))  # Adjust if necessary for more columns

# View the results
head(ligand_receptor_signature)


# Convert to matrix for heatmap
signature_matrix <- as.matrix(ligand_receptor_signature[,-1])  # Remove the Group column for heatmap
rownames(signature_matrix) <- ligand_receptor_signature$Group2  # Set Group as row names
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-2, 2, length.out = 101)  # Adjust breaks for the color scale
# Plot heatmap
pheatmap(signature_matrix, 
         main = "Ligand-Receptor Signatures Across Groups", 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = color_palette,
         breaks = breaks)  # Specify breaks for the legend)




# Load necessary libraries
library(ggplot2)
user_colors <- c(
  "Acinar Cells"="cyan", "B Cells"="#FFA500", "Basal Keratincytes"="lightgrey",
  "CD4 T Cells"="#FF0000", "CD8 T Cells"="#FF6347", "Dendritic Cells"="#FFD700",
  "Ductal Epithelial Cells"="#00FF00", "Fibroblasts"="#f032e6", "gd T Cells"="#FF6347",
  "Glial/Neuron"="lightgrey", "Ionocytes"="lightgreen", "Langerhans Cells"="#FFD700",
  "LECs"="orange", "Mast Cells"="gold", "Merkel Cells"="lightgrey",
  "Monocyte-Macrophage"="goldenrod", "Mural Cells"="lightgrey", "Myoepithelial Cells"="blue",
  "NK Cells"="darkred", "Others"="lightgrey", "Plasma Cells"="purple", 
  "Skeletal Myocytes"="lightgrey", "Suprabasal Keratinocytes"="powderblue", "VECs"="orange3"
)
# Get unique groups
unique_groups <- unique(data1$Group)

# Loop through each group and generate the plot
for (group in unique_groups) {
  
  # Subset the data for the current group
  data_sub <- data1[which(data1$Group == group),]
  
  # Create the plot
  p <- ggplot(data_sub, aes(x = X, y = Y, color = as.factor(data_sub$TACIT))) +
    geom_point(size = 0.1) +
    scale_color_manual(values = user_colors) +
    theme_classic(base_size = 15) +
    guides(color = guide_legend(override.aes = list(size = 3), title = "Cell Type")) +
    labs(title = paste("Plot for", group))
  
  # Display the plot for the current group (optional)
  print(p)
  
  # Export the plot as an image (PNG format) for the current group
  ggsave(paste0("C:/Users/huynhk4/Downloads/OCF_MERSCOPE_0409/Healthy/n=15/slide/",group, "_plot.png"), plot = p, width = 12, height = 8, dpi = 300)
}



















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
  ggsave(filename = paste0("TACIT_annotate_MERSCOPE/Gingiva_Fibroblast/","cluster_top_", cluster, ".png"), plot = p, width = 14, height = 8)
}














































data_final$TACIT_Ligands=data1$TACIT[data_final$cell_id]
data_final$TACIT_Receptors=data1$TACIT[data_final$Neighbor_ID]


tail(colnames(interaction_matrix))

vec=colSums(interaction_matrix[,3:166])



df <- data.frame(ColumnSum = vec)


# Plot the density plot using ggplot2
ggplot(df, aes(x = ColumnSum)) +
  geom_density(color = "blue", fill = "skyblue", alpha = 0.5) +
  labs(
    title = "Density Plot of LR",
    x = "Sum of each LR",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )+xlim(0,2000)


library(dplyr)
library(reshape2)

data_plot=data.frame(Cluster=data_final$Cluster,Ligands=data_final$TACIT_Ligands,Receptors=data_final$TACIT_Receptors)

# Calculate the proportion of Ligands for each cluster
ligands_proportion <- data_plot %>%
  group_by(Cluster, Ligands) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Normalize proportions by Ligands (cell type)
ligands_proportion <- ligands_proportion %>%
  group_by(Ligands) %>%
  mutate(Normalized_Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

# Convert to wide format for pheatmap
ligands_matrix <- dcast(ligands_proportion, Ligands ~ Cluster, value.var = "Normalized_Proportion", fill = 0)
ligands_matrix <- as.matrix(ligands_matrix[,-1])  # Remove the Ligands column to get only the matrix

# Set row names to Ligands
rownames(ligands_matrix) <- ligands_proportion$Ligands[!duplicated(ligands_proportion$Ligands)]

# Plot heatmap for Ligands using pheatmap
pheatmap(ligands_matrix, 
         main = "Heatmap of Normalized Ligands Proportion within Each Cluster",
         cluster_rows = TRUE, 
         cluster_cols = TRUE)



# Calculate the proportion of Receptors for each cluster
receptors_proportion <- data_plot %>%
  group_by(Cluster, Receptors) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Normalize proportions by Receptors (cell type)
receptors_proportion <- receptors_proportion %>%
  group_by(Receptors) %>%
  mutate(Normalized_Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

# Convert to wide format for pheatmap
receptors_matrix <- dcast(receptors_proportion, Receptors ~ Cluster, value.var = "Normalized_Proportion", fill = 0)
receptors_matrix <- as.matrix(receptors_matrix[,-1])  # Remove the Receptors column to get only the matrix

# Set row names to Receptors
rownames(receptors_matrix) <- receptors_proportion$Receptors[!duplicated(receptors_proportion$Receptors)]

# Plot heatmap for Receptors using pheatmap
pheatmap(receptors_matrix, 
         main = "Heatmap of Normalized Receptors Proportion within Each Cluster",
         cluster_rows = TRUE, 
         cluster_cols = TRUE)





data_final_sub=data_final[which(data_final$Cluster==16),]

ggplot() +
  geom_point(data = data_final,size = 0.5,aes(x=X_cell,Y_cell,color=TACIT_Ligands))+
  geom_point(data = data_final,size = 0.5,aes(x=X_neighbor,Y_neighbor,color=TACIT_Receptors))+
  geom_segment(data = data_final_sub, aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor,color = as.factor(t(data_final_sub$Cluster))), na.rm = TRUE,size=0.1) +
  theme_classic(base_size = 15)+
  scale_color_manual(values = user_colors)  +
  labs(x = "X", y = "Y", title = "") + guides(color = guide_legend(override.aes = list(size = 3),title = "Cluster"))





data_plot=data.frame(Cluster=data_final$Cluster,Ligands=data_final$TACIT_Ligands,Receptors=data_final$TACIT_Receptors)



data_plot=data_plot[which(data_plot$Ligands=="B Cells"),]

interaction_matrix_CX3CL1_CX3CR1_plot=data_plot[,c("Ligands","Receptors")]
# Calculate interaction counts
interaction_counts <- interaction_matrix_CX3CL1_CX3CR1_plot %>%
  group_by(Ligands, Receptors) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(transformed_percentile = log10(count + 1))


top_5_interaction_counts <- interaction_counts %>%
  slice_max(order_by = count, n = 20)

# Create a graph from the interaction data
g <- graph_from_data_frame(top_5_interaction_counts, directed = TRUE)
# Create a color palette
# Create a large color palette
color_palette_large <- c(brewer.pal(8, "Set1"),     # Red, blue, green, etc.
                         brewer.pal(8, "Set2"),     # Lighter and diverse colors
                         brewer.pal(8, "Set3"),     # Mixed colors
                         brewer.pal(8, "Dark2"),    # Dark colors
                         brewer.pal(8, "Pastel1"),  # Pastel colors
                         brewer.pal(8, "Pastel2"),  # More pastel colors
                         "red", "green", "blue", "orange", "purple", "cyan", "magenta")  # Specific colors

# Function to get a specific number of random colors from the palette
get_random_colors <- function(number_of_colors) {
  if (number_of_colors > length(color_palette_large)) {
    stop("Requested number of colors exceeds the prepared palette size.")
  }
  sample(color_palette_large, number_of_colors)
}
unique_types <- unique(c(data_final$TACIT_Ligands, data_final$TACIT_Receptors))

# Example usage: Get 10 random colors
n <- length(unique_types)  # Assume you have 10 types this time
color_palette <- get_random_colors(n)
names(color_palette) <- unique_types
# Assign colors to vertices based on cell type
V(g)$color <- color_palette[V(g)$name]
# Assign colors to edges based on interaction strength
# Use a continuous color scale or manually define breaks for clarity
E(g)$color <- ifelse(E(g)$transformed_percentile > median(E(g)$transformed_percentile, na.rm = TRUE), "grey90", "grey90")


#---------------------------------------------

# Graph cell type between ligands and receptors

# Plotting
plot(g, layout = layout_in_circle(g),
     edge.width = E(g)$transformed_percentile * 5,  # Scale edge width
     edge.arrow.size = 0.5,
     vertex.size = 15,
     vertex.label.cex = 2,
     vertex.color = V(g)$color,
     edge.color = E(g)$color,
     main = "Interaction Network")





interaction_matrix_pair=count[which(kmeans_result$cluster==4),]

# Convert the interaction matrix to a data frame of connection counts
interaction_counts_v2 <- interaction_matrix_pair %>%
  summarise_all(sum) %>% # Sum the connections for each ligand-receptor pair
  pivot_longer(everything(), names_to = "ligand_receptor", values_to = "count") %>%
  separate(ligand_receptor, into = c("ligand", "receptor"), sep = "_") %>%
  filter(count > 0) # Optional: filter out pairs with zero connections

# Assuming cell types are somehow encoded or need to be mapped:
# You will need additional data or logic to map ligands and receptors to specific cell types
# if you want to incorporate 'cell_type_ligands' and 'cell_type_receptors'
# For now, this code simply summarizes counts of connections for each ligand-receptor pair

# Optionally, add percentile transformation if needed:
interaction_counts_v2 <- interaction_counts_v2 %>%
  mutate(transformed_percentile = log10(count + 1))


#---------------------------------------------

# Top ligands and receptor (n=5)

top_5_interaction_counts <- interaction_counts_v2 %>%
  slice_max(order_by = count, n = 5)


# Create a graph from the interaction data
g <- graph_from_data_frame(top_5_interaction_counts, directed = TRUE)
# Create a color palette
# Create a large color palette
color_palette_large <- c(brewer.pal(8, "Set1"),     # Red, blue, green, etc.
                         brewer.pal(8, "Set2"),     # Lighter and diverse colors
                         brewer.pal(8, "Set3"),     # Mixed colors
                         brewer.pal(8, "Dark2"),    # Dark colors
                         brewer.pal(8, "Pastel1"),  # Pastel colors
                         brewer.pal(8, "Pastel2"),  # More pastel colors
                         "red", "green", "blue", "orange", "purple", "cyan", "magenta")  # Specific colors

# Function to get a specific number of random colors from the palette
get_random_colors <- function(number_of_colors) {
  if (number_of_colors > length(color_palette_large)) {
    stop("Requested number of colors exceeds the prepared palette size.")
  }
  sample(color_palette_large, number_of_colors)
}
unique_types <- unique(c(data_final$TACIT_Ligands, data_final$TACIT_Receptors))

# Example usage: Get 10 random colors
n <- length(unique_types)  # Assume you have 10 types this time
color_palette <- get_random_colors(n)
names(color_palette) <- unique_types
# Assign colors to vertices based on cell type
V(g)$color <- color_palette[V(g)$name]
# Assign colors to edges based on interaction strength
# Use a continuous color scale or manually define breaks for clarity
E(g)$color <- ifelse(E(g)$transformed_percentile > median(E(g)$transformed_percentile, na.rm = TRUE), "grey90", "grey90")


#---------------------------------------------

# Graph cell type between ligands and receptors

# Plotting
plot(g, layout = layout_in_circle(g),
     edge.width = E(g)$transformed_percentile * 5,  # Scale edge width
     edge.arrow.size = 0.1,
     vertex.size = 15,
     vertex.label.cex = 0.8,
     vertex.color = V(g)$color,
     edge.color = E(g)$color,
     main = "Top ligands and receptors")






TACIT_Ligands=data1$TACIT[interaction_matrix$cell_id]
TACIT_Receptors=data1$TACIT[interaction_matrix$Neighbor_ID]
data_plot=data.frame(TACIT_Ligands,TACIT_Receptors,count)









# Function to calculate top interactions for each cell type
get_top_interactions <- function(cell_type, data_plot, top_n = 40) {
  filtered_data <- data_plot %>%
    filter(TACIT_Ligands == cell_type)
  
  interaction_sums <- colSums(filtered_data[,-c(1,2)])
  top_interactions <- sort(interaction_sums, decreasing = TRUE)[1:top_n]
  
  top_df <- data.frame(
    Interaction = names(top_interactions),
    Count = as.numeric(top_interactions),
    CellType = cell_type
  )
  
  return(top_df)
}

# Calculate top interactions for each cell type
cell_types <- unique(data_plot$TACIT_Ligands)
top_interactions_list <- lapply(cell_types, get_top_interactions, data_plot = data_plot)
top_interactions_df <- do.call(rbind, top_interactions_list)



# Plot the top interactions using ggplot2
ggplot(top_interactions_df, aes(x = reorder(Interaction, -Count), y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Top Interactions for Each Cell Type",
    x = "Interaction",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = user_colors) # Customize colors as needed


# Convert to wide format for heatmap
heatmap_data <- dcast(top_interactions_df, Interaction ~ CellType, value.var = "Count", fill = 0)
rownames(heatmap_data) <- heatmap_data$Interaction
heatmap_matrix <- as.matrix(heatmap_data[,-1])

# Plot heatmap using pheatmap
library(pheatmap)
pheatmap(heatmap_matrix[which(rowSums(heatmap_matrix)>0),], 
         main = "Heatmap of Top Interactions for Each Cell Type",
         cluster_rows = TRUE, 
         cluster_cols = TRUE,scale="row")












data_plot=data.frame(Cluster=data_final$Cluster,X_cell=data_final$X_cell,Y_cell=data_final$Y_cell,
                     X_neighbor=data_final$X_neighbor,Y_neighbor=data_final$Y_neighbor,data_final[,colnames(count)])



# Define colors
user_colors2 <- c("0" = "grey90", "1" = "red")



library(dplyr)
# Loop through each cluster and generate the plot
for (cluster in unique(data_plot$Cluster)) {
  # Subset data for the specific cluster
  data_final_sub <- data_plot %>% filter(Cluster == cluster)
  
  # Calculate the sums of the columns and select the top 6
  column_sums <- colSums(data_final_sub[, colnames(count)])
  top_columns <- names(sort(column_sums, decreasing = TRUE))[1:6]
  
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
  ggsave(filename = paste0("cluster_xenium/","cluster_top_", cluster, ".png"), plot = p, width = 14, height = 8)
}






# Load necessary libraries
library(dplyr)
library(tidyr)
library(circlize)
library(pheatmap)

# Create the data_plot data frame
data_plot <- data.frame(TACIT_Ligands, TACIT_Receptors, Cluster = data_final$Cluster, data_final[, colnames(count)])

# Generate a proportion table for heatmap
proportion_table <- prop.table(table(data_plot$TACIT_Receptors, data_plot$Cluster), margin = 2)
pheatmap(proportion_table, scale = "row")

# Filter the data for Cluster 17
data_plot <- data_plot %>% filter(Cluster == 17)
final_results_sub <- top5 %>% filter(cluster == 17)

# Select the ligand-receptor columns
ligand_receptor_columns <- final_results_sub$gene
filtered_data <- data_plot

# Gather the ligand-receptor pairs and calculate their frequency
ligand_receptor_pairs <- filtered_data[,ligand_receptor_columns] %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "Ligand_Receptor", values_to = "Frequency") %>%
  arrange(desc(Frequency))

# Select the top 10 ligand-receptor pairs
top10_ligand_receptor_pairs <- head(ligand_receptor_pairs, 10)

# Prepare data for circos plot
circos_data <- filtered_data %>%
  gather(key = "Ligand_Receptor", value = "Value", -TACIT_Ligands, -TACIT_Receptors) %>%
  filter(Value > 0 & Ligand_Receptor != "Cluster")

# Summarize the number of connections for each ligand-receptor pair between cell types
summarized_connections <- circos_data %>%
  group_by(TACIT_Ligands, TACIT_Receptors, Ligand_Receptor) %>%
  summarise(Value = sum(Value), .groups = 'drop') %>%
  filter(Ligand_Receptor %in% top10_ligand_receptor_pairs$Ligand_Receptor) %>%
  group_by(Ligand_Receptor) %>%
  mutate(Proportion = Value / sum(Value)) %>%
  ungroup()

# Keep only the top 5 connections for each Ligand_Receptor
summarized_connections <- summarized_connections %>%
  group_by(Ligand_Receptor) %>%
  slice_max(order_by = Value, n = 10) %>%
  ungroup()



summarized_connections <- summarized_connections %>%
  group_by(TACIT_Ligands) %>%
  mutate(Total_Value_Ligands = sum(Value)) %>%
  ungroup() %>%
  mutate(Proportion_from_TACIT_Ligands = Value / Total_Value_Ligands)

summarized_connections <- summarized_connections %>%
  group_by(TACIT_Receptors) %>%
  mutate(Total_Value_Receptors = sum(Value)) %>%
  ungroup() %>%
  mutate(Proportion_from_TACIT_Receptors = Value / Total_Value_Receptors)


top_cell_types=data.frame(TACIT_Ligands=names(table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors)))),
                          count=as.numeric((table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors))))),
                          proportion=as.numeric((table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors)))))/sum(as.numeric((table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors)))))))


# Reassign cell types with proportions less than 0.02 to "Other"
top_cell_types$TACIT_Ligands[top_cell_types$proportion < 0.05] <- "Other"

# Sum up the counts and proportions for "Other"
other_count <- sum(top_cell_types$count[top_cell_types$TACIT_Ligands == "Other"])
other_proportion <- sum(top_cell_types$proportion[top_cell_types$TACIT_Ligands == "Other"])

# Remove individual rows with "Other" and add a single row
top_cell_types <- top_cell_types %>%
  filter(TACIT_Ligands != "Other") %>%
  bind_rows(tibble(TACIT_Ligands = "Other", count = other_count, proportion = other_proportion))




summarized_connections$TACIT_Ligands=ifelse(summarized_connections$TACIT_Ligands%in%top_cell_types$TACIT_Ligands,summarized_connections$TACIT_Ligands,"Other")
summarized_connections$TACIT_Receptors=ifelse(summarized_connections$TACIT_Receptors%in%top_cell_types$TACIT_Ligands,summarized_connections$TACIT_Receptors,"Other")



# Prepare the sector widths for circos plot
sector_widths <- top_cell_types %>%
  mutate(sector = TACIT_Ligands)

# Initialize circos plot with sectors and their widths
circos.clear()
circos.par(gap.degree = 2)
circos.initialize(factors = sector_widths$TACIT_Ligands, xlim = c(0, 1), sector.width = sector_widths$proportion)

# Add segments for cell types with background colors to form a circle
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.1, bg.border = "black", bg.col = user_colors[sector_widths$TACIT_Ligands], panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.5, CELL_META$sector.index, cex = 0.8, facing = "clockwise")
})

# Create a unique color for each ligand-receptor pair
ligand_receptor_colors <- setNames(rainbow(length(unique(summarized_connections$Ligand_Receptor))),
                                   unique(summarized_connections$Ligand_Receptor))

# Draw the links between TACIT_Ligands and TACIT_Receptors
unique_ligands <- unique(summarized_connections$TACIT_Ligands)
for (ligand in unique_ligands) {
  ligand_data <- summarized_connections %>% filter(TACIT_Ligands == ligand)
  for (i in 1:nrow(ligand_data)) {
    receptor <- ligand_data$TACIT_Receptors[i]
    lig_rec <- ligand_data$Ligand_Receptor[i]
    color <- ligand_receptor_colors[lig_rec]
    
    width_L <- ligand_data$Proportion[i] / sum(ligand_data$Proportion)
    receptor_data <- summarized_connections %>% filter(TACIT_Receptors == receptor)
    width_R <- ligand_data$Proportion[i] / sum(receptor_data$Proportion)
    
    circos.link(sector.index1 = ligand_data$TACIT_Ligands[i],
                point1 = c((i-1)/nrow(ligand_data), i/nrow(ligand_data)),
                sector.index2 = ligand_data$TACIT_Receptors[i],
                point2 = c((i-1)/nrow(ligand_data), i/nrow(ligand_data)),
                col = color, directional = 1)
  }
}

# Add a legend
legend("topright", legend = names(ligand_receptor_colors), col = as.character(ligand_receptor_colors), pch = 19)
# Add a title
title(main = "Chord Diagram of Ligand-Receptor Interactions in Cluster 17")











TACIT_Ligands=data1$TACIT[data_final$cell_id]
TACIT_Receptors=data1$TACIT[data_final$Neighbor_ID]





# Load necessary libraries
library(dplyr)
library(tidyr)
library(circlize)
library(pheatmap)
library(ggplot2)

# Define a function to generate and save chord diagrams for each cluster
generate_chord_diagram <- function(cluster_id) {
  # Filter the data for the current cluster
  data_plot <- data.frame(TACIT_Ligands, TACIT_Receptors, Cluster = data_final$Cluster, data_final[, colnames(count)])
  data_plot <- data_plot %>% filter(Cluster == cluster_id)
  final_results_sub <- top5 %>% filter(cluster == cluster_id)
  
  # Select the ligand-receptor columns
  ligand_receptor_columns <- final_results_sub$gene
  filtered_data <- data_plot
  
  # Gather the ligand-receptor pairs and calculate their frequency
  ligand_receptor_pairs <- filtered_data[,ligand_receptor_columns] %>%
    summarise(across(everything(), sum)) %>%
    pivot_longer(cols = everything(), names_to = "Ligand_Receptor", values_to = "Frequency") %>%
    arrange(desc(Frequency))
  
  # Select the top 10 ligand-receptor pairs
  top10_ligand_receptor_pairs <- head(ligand_receptor_pairs, 5)
  
  # Prepare data for circos plot
  circos_data <- filtered_data %>%
    gather(key = "Ligand_Receptor", value = "Value", -TACIT_Ligands, -TACIT_Receptors) %>%
    filter(Value > 0 & Ligand_Receptor != "Cluster")
  
  # Summarize the number of connections for each ligand-receptor pair between cell types
  summarized_connections <- circos_data %>%
    group_by(TACIT_Ligands, TACIT_Receptors, Ligand_Receptor) %>%
    summarise(Value = sum(Value), .groups = 'drop') %>%
    filter(Ligand_Receptor %in% top10_ligand_receptor_pairs$Ligand_Receptor) %>%
    group_by(Ligand_Receptor) %>%
    mutate(Proportion = Value / sum(Value)) %>%
    ungroup()
  
  # Keep only the top 5 connections for each Ligand_Receptor
  summarized_connections <- summarized_connections %>%
    group_by(Ligand_Receptor) %>%
    slice_max(order_by = Value, n = 5) %>%
    ungroup()
  
  # Calculate proportions
  summarized_connections <- summarized_connections %>%
    group_by(TACIT_Ligands) %>%
    mutate(Total_Value_Ligands = sum(Value)) %>%
    ungroup() %>%
    mutate(Proportion_from_TACIT_Ligands = Value / Total_Value_Ligands)
  
  summarized_connections <- summarized_connections %>%
    group_by(TACIT_Receptors) %>%
    mutate(Total_Value_Receptors = sum(Value)) %>%
    ungroup() %>%
    mutate(Proportion_from_TACIT_Receptors = Value / Total_Value_Receptors)
  
  top_cell_types <- data.frame(
    TACIT_Ligands = names(table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors)))),
    count = as.numeric((table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors))))),
    proportion = as.numeric((table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors)))))/sum(as.numeric((table((c(summarized_connections$TACIT_Ligands, summarized_connections$TACIT_Receptors))))))
  )
  
  # Reassign cell types with proportions less than 0.02 to "Other"
  top_cell_types$TACIT_Ligands[top_cell_types$proportion < 0.05] <- "Other"
  
  # Sum up the counts and proportions for "Other"
  other_count <- sum(top_cell_types$count[top_cell_types$TACIT_Ligands == "Other"])
  other_proportion <- sum(top_cell_types$proportion[top_cell_types$TACIT_Ligands == "Other"])
  
  # Remove individual rows with "Other" and add a single row
  top_cell_types <- top_cell_types %>%
    filter(TACIT_Ligands != "Other") %>%
    bind_rows(tibble(TACIT_Ligands = "Other", count = other_count, proportion = other_proportion))
  
  top_cell_types=top_cell_types[which(top_cell_types$proportion>0),]
  
  
  summarized_connections$TACIT_Ligands <- ifelse(summarized_connections$TACIT_Ligands %in% top_cell_types$TACIT_Ligands, summarized_connections$TACIT_Ligands, "Other")
  summarized_connections$TACIT_Receptors <- ifelse(summarized_connections$TACIT_Receptors %in% top_cell_types$TACIT_Ligands, summarized_connections$TACIT_Receptors, "Other")
  
  # Prepare the sector widths for circos plot
  sector_widths <- top_cell_types %>%
    mutate(sector = TACIT_Ligands)
  
  
  
  #  file_name <- paste0("cluster_xenium/Chord_Diagram_Cluster_", cluster_id, ".png")
  pdf(file = paste0("TACIT_annotate_MERSCOPE/Gingiva_Fibroblast/Chord_Diagram_Cluster_", cluster_id, ".pdf"), width = 8, height = 6)
  
  
  # Initialize circos plot with sectors and their widths
  circos.clear()
  circos.par(gap.degree = 2)
  circos.initialize(factors = sector_widths$TACIT_Ligands, xlim = c(0, 1), sector.width = sector_widths$proportion)
  
  # Add segments for cell types with background colors to form a circle
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.1, bg.border = "black", bg.col = user_colors[sector_widths$TACIT_Ligands], panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.5, CELL_META$sector.index, cex = 0.8, facing = "clockwise")
  })
  
  # Create a unique color for each ligand-receptor pair
  ligand_receptor_colors <- setNames(rainbow(length(unique(summarized_connections$Ligand_Receptor))),
                                     unique(summarized_connections$Ligand_Receptor))
  
  # Draw the links between TACIT_Ligands and TACIT_Receptors
  unique_ligands <- unique(summarized_connections$TACIT_Ligands)
  for (ligand in unique_ligands) {
    ligand_data <- summarized_connections %>% filter(TACIT_Ligands == ligand)
    for (i in 1:nrow(ligand_data)) {
      receptor <- ligand_data$TACIT_Receptors[i]
      lig_rec <- ligand_data$Ligand_Receptor[i]
      color <- ligand_receptor_colors[lig_rec]
      
      width_L <- ligand_data$Proportion[i] / sum(ligand_data$Proportion)
      receptor_data <- summarized_connections %>% filter(TACIT_Receptors == receptor)
      width_R <- ligand_data$Proportion[i] / sum(receptor_data$Proportion)
      
      circos.link(sector.index1 = ligand_data$TACIT_Ligands[i],
                  point1 = c((i-1)/nrow(ligand_data), i/nrow(ligand_data)),
                  sector.index2 = ligand_data$TACIT_Receptors[i],
                  point2 = c((i-1)/nrow(ligand_data), i/nrow(ligand_data)),
                  col = color, directional = 1)
    }
  }
  
  # Add a legend
  legend("topright", legend = names(ligand_receptor_colors), col = as.character(ligand_receptor_colors), pch = 19)
  
  # Add a title
  title(main = paste("Chord Diagram of Ligand-Receptor Interactions in Cluster", cluster_id))
  
  # Save the plot
  dev.off()
}

# Get all unique clusters
all_clusters <- unique(data_final$Cluster)

# Loop over all clusters and generate the chord diagrams
for (cluster_id in all_clusters[2:20]) {
  generate_chord_diagram(cluster_id)
}




data_final_sub=data_final[which(data_final$TACIT_Ligands%in%c("T Cell Progenitors")|
                                  data_final$TACIT_Receptors%in%c("T Cell Progenitors")),]
ggplot() +
  geom_point(data = data_final_sub,size = 0.1,aes(x=X_cell,Y_cell,color=TACIT_Ligands))+
  geom_point(data = data_final_sub,size = 0.1,aes(x=X_neighbor,Y_neighbor,color=TACIT_Receptors))+
  geom_segment(data = data_final_sub, aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor,color = as.factor(t(data_final_sub$Cluster))), na.rm = TRUE,size=0.1) +
  theme_classic(base_size = 15)+
  scale_color_manual(values = user_colors)  +
  labs(x = "X", y = "Y", title = "") + guides(color = guide_legend(override.aes = list(size = 3),title = "Cluster"))




data_cluster=data.frame(grid_id=density_list_final$grid_id,Cluster=kmeans_result$cluster)
data_cluster=data_cluster[unique(data_cluster$grid_id),]

data_final_sub=merge(data_cluster,data1_sub_plot,"grid_id")
data_final_sub=data_final_sub[unique(data_final_sub$ID),]

data_final_sub=data_final_sub[,c("ID","Cluster","grid_id")]

interaction_matrix$ID=1:nrow(interaction_matrix)
data_final=merge(data_final_sub,interaction_matrix,"ID")

data_final$Group=data1$Group[data_final$cell_id]

data_final=data_final[which(data_final$Cluster==11),]

data_final$Group2=ifelse(data_final$Group%in%c("Tongue12","Tongue62","MSG8","MSG9",
                                               "Gingiva2","Parotid39","Submandibular53_56")&data_final$Cluster=="11","11_high_density",data_final$Cluster)

library(dplyr)

# Load necessary library
library(dplyr)

# Identify gene-related columns (all columns that contain an underscore)
gene_columns <- grep("_", names(data_final), value = TRUE)

# Check for and exclude 'grid_id' from gene columns if mistakenly included
gene_columns <- setdiff(gene_columns, "grid_id")

# Summarize data by grid_id: sum gene-related columns and handle Group2
summarized_data <- data_final %>%
  group_by(grid_id) %>%
  summarise(
    across(all_of(gene_columns), sum, na.rm = TRUE),  # Sum gene-related columns
    Group2 = ifelse(length(unique(Group2)) == 1, unique(Group2), first(Group2))  # Handle non-unique Group2
  )

# Print the summarized data
print(summarized_data)


# Print the summarized data
print(summarized_data)
library(dplyr)
# Group the data by the predicted label and calculate the mean for each column
data_plot=data.frame(TACIT=summarized_data$Group2,summarized_data[,c(4:54)])




mean_values_TACIT <- data_plot %>%
  group_by(TACIT) %>%
  summarise_all(~quantile(., 0.999))
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
# Custom color palette from blue to white to red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the color scale
breaks <- seq(-3, 3, length.out = 101)  # Adjust to match your data range
library(pheatmap)
aa=aa[c(3,4),]
pheatmap(aa,
         cluster_cols = T,
         cluster_rows = T,
         show_colnames = TRUE,
         fontsize_col = 15,
         fontsize_row = 15,clustering_method="ward.D2",color = my.colors,breaks = my.breaks)



library(umap)
library(ggplot2)

# Run UMAP
umap_result <- umap(data_plot[,-1], n_neighbors = 15, min_dist = 0.1, metric = "euclidean")

# Convert UMAP result to a data frame
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Add additional information, e.g., Group2
umap_df$Group2 <- data_final$Group2






data_plot=data.frame(Group=data_final$GroupMG,TACIT_Ligands=data_final$TACIT_Ligands,
                     TACIT_Receptors=data_final$TACIT_Receptors,Cluster=data_final$Cluster)

data_plot_subset=data_plot[which(data_plot$Group=="BuccalMucosa"),-1]


data_combined <- data_plot_subset %>%
  pivot_longer(cols = c(TACIT_Ligands, TACIT_Receptors), names_to = "Type", values_to = "Cell_Type") %>%
  group_by(Cluster, Cell_Type) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate the proportion of each cell type within each cluster
data_prop <- data_combined %>%
  group_by(Cluster) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

data_prop=data_prop[,-3]


data_matrix <- data_prop %>%
  pivot_wider(names_from = Cell_Type, values_from = Proportion, values_fill = list(Proportion = 0))

data_matrix=as.data.frame(data_matrix)
# Set row names to clusters for the heatmap
rownames(data_matrix) <- data_matrix$Cluster
data_matrix <- data_matrix %>% select(-Cluster)  # Remove Cluster column after setting row names

data_matrix=data_matrix[,-1]


color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the color scale
breaks <- seq(-2, 2, length.out = 101)  # Adjust to match your data range
pheatmap(as.matrix(data_matrix),
         cluster_rows = F,
         cluster_cols = TRUE,
         main = "Heatmap of Cell Type Proportions by BuccalMucosa",
         border_color = NA,
         scale = "row",color = color_palette,breaks = breaks)  # Adjust scale if necessary




# Get unique groups
groups <- unique(data_plot$Group)

# Loop through each group to generate the heatmap
for (group in groups) {
  # Subset the data for the current group
  data_plot_subset <- data_plot[data_plot$Group == group, ]
  data_plot_subset <- data_plot_subset[ , !(names(data_plot_subset) %in% "Group")] # Remove the 'Group' column
  
  # Combine TACIT_Ligands and TACIT_Receptors into one column 'Cell_Type'
  data_combined <- data_plot_subset %>%
    pivot_longer(cols = c(TACIT_Ligands, TACIT_Receptors), names_to = "Type", values_to = "Cell_Type") %>%
    group_by(Cluster, Cell_Type) %>%
    summarise(Count = n(), .groups = 'drop')
  
  # Calculate the proportion of each cell type within each cluster
  data_prop <- data_combined %>%
    group_by(Cluster) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
  
  # Remove the 'Count' column using base R
  data_prop <- data_prop[ , !(names(data_prop) %in% "Count")]
  
  # Pivot to create a matrix format for the heatmap
  data_matrix <- data_prop %>%
    pivot_wider(names_from = Cell_Type, values_from = Proportion, values_fill = list(Proportion = 0))
  
  # Convert to data frame and set row names
  data_matrix <- as.data.frame(data_matrix)
  rownames(data_matrix) <- data_matrix$Cluster
  data_matrix <- data_matrix[ , !(names(data_matrix) %in% "Cluster")] # Remove Cluster column after setting row names
  
  # Define color palette and breaks for the heatmap
  color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
  breaks <- seq(-2, 2, length.out = 101)  # Adjust to match your data range
  
  # Set the file name for the heatmap SVG output
  heatmap_filename <- paste0("C:/Users/huynhk4/Downloads/OCF_MERSCOPE_0409/Healthy/Heatmap_of_Cell_Type_Proportions_", group, ".svg")
  
  # Open SVG device
  svg(filename = heatmap_filename, width = 8, height = 6)
  
  # Plot heatmap
  pheatmap(as.matrix(data_matrix),
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           main = paste("Heatmap of Cell Type Proportions by", group),
           border_color = NA,
           scale = "row",  # Adjust scale if necessary
           color = color_palette,
           breaks = breaks)
  
  # Close SVG device
  dev.off()
}
