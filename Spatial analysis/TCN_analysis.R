
library(ggplot2)
library(readr)
data=data.table::fread("C:/Users/huynhk4/Downloads/combined_spatialLDA_PCF_v3 (1).csv")
data=as.data.frame(data)

data_final=data
user_colors <- c(
  "Acini" = "cyan",                    # Acinar Cells
  "B cells" = "#FFA500",               # B Cells
  "DC cells" = "#FFD700",              # Dendritic Cells
  "Ducts" = "#00FF00",                 # Ductal Cells
  "Epithelial" = "powderblue",         # Intermediate Epithelium
  "Fibroblast" = "#f032e6",            # Fibroblasts
  "LECs" = "orange",                   # LECs
  "Macrophage" = "gold3",              # M2 Macrophages
  "Melanocyte" = "purple",             # IgG Plasma Cells
  "Myfibroblast" = "#f032e6",          # Fibroblast Progen
  "Myoepithelial" = "blue",            # Myoepithelial Cells
  "neural/cres" = "brown3",            # Memory T cells
  "Neutrophil" = "yellow3",            # Neutrophils
  "NK cells" = "darkred",              # NK Cells
  "Tc" = "#CC0000",                    # CD8+ Effector T Cells
  "Th" = "#FF0000",                    # CD4+ T Cells / T Helper
  "Treg" = "brown2",                   # Regulatory T Cells
  "VEC" = "orange3",                   # VEC
  "VEC Progen" = "orange1",            # VEC Progen
  "Others" = "grey90"                  # Others
)




data_final$Group=ifelse(data_final$Group%in%c("Buccal_Mucosa_1_5_8_9","Gingiva_1","Gingiva_2","Tongue_10",
                                                "Tongue_12","Tongue_62"),"mucosal","glands")

data_final$TACIT=ifelse(data_final$TACIT=="neural/cres","Others",data_final$TACIT)

#colnames(data_mucosal)


table(data_final$TACIT)











data_final$CT_Group=paste0(data_final$TACIT,"::",data_final$Group,sep="")


# Define immune-related cell types
immune_cell_types <- c("B cells", "DC cells", "Macrophage", "NK cells", "Tc", "Th", "Treg", "Neutrophil")

# Filter data_final for immune cell types
data_final <- data_final[data_final$TACIT %in% immune_cell_types, ]


library(ggplot2)
library(pheatmap)
library(dplyr)




library(dplyr)
# Group the data by the predicted label and calculate the mean for each column
data_plot=data.frame(TACIT=data_final$CT_Group,data_final[,names(thresholds)])
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



# Reorder the row names of the marker_proportion dataframe
desired_order <- c(
  "B cells::mucosal", "NK cells::mucosal", "Macrophage::mucosal", "Neutrophil::mucosal",
  "DC cells::mucosal", "Th::mucosal", "Tc::mucosal", "Treg::mucosal",
  "B cells::glands", "NK cells::glands", "Macrophage::glands", "Neutrophil::glands",
  "DC cells::glands", "Th::glands", "Tc::glands", "Treg::glands"
)

# Reorder the marker_proportion dataframe by the desired order
marker_proportion_ordered <- aa[desired_order, ]

# Check the result
marker_proportion_ordered
pheatmap(as.matrix((marker_proportion_ordered)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color =color_palette,
         main = "",
         fontsize_row = 15,
         fontsize_col = 20,show_rownames = T,scale = "row",breaks = breaks)


pheatmap(as.matrix((marker_proportion_ordered)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color =color_palette,
         main = "",
         fontsize_row = 15,
         fontsize_col = 20,show_rownames = T,scale = "column",breaks = breaks)












# Define the thresholds
thresholds <- c(IDO1 = 50, HLA_A = 65, CD45RO=40,HLA_DR = 110, PD_L1 = 20, PD_1 = 15, IFNG = 20, Bcl_2 = 110, ICOS = 50)

thresholds <- c(IDO1 = 50, HLA_A = 65, CD45RO=40, PD_L1 = 20, PD_1 = 15, IFNG = 20, Bcl_2 = 110, ICOS = 50)

# Check for duplicate column names
colnames(data_final)

# Rename the duplicate columns to ensure they are unique
colnames(data_final)[1] <- "V1_1"
colnames(data_final)[2] <- "V1_2"




# Create a binary matrix: 1 if the marker is greater than the threshold, 0 otherwise
data_mucosal_binary <- data_final %>%
  mutate(IDO1 = ifelse(IDO1 > thresholds["IDO1"], 1, 0),
         HLA_A = ifelse(HLA_A > thresholds["HLA_A"], 1, 0),
#         HLA_DR = ifelse(HLA_DR > thresholds["HLA_DR"], 1, 0),
         PD_L1 = ifelse(PD_L1 > thresholds["PD_L1"], 1, 0),
         PD_1 = ifelse(PD_1 > thresholds["PD_1"], 1, 0),
         IFNG = ifelse(IFNG > thresholds["IFNG"], 1, 0),
         CD45RO = ifelse(CD45RO > thresholds["CD45RO"], 1, 0),
         ICOS = ifelse(ICOS > thresholds["ICOS"], 1, 0))

# Group by CT_Group and calculate the proportion of cells where each marker > threshold
marker_proportion <- data_mucosal_binary %>%
  group_by(CT_Group) %>%
  summarise(across(c(IDO1, HLA_A, PD_L1, PD_1, IFNG, ICOS,CD45RO), mean))

marker_proportion=as.data.frame(marker_proportion)
rownames(marker_proportion)=marker_proportion$CT_Group

# Reorder CT_Group: Buccal Mucosa, Tongue, Gingiva, Immune, Structure
# order_CT_Group <- c(
#   "Buccal Mucosa::B cells", "Buccal Mucosa::DC cells", "Buccal Mucosa::Macrophage", 
#   "Buccal Mucosa::Neutrophil", "Buccal Mucosa::NK cells", "Buccal Mucosa::Tc", 
#   "Buccal Mucosa::Th", "Buccal Mucosa::Treg", 
#   "Gingiva::B cells", "Gingiva::DC cells", "Gingiva::Macrophage", "Gingiva::Neutrophil", 
#   "Gingiva::NK cells", "Gingiva::Tc", "Gingiva::Th", "Gingiva::Treg", 
#   "Tongue::B cells", "Tongue::DC cells", "Tongue::Macrophage", "Tongue::Neutrophil", 
#   "Tongue::NK cells", "Tongue::Tc", "Tongue::Th", "Tongue::Treg", 
#   "Buccal Mucosa::Epithelial", "Buccal Mucosa::Fibroblast", "Buccal Mucosa::LECs", 
#   "Buccal Mucosa::VEC", "Buccal Mucosa::Others", 
#   "Gingiva::Epithelial", "Gingiva::Fibroblast", "Gingiva::LECs", "Gingiva::VEC", "Gingiva::Others", 
#   "Tongue::Epithelial", "Tongue::Fibroblast", "Tongue::LECs", "Tongue::VEC", "Tongue::Others"
# )

# Extract the rows in the desired order
# marker_proportion$CT_Group <- factor(marker_proportion$CT_Group, levels = order_CT_Group)
# marker_proportion <- marker_proportion %>%
#   arrange(CT_Group)

# Set the row names to CT_Group and remove the CT_Group column for the heatmap
rownames(marker_proportion) <- marker_proportion$CT_Group
marker_proportion <- marker_proportion %>% select(-CT_Group)


# Custom color palette from blue to white to red
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the color scale
breaks <- seq(-2, 2, length.out=100)  # Adjust to match your data range

# Plot the heatmap
pheatmap(as.matrix((marker_proportion)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color =color_palette,
         main = "",
         fontsize_row = 15,
         fontsize_col = 20,show_rownames = T,scale = "column",breaks = breaks)




# Reorder the row names of the marker_proportion dataframe
desired_order <- c(
  "B cells::mucosal", "NK cells::mucosal", "Macrophage::mucosal", "Neutrophil::mucosal",
  "DC cells::mucosal", "Th::mucosal", "Tc::mucosal", "Treg::mucosal",
  "B cells::glands", "NK cells::glands", "Macrophage::glands", "Neutrophil::glands",
  "DC cells::glands", "Th::glands", "Tc::glands", "Treg::glands"
)

# Reorder the marker_proportion dataframe by the desired order
marker_proportion_ordered <- marker_proportion[desired_order, ]
breaks <- seq(-2, 2, length.out=100) 
# Check the result
marker_proportion_ordered
pheatmap(as.matrix((marker_proportion_ordered)), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color =color_palette,
         main = "",
         fontsize_row = 15,
         fontsize_col = 20,show_rownames = T,scale = "column",breaks = breaks)




test=data_mucosal_binary[,c(names(thresholds)[-7],"CT_Group")]

# List of unique cell types for mucosa and glands
cell_types <- c("B cells", "DC cells", "Macrophage", "Neutrophil", "NK cells", "Tc", "Th", "Treg")

# Initialize a result data frame to store the mean proportions and p-values
results <- data.frame(
  Cell_Type = character(),
  Marker = character(),
  Mean_Mucosa = numeric(),
  Mean_Glands = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Function to run the chi-square test for each marker and cell type
run_test <- function(cell_type, marker) {
  # Subset data for the specific cell type and marker
  mucosa_data <- test[test$CT_Group == paste0(cell_type, "::mucosal"), marker]
  glands_data <- test[test$CT_Group == paste0(cell_type, "::glands"), marker]
  
  # Create a contingency table
  contingency_table <- table(c(mucosa_data, glands_data), 
                             rep(c("mucosa", "glands"), times = c(length(mucosa_data), length(glands_data))))
  
  # Perform chi-square test (you can also use fisher.test for smaller sample sizes)
  test_result <- chisq.test(contingency_table)
  
  # Calculate mean proportion for each group
  mean_mucosa <- mean(mucosa_data)
  mean_glands <- mean(glands_data)
  
  # Store the results in the results data frame
  results <<- rbind(results, data.frame(
    Cell_Type = cell_type,
    Marker = marker,
    Mean_Mucosa = mean_mucosa,
    Mean_Glands = mean_glands,
    P_Value = test_result$p.value
  ))
}

# List of markers
markers <- colnames(test)[1:7]  # First 7 columns are markers (IDO1, HLA_A, CD45RO, etc.)

# Loop through each cell type and marker, and run the test
for (cell_type in cell_types) {
  for (marker in markers) {
    run_test(cell_type, marker)
  }
}

# View the results
print(results)
