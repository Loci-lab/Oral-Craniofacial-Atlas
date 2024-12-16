library(readr)
LR_CT_ST_IM_1_ <- read_csv("C:/Users/huynhk4/Downloads/LR_CT_ST_IM (1).csv")
library(dplyr)
library(ggplot2)

# 1. Define groups for Var2: mucosa, glands, and disease (perio)
LR_CT_ST_IM_1_ <- LR_CT_ST_IM_1_ %>%
  mutate(Group = case_when(
    Var2 %in% c("BuccalMucosa", "Gingiva1", "Gingiva2", "Tongue12", "Tongue62") ~ "Mucosa",
    Var2 %in% c("MSG1", "MSG8", "MSG9", "Parotid20", "Parotid35", "Parotid39", "Parotid40", "Submandibular53_56") ~ "Glands",
    Var2 %in% c("perio_11", "perio_12", "perio_13") ~ "Disease"
  ))

# 2. Define interactions as homotypic or heterotypic
LR_CT_ST_IM_1_ <- LR_CT_ST_IM_1_ %>%
  mutate(Interaction_Type = case_when(
    Var1 %in% c("Immune to Immune", "Structure to Structure") ~ "Homotypic",
    Var1 %in% c("Immune to Structure", "Structure to Immune") ~ "Heterotypic"
  ))

# 3. Create a boxplot to compare the frequency of homotypic vs. heterotypic interactions by group (mucosa, glands, disease)
ggplot(LR_CT_ST_IM_1_, aes(x = Interaction_Type, y = Freq, fill = Group)) +
  geom_boxplot() +
  labs(title = "Comparison of Homotypic and Heterotypic Interactions Across Niches",
       x = "Interaction Type",
       y = "Frequency",
       fill = "Group") +
  theme_minimal()



# 1. Calculate proportion of Freq within each Group and Interaction_Type
LR_CT_ST_IM_1_ <- LR_CT_ST_IM_1_ %>%
  group_by(Group) %>%
  mutate(Prop_Freq = Freq / sum(Freq)) %>%
  ungroup()

# 2. Create a boxplot of the proportions
ggplot(LR_CT_ST_IM_1_, aes(x = Interaction_Type, y = Prop_Freq, fill = Group)) +
  geom_boxplot() +
  labs(title = "Comparison of Homotypic and Heterotypic Interaction Proportions Across Niches",
       x = "Interaction Type",
       y = "Proportion of Frequency",
       fill = "Group") +
  scale_y_continuous(labels = scales::percent_format()) +  # To show proportions as percentages
  theme_minimal()





library(readr)
data_final <- read_csv("C:/Users/huynhk4/Downloads/MERSCOPE_LR_2309 (3)/MERSCOPE_LR_2309 (3)/MERSCOPE_LR_2309/Mucosal vs Glands/data_final.csv")

data_final$TACIT_Ligands=data1$TACIT[data_final$cell_id]
data_final$TACIT_Receptors=data1$TACIT[data_final$Neighbor_ID]


head(data_final$TACIT_Ligands)
head(data_final$TACIT_Receptors)



# Define immune and structure cell types
immune_cells <- c("B Cells", "CD4 T Cells", "CD8 T Cells", "Dendritic Cells", "gd T Cells", 
                  "Langerhans Cells", "Mast Cells", "Monocyte-Macrophage", "NK Cells", "Plasma Cells")

structure_cells <- c("Acinar Cells", "Basal Keratincytes", "Ductal Epithelial Cells", "Fibroblasts", 
                     "Glial/Neuron", "Ionocytes", "LECs", "Mural Cells", "Myoepithelial Cells", 
                     "Skeletal Myocytes", "Suprabasal Keratinocytes", "VECs")

# Function to classify interactions
interaction_type <- function(ligand, receptor) {
  if (ligand %in% immune_cells && receptor %in% immune_cells) {
    return("Immune to Immune")
  } else if (ligand %in% structure_cells && receptor %in% structure_cells) {
    return("Structure to Structure")
  } else if (ligand %in% immune_cells && receptor %in% structure_cells) {
    return("Immune to Structure")
  } else if (ligand %in% structure_cells && receptor %in% immune_cells) {
    return("Structure to Immune")
  } else {
    return(NA)
  }
}

# Apply the function to create the interaction vector
data_final$Interaction_Type <- mapply(interaction_type, data_final$TACIT_Ligands, data_final$TACIT_Receptors)

# Create a Homotypic or Heterotypic vector
data_final$Homotypic_Heterotypic <- ifelse(data_final$Interaction_Type %in% c("Immune to Immune", "Structure to Structure"),
                                           "Homotypic", "Heterotypic")

# View the result
head(data_final[, c("TACIT_Ligands", "TACIT_Receptors", "Interaction_Type", "Homotypic_Heterotypic")])



# Recreate proportional_data with Group_Category included
proportional_data <- data_final %>%
  mutate(Group_Category = case_when(
    Group %in% mucosa_groups ~ "Mucosa",
    Group %in% glands_groups ~ "Glands"
  )) %>%
  group_by(Group, Homotypic_Heterotypic, Group_Category) %>%
  summarise(Count = n()) %>%
  group_by(Group) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Now plot the boxplot for proportions between Homotypic and Heterotypic interactions by Group Category
ggplot(proportional_data, aes(x = Homotypic_Heterotypic, y = Proportion, fill = Group_Category)) +
  geom_boxplot() +
  labs(title = "Proportion of Homotypic and Heterotypic Interactions by Group",
       x = "Interaction Type",
       y = "Proportion of Frequency",
       fill = "Group Category") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show proportions as percentages
  theme_minimal()



# Install and load ggpubr if not already installed
# install.packages("ggpubr")
library(ggpubr)

# Create the boxplot with statistical comparison
ggplot(proportional_data, aes(x = Homotypic_Heterotypic, y = Proportion, fill = Group_Category)) +
  geom_boxplot() +
  labs(title = "Proportion of Homotypic and Heterotypic Interactions by Group",
       x = "Interaction Type",
       y = "Proportion of Frequency",
       fill = "Group Category") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show proportions as percentages
  theme_classic(base_size = 15) +
  stat_compare_means(aes(group = Group_Category), method = "wilcox.test", label = "p.signif")  # Add p-values







library(dplyr)
library(ggplot2)

# Define immune and structure cell types
immune_cells <- c("B Cells", "CD4 T Cells", "CD8 T Cells", "Dendritic Cells", "gd T Cells", 
                  "Langerhans Cells", "Mast Cells", "Monocyte-Macrophage", "NK Cells", "Plasma Cells")

structure_cells <- c("Acinar Cells", "Basal Keratincytes", "Ductal Epithelial Cells", "Fibroblasts", 
                     "Glial/Neuron", "Ionocytes", "LECs", "Mural Cells", "Myoepithelial Cells", 
                     "Skeletal Myocytes", "Suprabasal Keratinocytes", "VECs")

# Function to classify interactions
interaction_type <- function(ligand, receptor) {
  if (ligand %in% immune_cells && receptor %in% immune_cells) {
    return("Immune to Immune")
  } else if (ligand %in% structure_cells && receptor %in% structure_cells) {
    return("Structure to Structure")
  } else if (ligand %in% immune_cells && receptor %in% structure_cells) {
    return("Immune to Structure")
  } else if (ligand %in% structure_cells && receptor %in% immune_cells) {
    return("Structure to Immune")
  } else {
    return(NA)
  }
}

# Apply the interaction_type function to the dataset
data_final$Interaction_Type <- mapply(interaction_type, data_final$TACIT_Ligands, data_final$TACIT_Receptors)

# Define Mucosa and Glands groups
mucosa_groups <- c("BuccalMucosa", "BuccalMucosa_v1", "BuccalMucosa_v2", "BuccalMucosa_v3", 
                   "Gingiva1", "Gingiva2", "Tongue12", "Tongue62")
glands_groups <- c("MSG_v1", "MSG9", "Parotid20", "Parotid35", "Parotid39", "Parotid40", "Submandibular53_56")

# Add a new column to classify Groups as Mucosa or Glands
data_final <- data_final %>%
  mutate(Group_Category = case_when(
    Group %in% mucosa_groups ~ "Mucosa",
    Group %in% glands_groups ~ "Glands"
  ))

# Calculate proportion of interaction types within each group
proportional_data <- data_final %>%
  group_by(Group, Interaction_Type, Group_Category) %>%
  summarise(Count = n()) %>%
  group_by(Group) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportional_data=proportional_data[which(is.na(proportional_data$Interaction_Type)==F),]
# Plot the boxplot for proportions of interaction types by Group Category (Mucosa, Glands)
ggplot(proportional_data, aes(x = Interaction_Type, y = Proportion, fill = Group_Category)) +
  geom_boxplot() +
  labs(title = "Proportion of Interaction Types by Group (Mucosa vs Glands)",
       x = "",
       y = "Proportion of Frequency",
       fill = "Group Category") +
  scale_y_continuous(labels = scales::percent_format()) +  # Show proportions as percentages
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels by 45 degrees
  stat_compare_means(aes(group = Group_Category), method = "wilcox.test", label = "p.signif")  # Add p-values



# Consolidate the groups into broader categories
proportional_data <- proportional_data %>%
  mutate(Consolidated_Group = case_when(
    grepl("BuccalMucosa", Group) ~ "BuccalMucosa",
    grepl("Tongue", Group) ~ "Tongue",
    grepl("Gingiva", Group) ~ "Gingiva",
    grepl("Parotid", Group) ~ "Parotid",
    grepl("MSG", Group) ~ "MSG",
    grepl("Submandibular", Group) ~ "Submandibular",
    TRUE ~ "Other"
  ))
# Plot the boxplot with consolidated groups
ggplot(proportional_data, aes(x = Consolidated_Group, y = Proportion, fill = Interaction_Type)) +
  geom_boxplot() +
  labs(title = "Proportion of Interaction Types by Consolidated Group (e.g., Tongue, BuccalMucosa)",
       x = "Consolidated Group",
       y = "Proportion of Frequency",
       fill = "Interaction Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

proportional_data$Group_HH=ifelse(proportional_data$Interaction_Type%in%c("Immune to Immune","Structure to Structure"),"Homotypic","Heterotypic")

ggplot(proportional_data, aes(fill = Consolidated_Group, y = Proportion, x = Group_HH)) +
  geom_boxplot() +
  labs(title = "Proportion of Interaction Types by Niches Group",
       x = "",
       y = "Proportion of Frequency",
       fill = "Interaction Type") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
