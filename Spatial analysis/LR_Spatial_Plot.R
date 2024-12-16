library(readr)
library(readxl)
library(ggplot2)
user_colors <- c(
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

data=read_csv("C:/Users/huynhk4/Downloads/MERSCOPE_LR_2309 (3)/MERSCOPE_LR_2309/Mucosal vs Glands/data_final.csv")

table(data$Group)
#choose sample here
data_sub=data[which(data$Group=="BuccalMucosa_v1"),]
ggplot() +
  geom_segment(data = data_sub, 
               aes(x = X_cell, y = Y_cell, xend = X_neighbor, yend = Y_neighbor, 
                   color = as.factor(Cluster)),  # Correct reference to Cluster
               na.rm = TRUE, size = 0.1) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = user_colors) +  # Ensure user_colors matches the factors in Cluster
  labs(x = "X", y = "Y", title = "") +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Cluster"))

#proportion of cluster LR for each sample


round(table(data_sub$Cluster)/nrow(data_sub),2)










