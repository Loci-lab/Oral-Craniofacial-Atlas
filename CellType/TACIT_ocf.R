#-----------------------------------------------------------------------------
###Goal: Obtain TACIT: Threshold-based Assignment for Cell Type Identifying in Spatial Omics
#-----------------------------------------------------------------------------
###Cell type annotation Step: 
#### 1) Load the data
#### 2) Load signature matrix
#### 3) Quality Control
#### 4) Run TACIT
#-----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(segmented)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

###Colon datasets (Figure 2 top)

#-----------------------------------------------------------------------------


#### 1) Load the data
data <-  readr::read_csv("OCF_expression.csv")
#### 2) Load signature matrix
Signature <- readr::read_csv("Signature_OCF.csv")

#### 3) Run TACIT
CRC_TACIT=TACIT_update(data_anb = (data[,colnames(Signature)[-c(1)]]),r = 200,p=20,Signature = Signature)

#Get annotation
data_predicted=CRC_TACIT[3]
#Get annotation with cell ID
data_predicted=data.frame(cell_ID=data$CellID,TACIT=data_predicted$mem)
#Validate TACIT
confusionMatrix(data=factor(data_predicted$TACIT,levels = c(unique(data$ClusterName))),
                reference = factor(data$ClusterName,levels = c(unique(data$ClusterName))))

write.csv(data_predicted,"result_TACIT.csv")

