##Analysis of Cao data larva Cluster3-33 anterior non neural ectorderm

#Retrieved the data from geo repository with HDF5

# Loaded the dataset from the Cao paper

#If not installed
#install.packages("hdf5r")

library(rlang)
library(Seurat)
library(hdf5r)
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(sctransform)

##Load the R object if needed

larva.combined <- readRDS("C:/Users/frazy3/Documents/Seurat_scRNAseq/filtered_gene_bc_matrices/ CionaCaolarva_combined.rds")


#get the clusters 3 33 as a new subset object
##larva.combined.cluster3_33 <- subset(x = larva.combined, idents = c(3,33))

larva.combined.cluster3_33.larva1 <- subset(x = larva.combined, idents = c(3,33), subset = orig.ident == "CionaCaolarva1")
larva.combined.cluster3_33.larva2 <- subset(x = larva.combined, idents = c(3,33), subset = orig.ident == "CionaCaolarva2")
larva.combined.cluster3_33.larva3 <- subset(x = larva.combined, idents = c(3,33), subset = orig.ident == "CionaCaolarva3")





## run sctransform (instead of scaling and normalizing and find variable features)
#there is an option I do not use  to regress out variables : vars.to.regress = "percent.mt"

larva.combined.cluster3_33.larva1 <- SCTransform(larva.combined.cluster3_33.larva1, verbose = FALSE)
larva.combined.cluster3_33.larva2 <- SCTransform(larva.combined.cluster3_33.larva2, verbose = FALSE)
larva.combined.cluster3_33.larva3 <- SCTransform(larva.combined.cluster3_33.larva3, verbose = FALSE)




##Perform Integration on the subset
larva.cluster3_33.list <- list(larva.combined.cluster3_33.larva1, larva.combined.cluster3_33.larva2, larva.combined.cluster3_33.larva3)


larva.cluster3_33.features <- SelectIntegrationFeatures(object.list = larva.cluster3_33.list, nfeatures = 3000)
larva.cluster3_33.list <- PrepSCTIntegration(object.list = larva.cluster3_33.list, anchor.features = larva.cluster3_33.features)

larva.cluster3_33.anchors <- FindIntegrationAnchors(object.list = larva.cluster3_33.list, normalization.method = "SCT", anchor.features = larva.cluster3_33.features)
larva.cluster3_33.combined <- IntegrateData(anchorset = larva.cluster3_33.anchors, normalization.method = "SCT")

#PCA clustering
DefaultAssay(larva.cluster3_33.combined) <- "integrated"

larva.cluster3_33.combined <- RunPCA(larva.cluster3_33.combined, verbose = FALSE)
# UMAP and Clustering
larva.cluster3_33.combined <- RunUMAP(larva.cluster3_33.combined, reduction = "pca", dims = 1:30)
larva.cluster3_33.combined <- FindNeighbors(larva.cluster3_33.combined, reduction = "pca", dims = 1:30)
larva.cluster3_33.combined <- FindClusters(larva.cluster3_33.combined, resolution = 0.5)

plot1 <- DimPlot(larva.cluster3_33.combined, reduction = "umap", label = TRUE, label.size = 3)
plot1+ labs(title = "larva", subtitle = "Clusters 3 and 33 using 30 dims", caption ="scRNAseq from Cao et al")


plot2 <- DimPlot(larva.cluster3_33.combined, reduction = "umap", group.by="orig.ident", label = TRUE, label.size = 3)

plot2

DefaultAssay(larva.cluster3_33.combined) <- "RNA"

AllMarkerslarva.cluster3_33 <- FindAllMarkers(larva.cluster3_33.combined)
write.csv(AllMarkerslarva.cluster3_33, file = ("C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//CionaCaolarva.cluster3_33_SCtransform.csv"), sep="/t")
larva.cluster3_33.combined@misc$markers <- AllMarkerslarva.cluster3_33

# You can save the final file (with UMAP)
saveRDS(larva.cluster3_33.combined, file = "C:/Users/frazy3/Documents/Seurat_scRNAseq/filtered_gene_bc_matrices/CionaCaolarva.cluster3_33_combined.rds")

##Load the R object if needed
larva.cluster3_33.combined <- readRDS("C:/Users/frazy3/Documents/Seurat_scRNAseq/filtered_gene_bc_matrices/CionaCaolarva.cluster3_33_combined.rds")

##Load the marker list and copy it in the seurat object
AllMarkerslarva.cluster3_33 <-read.table(file="/Users/frazy3/Documents/Seurat_scRNAseq/filtered_gene_bc_matrices/CionaCaolarva.cluster3_33_SCtransform.csv", header=TRUE, row.names=1, sep = ",")
larva.cluster3_33.combined@misc$markers <- AllMarkerslarva


##Append the tables withnames and insitu data
KH_name_insitu_data <- read_excel("C:/Users/frazy3/Documents/Seurat_scRNAseq/KHID-UniqueName-URLs-InSitu-COMPLETE_KH2012.xlsx")

AllMarkerslarva.cluster3_33 <-read.table(file="/Users/frazy3/Documents/Seurat_scRNAseq/filtered_gene_bc_matrices/CionaCaolarva.cluster3_33_SCtransform.csv", header=TRUE, row.names=1, sep = ",")

AllMarkerslarva.cluster3_33_append <- left_join(AllMarkerslarva.cluster3_33,KH_name_insitu_data, by = c("gene" = "KHID")  )



write.csv(AllMarkerslarva.cluster3_33_append, file = ("C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//larva.cluster3_33_SCtransform_append.csv"), sep="/t")








#number of cells per cluster
table(Idents(larva.cluster3_33.combined))



DefaultAssay(larva.cluster3_33.combined) <- "RNA"
DefaultAssay(larva.cluster3_33.combined) <- "integrated"
DefaultAssay(larva.combined) <- "RNA"
DefaultAssay(larva.combined) <- "integrated"


##Distribution of expression of the top markers of the subclusters

larvatop10cluster33sub0 <- c("KH2012:KH.S605.3", "KH2012:KH.C5.227", "KH2012:KH.L12.5", "KH2012:KH.C7.391", "KH2012:KH.C4.615", "KH2012:KH.C1.525", "KH2012:KH.C11.667", "KH2012:KH.C7.154", "KH2012:KH.S164.24")
larvatop10cluster33sub1 <- c("KH2012:KH.L116.94", "KH2012:KH.L116.84", "KH2012:KH.L92.8", "KH2012:KH.C4.317", "KH2012:KH.L4.43", "KH2012:KH.C1.337", "KH2012:KH.C10.335", "KH2012:KH.C9.196", "KH2012:KH.C3.781")
larvatop10cluster33sub2 <- c("KH2012:KH.C1.1125", "KH2012:KH.C3.448", "KH2012:KH.C2.322", "KH2012:KH.L116.85")

larvatop10cluster33sub3 <-  c("KH2012:KH.C5.227", "KH2012:KH.C10.516", "KH2012:KH.C8.420", "KH2012:KH.C5.120", "KH2012:KH.S1555.2", "KH2012:KH.C1.486", "KH2012:KH.L141.45")

larvatop10cluster33sub4


larva.cluster3_33.combined.sub8 <- subset(x = larva.cluster3_33.combined, idents = c(8))

larva.cluster3_33.combined.sub9 <- subset(x = larva.cluster3_33.combined, idents = c(9))

larva.cluster3_33.combined.sub8_9 <- subset(x = larva.cluster3_33.combined, idents = c(8,9))


larva.cluster3_33.combined.sub7 <- subset(x = larva.cluster3_33.combined, idents = c(7))


larva.cluster3_33.combined.sub8.NoBig <- subset(x = larva.cluster3_33.combined.sub8, nFeature_RNA > 200 & nFeature_RNA < 2500)


FeaturePlot(larva.cluster3_33.combined, features = larvatop10cluster33sub3, cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 0.5)



##Finding the expression of a few classic markers
#Where is FoxG ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C8.774"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)

#Where is Islet ?
FeaturePlot(larva.cluster3_33.combinedE, features = c("KH2012:KH.L152.2"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)

#Where is Six12 ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C3.553"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 1)

#Where is EMX ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.L142.14"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)

#Where is Pou4 ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C2.42"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)


#Where is SP8 ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C13.22"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)


#Where is ZFP69B; ZNF470; ZNF713 ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.L40.18"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)


#Where are the calponins domain prot ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C5.56","KH2012:KH.C12.409", "KH2012:KH.C4.559"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)

#Where is RSPO3 ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C11.249"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)


#Where is CryBG ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.S605.3"), cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 0.5)

#Where is MRF ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C14.307"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)


#Where is COE (ASM OSM, some nervous system) ?
FeaturePlot(larva.rep_TSNE.Cluster8_10_18_TSNE, features = c("KH2012:KH.L24.10"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)

#Where is GnRH2 (ASM OSM, some nervous system) ?
FeaturePlot(larva.cluster3_33.combined, features = c("KH2012:KH.C9.484"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)


# find all markers of the clusters and list them in csv files
for (i in 0:11){
  nam <- paste0("CionaCaolarvaCluster8_10_18sub",i,".markers")
  assign(nam, FindMarkers(larva.rep_TSNE.Cluster8_10_18_TSNE, ident.1 = i, min.pct = 0.25))
  write.csv(x = get(paste0("CionaCaolarvaCluster8_10_18sub",i,".markers")) , file = paste0("C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//CionaCaolarvafilter500_Cluster8_10_18sub", i,"markers.csv"), sep="/t")
}



Cluster6markersanti10 <- FindMarkers(larva.rep_TSNE.Cluster8_10_18_TSNE, ident.1 = 6, ident.2 = 10, min.pct = 0.25)

head(Cluster6markersanti10)
write.csv(x = Cluster6markersanti10 , file = paste0("C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//CionaCaolarvafilter500_Cluster8_10_18sub6anti10markers.csv"), sep="/t")


##Select manually the small group with high expression of MYLCKC9.384 Cells using overlocator in a dimplot or a feature plot
library(ggplot2)

## To interact with the plot to check the features of the cells you are passing over
#plot <- FeaturePlot(pbmc, features = "MS4A1")
#HoverLocator(plot = plot, information = FetchData(pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))

# to select specific cellc from a scatterplot
plot <- FeaturePlot(larva.rep_TSNE.Cluster8_10_18_TSNE, features = c("KH2012:KH.C9.384"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 1)

select.cells <- CellSelector(plot = plot)
head(select.cells)

Idents(larva.rep_TSNE.Cluster8_10_18_TSNE, cells = select.cells) <- "KH.C9.384strong_Palpcluster"

# Now, we find markers that are specific to the new cells, and find clear DC markers
StrongC9.384_larva.markers <- FindMarkers(larva.rep_TSNE.Cluster8_10_18_TSNE, ident.1 = "KH.C9.384strong_Palpcluster", min.diff.pct = 0.3, 
                                              only.pos = TRUE)

head(SP8cellsCluster2_larva.markers)
write.csv(x = StrongC9.384_larva.markers , file = "C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//CionaCaolarvafilter500_Cluster8_10_18strongC9.384.csv", sep="/t")

LarvaCluster2SP8cells_top6markers <- head(row.names(x = SP8cellsCluster2_larva.markers), 6)


FeaturePlot(CionaStolfi_larva_TSNE_Cluster2_TSNE, features = c(LarvaCluster2SP8cells_top6markers), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)

FeaturePlot(CionaStolfi_larva_TSNE, features = c(LarvaCluster2SP8cells_top6markers), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 0.5)




pbmc <- CellSelector(plot = plot, object = pbmc, ident = "selected")

levels(CionaStolfi_larva)


usual <-c("KH2012:KH.C2.42", "KH2012:KH.C3.724", "KH2012:KH.S605.3", "KH2012:KH.C11.360")
                 

  
DefaultAssay(larva.cluster3_33.combined) <- "RNA"
DefaultAssay(larva.cluster3_33.combined) <- "integrated"

larva.combined.cluster3_33.sub8.combined <- subset(x = larva.cluster3_33.combined, idents = c(8))
larva.combined.cluster3_33.sub8_9.combined <- subset(x = larva.cluster3_33.combined, idents = c(8, 9))



DefaultAssay(larva.combined.cluster3_33.sub8.combined) <- "RNA"
DefaultAssay(larva.combined.cluster3_33.sub8.combined) <- "integrated"
DefaultAssay(larva.combined.cluster3_33.sub8_9.combined) <- "RNA"
DefaultAssay(larva.combined.cluster3_33.sub8_9.combined) <- "integrated"




FeaturePlot(larva.cluster3_33.combined, features = usual, cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 0.5)
FeaturePlot(larva.combined.cluster3_33.sub8.combined, features = sub6, cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 1)

FeaturePlot(larva.combined.cluster3_33.sub8.combined, features ="KH2012:KH.C11.360" , cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 1)

FeaturePlot(larva.combined.cluster3_33.sub8.combined, features =c ("KH2012:KH.S605.3", "KH2012:KH.C11.360", "KH2012:KH.C3.724", "nFeature_RNA", "nCount_RNA"), cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 1)

nFeaturePlotDetection <- FeatureScatter(object = larva.combined.cluster3_33.sub8.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nFeaturePlotDetection

sub

larva.combined.cluster3_33.sub8.combined_lowCount <- subset(x = larva.combined.cluster3_33.sub8.combined, subset = nCount_RNA < (8000))
larva.combined.cluster3_33.sub8_9.combined_lowCount <- subset(x = larva.combined.cluster3_33.sub8_9.combined, subset = nCount_RNA < (6000))


FeatureScatter(object = larva.combined.cluster3_33.sub8.combined_lowCount, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeaturePlot(larva.combined.cluster3_33.sub8.combined_lowCount, features =c ("KH2012:KH.S605.3", "KH2012:KH.C11.360", "KH2012:KH.C3.724", "KH2012:KH.S1051.1", "nFeature_RNA", "nCount_RNA"), cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 1)

FeatureScatter(object = larva.cluster3_33.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

par(mfrow = c(2, 2))
for (i in 1:4) {
  FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(i)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
}
FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(6, 7)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(6, 8)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(6, 1)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(6, 1)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="orig.ident")

group.by="orig.ident"



plot1 <-  FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(8)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <-    FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(2)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = subset(larva.cluster3_33.combined,idents = c(8)) , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


CombinePlots(plot1, plot2)




FeatureScatter(object = larva.combined.cluster3_33.sub8.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeatureScatter(object = larva.combined.cluster3_33.sub8_9.combined_lowCount, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

FeaturePlot(larva.combined.cluster3_33.sub8_9.combined_lowCount, features =c ("KH2012:KH.S605.3", "KH2012:KH.C11.360", "KH2012:KH.C3.724", "KH2012:KH.S1051.1", "nFeature_RNA", "nCount_RNA"), cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 1)

FeatureScatter(object = larva.combined.cluster3_33.sub8_9.combined, feature1 ="nCount_RNA", feature2 = "nFeature_RNA")








Dimplot1 <- DimPlot(larva.rep_TSNE.Cluster8_10_18_TSNE,, reduction = "tsne", label = TRUE, label.size = 3)
Dimplot1
HoverLocator(plot = nFeaturePlotDetection, information = FetchData(object = larva.combined.cluster3_33.sub8.combined, vars = "KH2012:KH.C3.724"))

VlnPlot(larva.cluster3_33.combined,usual)
VlnPlot(larva.cluster3_33.combined, "KH2012:KH.C4.626")

gene1 <- "KH2012:KH.L96.104" 
gene2 <-

FeatureScatter(object = larva.cluster3_33.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = larva.cluster3_33.combined, feature1 = `KH2012:KH.S605.3`, feature2 = "nFeature_RNA")


##subcluster5 are RTEN homogeneous ? it seems they are not


larva.cluster3_33.sub5.combined <- subset(x = larva.cluster3_33.combined, idents = c(5))


#PCA clustering



DefaultAssay(larva.cluster3_33.sub5.combined) <- "RNA"

larva.cluster3_33.sub5.combined <- RunPCA(larva.cluster3_33.sub5.combined, verbose = FALSE)
# UMAP and Clustering
larva.cluster3_33.sub5.combined <- RunUMAP(larva.cluster3_33.sub5.combined, reduction = "pca", dims = 1:30)
larva.cluster3_33.sub5.combined <- FindNeighbors(larva.cluster3_33.sub5.combined, reduction = "pca", dims = 1:30)
larva.cluster3_33.sub5.combined <- FindClusters(larva.cluster3_33.sub5.combined, resolution = 0.2)

plot1 <- DimPlot(larva.cluster3_33.sub5.combined, reduction = "umap", label = TRUE, label.size = 3)
plot1+ labs(title = "larva", subtitle = "Clusters 3 and 33 subcluster 5 using 30 dims", caption ="scRNAseq from Cao et al")


AllMarkerslarva.cluster3_33.sub5 <- FindAllMarkers(larva.cluster3_33.sub5.combined)



##subset 3_33 removing high RNA counts or high RNA features (both variables are correlated linearly)




#get the clusters 3 33 with cutoff for RNA counts     as a new subset object
# this was not posssible due to an issue with anchors (probably the samples are too small for this step)
# try to take a direct subset and recluster, this doesn't give better clusters (from 6, 8, 9, only two clusters are formed)
# try to manually select the couple of purest inner collocytes (with C11.360 strong with no CryBG)

##Create the subset of interest with the chosen cutoff for counts
larva.combined.cluster3_33.sub6_8_9 <- subset(x = larva.cluster3_33.combined, idents = c(6, 8, 9))
larva.combined.cluster3_33.sub6_8_9_MaxCount7000 <- subset(x = larva.combined.cluster3_33.sub6_8_9, nCount_RNA < 7000)
DimPlot(larva.combined.cluster3_33.sub6_8_9)
DimPlot(larva.combined.cluster3_33.sub6_8_9_MaxCount7000)

larva.cluster3_33.combined.sub8 <- subset(x = larva.cluster3_33.combined, idents = c(8))
larva.cluster3_33.combined.sub8_MaxCount7000 <- subset(x = larva.cluster3_33.combined.sub8, nCount_RNA < 7000)
DimPlot(larva.cluster3_33.combined.sub8)
DimPlot(larva.cluster3_33.combined.sub8_MaxCount7000)




#Select manually the small group with high expression of MYLCKC9.384 Cells using overlocator in a dimplot or a feature plot
library(ggplot2)

## To interact with the plot to check the features of the cells you are passing over
#plot <- FeaturePlot(pbmc, features = "MS4A1")
#HoverLocator(plot = plot, information = FetchData(pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))

# to select specific cellc from a scatterplot
DefaultAssay(larva.combined.cluster3_33.sub6_8_9) <- "RNA"

DefaultAssay(larva.combined.cluster3_33.sub6_8_9_MaxCount7000) <- "RNA"

plot <- FeaturePlot(larva.combined.cluster3_33.sub6_8_9_MaxCount7000, features = c("KH2012:KH.C11.360"), cols = c("lightgrey", "blue","green","yellow", "orange", "red"), pt.size = 1)
plot
select.cells <- CellSelector(plot = plot)
head(select.cells)

Idents(larva.combined.cluster3_33.sub6_8_9_MaxCount7000, cells = select.cells) <- "InnerCollocytes_larva"


# Now, we find markers that are specific to the new cells, and find clear markers
AllMarkerslarva.larva.combined.cluster3_33.sub6_8_9_MaxCount7000 <- FindAllMarkers(larva.combined.cluster3_33.sub6_8_9_MaxCount7000)

#Append annotations
KH_name_insitu_dataNEW <- read_excel("C:/Users/frazy3/Documents/Seurat_scRNAseq/ANISEEDID-KHID-UniqueName-Names-InSitus-URLS_KH2012.xlsx")


AllMarkerslarva.larva.combined.cluster3_33.sub6_8_9_MaxCount7000_append <- left_join(AllMarkerslarva.larva.combined.cluster3_33.sub6_8_9_MaxCount7000,KH_name_insitu_dataNEW, by = c("gene" = "KHID")  )
write.csv(AllMarkerslarva.larva.combined.cluster3_33.sub6_8_9_MaxCount7000_append, file = ("C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//CionaCaolarva.cluster6_8_9_InnerCollocytes_append.csv"), sep="/t")


##Append annotations new version

# Now, we find markers that are specific to the new cells, and find clear markers
DefaultAssay(larva.cluster3_33.combined) <- "RNA"
AllMarkerslarva.cluster3_33 <- FindAllMarkers(larva.cluster3_33.combined)


#Append annotations
KH_name_insitu_dataNEW <- read_excel("C:/Users/frazy3/Documents/Seurat_scRNAseq/ANISEEDID-KHID-UniqueName-Names-InSitus-URLS_KH2012.xlsx")
KH_name_TF_ZF_signaling <- read_excel("C:/Users/frazy3/Documents/Seurat_scRNAseq/TFandSignaling_ZF_KHmodels2012_2013.xlsx")

AllMarkerslarva.cluster3_33_append <- left_join(AllMarkerslarva.cluster3_33,KH_name_insitu_dataNEW, by = c("gene" = "KHID")  )
AllMarkerslarva.cluster3_33_append_TFsignaling <- left_join(AllMarkerslarva.cluster3_33_append,KH_name_TF_ZF_signaling, by = c("gene" = "KH2012model")  )

write.csv(AllMarkerslarva.cluster3_33_append_TFsignaling, file = ("C://Users//frazy3//Documents//Seurat_scRNAseq//filtered_gene_bc_matrices//CionaCaolarva.cluster3_33_append_TFZFsignaling.csv"), sep="/t")






















## the issue if you go only for theInnercolocyte markers here is that it doesn't give a gene colum
InnerCollocytes_larva.markers <- FindMarkers(larva.combined.cluster3_33.sub6_8_9_MaxCount7000, ident.1 = "InnerCollocytes_larva", min.diff.pct = 0.3, 
                                          only.pos = TRUE)

#visualization of the top Innercollocyte markers
InnercollocytesLarva1_9 <-c("KH2012:KH.C11.360", "KH2012:KH.C2.607", "KH2012:KH.C9.423", "KH2012:KH.C11.404", "KH2012:KH.C2.120", "KH2012:KH.S1090.3", "KH2012:KH.C9.780", "KH2012:KH.C11.249", "KH2012:KH.C3.665")
FeaturePlot(larva.combined.cluster3_33.sub6_8_9_MaxCount7000, features =c (InnercollocytesLarva1_9), cols = c("lightgrey", "blue","green", "orange", "red"), pt.size = 1, order = TRUE)




