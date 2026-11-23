
#This script covers Seurat object creation, sample merging, quality control, doublet detection, normalization, dimensionality reduction, clustering, and UMAP visualization.

#DoubletFinder is applied independently to each sample, with the expected multiplet rate estimated based on the number of recovered cells. Homotypic doublet proportions are taken into account when estimating the expected number of doublets.

#Cells are filtered based on the number of detected features and mitochondrial gene expression before downstream analysis.


library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyverse)
library(DoubletFinder)
library(qs)

#######creat seurat object

data_dir <- "/Users/yue/MLDS272d0_gene_bc_matrices"
list.files(data_dir)
data <- Read10X(data.dir = data_dir)
CB597  = CreateSeuratObject(counts = data, project = "CB597", min.cells = 3, min.features = 200)

#Merge object and join layer

alldata <-merge(x = MLDS272Pd0 , y=c(MLDS272PTAM))
alldata <- JoinLayers(alldata)



#######QC

#doublets remove

run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[['SampleID']]), '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.002,0.004, 0.008, 0.0120, 0.0160, 0.020, 0.024, 0.028, 0.032, 0.036, 0.040, 0.044, 0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080),
                                      'Loaded_cells' = c(725, 1450, 2900, 4350, 5800, 7250, 8700, 10150, 11600, 13050, 14500,15950,17400,18850,20300,21750,23200,24650,26100,27550,29000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,11000,12000,13000,14000,15000,16000,17000,18000,19000,20000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  ## Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # run DoubletFinder
  sample <- doubletFinder(seu = sample, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  # change name of metadata column with Singlet/Doublet information
  colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
  
  # Subset and save
  # head(sample@meta.data['doublet_finder'])
  # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
  # singlets$ident
  double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
  double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
  return(double_finder_res)
}

#load data
seu <- alldata
Idents(seu) <- "orig.ident"
seu$SampleID <- Idents(seu)
table(seu$orig.ident)

## DoubletFinder -------------------------------
# DoubletFinder should be run on a per sample basis
samp_split <- SplitObject(seu, split.by = "SampleID") 
# Get Doublet/Singlet IDs by DoubletFinder()
samp_split <- lapply(samp_split, run_doubletfinder_custom) # get singlet/doublet assigned to each of the cell IDs (each element of the list is a different sample)

# You can also set the multiplet rate manually:
# samp_split <- # samp_split <-lapply(samp_split, run_doubletfinder_custom, multiplet_rate = 0.01)

#Add DoubletFinder results to Seurat object
sglt_dblt_metadata <- data.frame(bind_rows(samp_split)) # merge to a single dataframe
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names # assign cell IDs to row names to ensure match
sglt_dblt_metadata$row_names <- NULL
head(sglt_dblt_metadata)
seu <- AddMetaData(seu, sglt_dblt_metadata, col.name = 'doublet_finder')

#Summary of doublet detection results
# Check how doublets singlets differ in QC measures per sample.
VlnPlot(seu, group.by = 'SampleID', split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo", "percent_hb"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')

#remove doublets
Idents(alldata) <- alldata$doublet_finder
alldata <- subset(alldata, idents = c("Singlet"))

#####filter low-quality cells based on gene counts and mitochondrial content
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")
VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE, pt.size = 0)

alldata<- subset(alldata, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 15)

alldata <- alldata %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.)) 

ElbowPlot(alldata)

alldata <- alldata %>%
  FindNeighbors(dims = 1:15, reduction = "pca") %>%
  FindClusters(resolution = 0.5, cluster.name = "seurat_clusters_0.5") %>%
  RunUMAP(dims = 1:15, reduction = "pca",
          reduction.name = "umap_15", return.model = TRUE)


qsave(alldata, "/Users/yue/MLDS.qs")


