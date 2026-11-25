
# This script integrates mutation information into the Seurat metadata,
# processes variant read-support ratios, and generates UMAP plots
# showing the distribution of mutant cells across samples.

library(qs)
library(dplyr)
library(Seurat)
library(ggplot2)

#load mutation information

mutations <- read.csv("/Users/yue/variants.table.out", sep = '\t')

names(mutations)[1] <- 'sm'
names(mutations)[2] <- 'bc'
mutations$bc <- gsub("-1$", "", mutations$bc)
mutations$bc <- paste(mutations$sm, mutations$bc, sep = "_")
head(mutations)

#load seurat object
alldata <- qread("/Users/yue/MLDS.qs")

# Extract current Seurat metadata
alldata[['bc']] <- colnames(alldata)
alldata_metadata <- alldata@meta.data
alldata_metadata$bc <- rownames(alldata_metadata)

# Merge metadata (left join to keep all Seurat barcodes)
merged_metadata <- left_join(alldata_metadata, mutations, by = "bc")

# Restore row names (barcodes) and remove redundant column
rownames(merged_metadata) <- merged_metadata$bc
merged_metadata$bc <- NULL
print("B")

# Add updated metadata back to Seurat object
alldata@meta.data <- merged_metadata
head(alldata@meta.data)

#mutation information precessing
sr_cols <- grep("sr_ratio", colnames(alldata@meta.data), value = TRUE)

for (col in sr_cols) {
  alldata[[col]] <- factor(
    ifelse(is.na(alldata[[col]]), "NA",
           ifelse(alldata[[col]] == 0, "0", ">0")),
    levels = c("0", ">0", "NA")
  )
}


#illustration
colors <- c(">0" = "#4DBBD5", "0" = "#9a9a9a", "NA" = "#d9d9d9")

mk_dimplot_split <- function(group.by) {
  
  qplist <- DimPlot(
    alldata,
    reduction = "umap_15",
    group.by  = group.by,
    pt.size   = 0.3,
    combine   = FALSE,
    cols      = colors,
    order = c(">0","0","NA"), raster = F, split.by = "sample_order"
  )
  qplist <- lapply(
    qplist,
    function(p) {
      p +
        theme_void() +
        theme(
          plot.title      = element_blank(),   
          legend.position = "none",           
          strip.text      = element_blank(),   
          strip.background= element_blank(),  
          plot.margin     = unit(c(0, 0, 0, 0), "cm")
        )
    }
  )
  
  cowplot::plot_grid(
    plotlist = qplist,
    nrow = 1
  )
}

outdir <- "/Users/yue/fig1/output"

cowplot::ggsave2(
  mk_dimplot_split("MLDS272Pd0_GATA1"),
  file   = file.path(outdir, "P272Pd0_GATA1.png"),
  width  = 2.4 * 3,
  height = 2.25,
  dpi    = 600
)



