
# This script defines gene signatures used for the main, extend and supplementary
# figures and generates a DotPlot to visualize the expression of selected genes.


library(Seurat)
library(ggplot2)



###Gene signatures for Fig1

malignancy_signatures <- c("GATA1", "CD33", "CD7", "KIT", "ITGA2B",
                            "ITGAL", "PRAME", "DUXAP8", "GABRG2", "MAGI1",
                            "ABCB1", "GABRE", "HDC", "IRX1", "HPSE2",
                            "MYO16", "APOE", "APOC2", "LRRC7", "KYNU",
                            "COL6A5", "NTRK1", "IGF2BP1", "GP9", "CLCA1",
                            "RELN", "CD9", "LIN28B")

###Gene signatures for Fig2

fetal_to_TAM_transition <- c("GATA1", "CD33", "CD7", "KIT", "ITGA2B", "ITGAL",
                              "IGF2BP2",  "HMGA2", "HOXA9", "DPPA4", "HSPA1A",
                              "HSPA1B", "IGF2BP1", "LIN28B", "H19", "CD4",
                              "CLCA1", "RELN", "IRX1", "PAPPA","HPSE2", "MYO16",
                              "MYO18B", "APOC2", "GABRE", "LRRC7", "KYNU", "COL6A5",
                              "DUXAP8", "PRAME", "GABRG2", "MAGI1", "ABCB1")

###Gene signatures for Fig3

bridge_vs_main_blast_signatures <- c("GATA1", "CD33", "CD7", "ITGAL", "GABRE",
                                     "HPSE2", "MYO16", "DUXAP8", "GABRG2", "PRAME",
                                     "MEIS2", "RUNX1T1", "PARD3B", "ZEB2", "ETV6",
                                     "ZNF521", "AUTS2", "ZBTB16", "CBFA2T3", "MSI2",
                                     "ERG", "DYRK1A", "RUNX1", "MIR99AHG", "MAGI1",
                                     "ROBO2", "ITGA9", "EPHA3", "ABCB1", "ABCC1",
                                     "ABCC4", "APOC1", "APOC2", "GYPB", "EPCAM",
                                     "HBG2", "GP9", "ITGA2B", "CD34", "LY6E",
                                     "IFI6", "HLA-A", "HLA-B", "HLA-C", "B2M",
                                     "CALR", "HLA-DRA", "CDK4", "PCNA", "MCM5",
                                     "MCM7", "COX4I1", "UQCRC2", "ATP5F1B")


####Gene signatures for Fig4

signatures_across_diagnosis_and_relapse <- c("GATA1", "CD33", "CD7", "ITGAL", "GABRE",
                                             "HPSE2", "MYO16", "DUXAP8", "GABRG2", "PRAME",
                                             "MEIS2", "RUNX1T1", "PARD3B", "ZEB2", "ETV6",
                                             "ZNF521", "AUTS2", "ZBTB16", "CBFA2T3", "MSI2",
                                             "ERG", "DYRK1A", "RUNX1", "MIR99AHG", "MAGI1",
                                             "ROBO2", "ITGA9", "EPHA3", "ABCB1", "ABCC1",
                                             "ABCC4", "APOC1", "APOC2", "ITGA2B", "CD34",
                                             "LY6E", "IFI6", "HLA-A", "HLA-B", "HLA-C",
                                             "B2M", "CALR", "HLA-DRA", "CDK4", "PCNA",
                                             "MCM5", "MCM7", "COX4I1", "UQCRC2", "ATP5F1B")

####Gene signatures for EDFig3e

Signature_genes_GATA1_mutant_erythroid_cells <- c("IGF2BP2", "HSPA1A", "HSPA1B", "MKI67", "PCNA",
                                                  "RELN", "MYO18B", "CLCA1", "MYO16", "KYNU",
                                                  "HPSE2", "LRRC7", "IRX1", "PAPPA", "IFI27",
                                                  "HLA-A", "C5AR2", "PDZK1IP1", "GDF15", "H19",
                                                  "HBE1", "HBZ", "HBG1", "HBG2", "HBB",
                                                  "HBA1")


####Gene signatures for EDFig7c

signature_gene_expression_comparison_to_developmental_hematopoiesis <- c("GATA1", "CD33", "CD7", "ITGAL", "LRRC7",
                                                                         "HPSE2", "MYO18B", "MYO16", "EPHA3", "KYNU",
                                                                         "COL6A5", "PAPPA", "GABRE", "IRX1", "CLCA1",
                                                                         "RELN", "LIN28B", "DUXAP8", "PRAME", "GABRG2",
                                                                         "RUNX1T1", "MEIS2", "ZNF521", "ITGA9", "MSI2",
                                                                         "CBFA2T3", "ETV6", "AUTS2", "MIR99AHG", "ERG",
                                                                         "DYRK1A", "RUNX1", "PARD3B", "NLGN1", "TFPI",
                                                                         "PAM", "NUBPL", "IGF2BP3", "IGF2BP1", "IGF1R",
                                                                         "ROBO2", "MAGI1", "APOC1", "APOC2", "APOE",
                                                                         "GYPB", "GYPA", "NFE2", "KLF1", "EPCAM",
                                                                         "ITGA2B", "CD9", "GP9")

###Gene signatures for SupFig2f,g,h


MLDS_celltype_signatures <- c("HLF", "AVP", "MECOM", "HOXA9", "MLLT3",
                         "SPINK2", "CD34", "BCL11A", "IL7R", "IL2RG",
                         "FLT3", "LTB", "DNTT", "MME", "CD19",
                         "EBF1", "IGLL1", "PAX5", "VPREB1", "RAG1",
                         "IGHM", "CD79B", "IGHD", "MS4A1", "GZMB",
                         "JCHAIN", "SERPINF1", "IRF8", "CD1C", "S100A10",
                         "FCER1A", "CD3D", "CD4", "CCR7", "CD8A",
                         "CD8B", "TRAC", "KLRB1", "NKG7", "GZMK",
                         "PCNA", "TRDC", "CD69", "DUSP2", "PRSS57",
                         "CALR", "AZU1", "PRTN3", "LYZ", "KLF4",
                         "FCN1", "CD14", "CD300E", "VCAN", "HLA-DRA",
                         "CD68", "ADGRE1", "ITGAM", "S100A8", "NCF1",
                         "CLC", "IL5RA", "TPSAB1", "HDC", "HPGDS",
                         "GATA2", "KIT", "KLF1", "GATA1", "EPCAM",
                         "APOC1", "GYPA", "PRDX2", "ALAS2", "HBA1",
                         "HBZ", "HBE1", "ITGA2B", "PF4", "CMTM5",
                         "PPBP", "TREML1", "DUXAP8", "GABRE", "LRRC7",
                          "ROBO2", "PRAME", "MAGI1", "IFI27", "MYO16",
                         "RELN", "CLCA1")


###Gene signatures for SupFig5d

MLDS_integrated_fetal_tissue_signatures <- c("HLF", "AVP", "MECOM", "HOXA9", "MLLT3",
                                            "SPINK2", "CD34", "BCL11A", "IL7R", "IL2RG",
                                            "FLT3", "LTB", "DNTT", "MME", "CD19",
                                            "EBF1", "IGLL1", "PAX5", "VPREB1", "RAG1",
                                            "IGHM", "CD79B", "IGHD", "MS4A1", "GZMB",
                                            "JCHAIN", "SERPINF1", "IRF8", "CD1C", "S100A10",
                                            "FCER1A", "CD3D", "CD4", "CCR7", "CD8A",
                                            "CD8B", "TRAC", "KLRB1", "NKG7", "GZMK",
                                            "PCNA", "TRDC", "CD69", "DUSP2", "PRSS57",
                                            "CALR", "AZU1", "PRTN3", "LYZ", "KLF4",
                                            "FCN1", "CD14", "CD300E", "VCAN", "HLA-DRA",
                                            "CD68", "ADGRE1", "ITGAM", "S100A8", "NCF1",
                                            "CLC", "IL5RA", "TPSAB1", "HDC", "HPGDS",
                                            "GATA2", "KIT", "KLF1", "GATA1", "EPCAM",
                                            "APOC1", "GYPA", "PRDX2", "ALAS2", "HBA1",
                                            "HBZ", "HBE1", "ITGA2B", "PF4", "CMTM5",
                                            "PPBP", "TREML1", "CDH5", "PECAM1", "PDGFRB",
                                            "NG2", "ENG", "THY1", "HEPN1", "ARG1",
                                            "DUXAP8", "GABRE", "LRRC7", "ROBO2", "PRAME",
                                            "MAGI1", "IFI27", "MYO16", "RELN", "CLCA1")

###Gene signatures for SupFig10g

Shared_DE_genes_in_bridge_progenitors_at_relapse_vs_diagnosis <- c("ACSL4", "BRWD1", "CAMK1D", "CD36", "CDC42SE2",
                                                                   "ERG", "GNG2", "GNPTAB", "PBX4", "RALA",
                                                                   "SEPTIN2", "SIPA1L1", "ZEB1", "AC020916.1", "FTX",
                                                                   "IGF1R", "KLF13", "MBD5", "MYO16", "RABGEF1")


###generate DotPlot

p <- DotPlot(
  alldata
  ,
  features = rev(malignancy_signatures),
  cols = c("grey90", "red3")
) +
  ggplot2::coord_flip() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = "italic"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::guides(
    size  = ggplot2::guide_legend(title = NULL),
    color = ggplot2::guide_colorbar(title = NULL)
  )

ggsave(p, file= file.path("/Users/yue/malignancy_signatures.png"),
       width  = 8,
       height = 6,
       bg = "white",
       dpi    = 600)





