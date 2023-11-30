#!/usr/bin/env Rscript

# This script performs HLCA v2 annotation

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Usage: Rscript sc_annotation_HLCAv2.R <object> <outDIR name>", call.=FALSE)
}

# load packages ----------------------------------------------------------------
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pals))
suppressPackageStartupMessages(library(glue))

# parse arguments  -------------------------------------------------------------
obj <- args[1]
dirName <- args[2]
dirSCE <- glue(dirName, '/data/SCEs/')
dirPLOT <- glue(dirName, '/data/figs/')
capture <- unlist(strsplit(basename(obj),'.',fixed=TRUE))[1]

if(!dir.exists(dirName)) {
  stop("Error: out directory not exist! Check if clustering has been done.")
}

# Read object ------------------------------------------------------------------
seu <- readRDS(obj)

# Slim down object -------------------------------------------------------------
DefaultAssay(seu) <- "RNA"
seu.diet <- DietSeurat(seu, assays = "RNA")

# run Azimuth ------------------------------------------------------------------
seu.diet <- RunAzimuth(seu.diet,
                       "/group/canc2/anson/reference/annotation/HLCA_v2",
                       umap.name = 'azimuth.umap',
                       verbose = TRUE,
                       assay = "RNA",
                       k.weight = 50,
                       n.trees = 20,
                       mapping.score.k = 100)
# rds format from the web Azimuth app:
# list(umap=umap, predicted.df=df)
umap <- seu.diet@reductions$azimuth.umap
df <- data.frame(row.names=rownames(seu.diet@meta.data),
                 cell = rownames(seu.diet@meta.data),
                 predicted.ann_level_1=seu.diet$predicted.ann_level_1,
                 predicted.ann_level_2=seu.diet$predicted.ann_level_2,
                 predicted.ann_level_3=seu.diet$predicted.ann_level_3,
                 predicted.ann_level_4=seu.diet$predicted.ann_level_4,
                 predicted.ann_level_5=seu.diet$predicted.ann_level_5,
                 predicted.ann_finest_level=seu.diet$predicted.ann_finest_level,
                 predicted.ann_level_1.score=seu.diet$predicted.ann_level_1.score,
                 predicted.ann_level_2.score=seu.diet$predicted.ann_level_2.score,
                 predicted.ann_level_3.score=seu.diet$predicted.ann_level_3.score,
                 predicted.ann_level_4.score=seu.diet$predicted.ann_level_4.score,
                 predicted.ann_finest_level.score=seu.diet$predicted.ann_finest_level.score,
                 mapping.score=seu.diet$mapping.score)

# generate the azimuth result (which is a list) --------------------------------
azimuth.result <- list(umap = umap, pred.df = df)

# save the list ----------------------------------------------------------------
saveRDS(azimuth.result, file = glue(dirSCE, 'azimuth.result.rds'))

# add annotation ---------------------------------------------------------------
seu <- AddAzimuthResults(seu, filename=glue(dirSCE, 'azimuth.result.rds'))
seu$predicted.ann_level_1 <- fct_drop(seu$predicted.ann_level_1)
seu$predicted.ann_level_2 <- fct_drop(seu$predicted.ann_level_2)
seu$predicted.ann_level_3 <- fct_drop(seu$predicted.ann_level_3)
seu$predicted.ann_level_4 <- fct_drop(seu$predicted.ann_level_4)
seu$predicted.ann_finest_level <- fct_drop(seu$predicted.ann_finest_level)

out <- glue(dirSCE, capture, '.HLCA.SEU.rds')

# save object ------------------------------------------------------------------
if(!file.exists(out)) {
  saveRDS(seu, out)
}

# Split and save object --------------------------------------------------------
labels <- levels(factor(seu$predicted.ann_level_3))
macrophages <- c("Macrophages","Monocytes")
tcells <- c("T cell lineage", "Innate lymphoid cell NK")
#myeloid <- c("Dendritic cells", "Mast cells", "Monocytes", "Macrophages")
#epithelial <- c("Basal", "Multiciliated lineage", "Rare", "Secretory")

# lung <- c("AT1", "Rare", "Secretory", "Basal", "EC capillary", #old
#           "EC venous", "Multiciliated lineage")

subList <- list(macrophages = macrophages,
                tcells = tcells,
                others = labels[!labels %in% c(macrophages, tcells)]
                )

Idents(seu) <- "predicted.ann_level_3"

for(sub in names(subList)){
  out <- here(glue(dirSCE, capture, ".HLCA.{sub}.SEU.rds"))
  message(sub)
  a <- subset(seu, idents = subList[[sub]])
  DefaultAssay(a) <- "RNA"
  if(!file.exists(out)){
    saveRDS(DietSeurat(a, assays = "RNA", dimreducs = NULL, graphs = NULL), out)
  }
}


# Visualize UMAP
colours <- c(brewer.paired(12), brewer.dark2(8), brewer.set2(8), rev(brewer.accent(8)))
# Visualize annotation
p1 <- DimPlot(seu, reduction = "umap",
              cols = colours,
              group.by = "predicted.ann_level_3", label = FALSE,
              label.size = 2.5, repel = TRUE, label.box = FALSE) +
  #NoLegend() +
  ggtitle("predicted.ann_level_3")  +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        title = element_text(size = 9))
p2 <- DimPlot(seu, reduction = "umap",
              cols = colours,
              group.by = "predicted.ann_level_4", label = FALSE,
              label.size = 2.5, repel = TRUE, label.box = FALSE) +
  #NoLegend() +
  ggtitle("predicted.ann_level_4")  +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        title = element_text(size = 9))
p3 <- DimPlot(seu, reduction = "umap",
              cols = colours,
              group.by = "predicted.ann_finest_level", label = FALSE,
              label.size = 2.5, repel = TRUE, label.box = FALSE) +
  #NoLegend() +
  ggtitle("predicted.ann_finest_level")  +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        title = element_text(size = 9))

p4 <- DimPlot(seu, reduction = "umap",
              cols = colours,
              group.by = "SCT_snn_res.0.9", label = TRUE,
              label.size = 2.5, repel = TRUE, label.box = TRUE) +
  NoLegend() +
  ggtitle("SCT_snn_res.0.9")  +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        title = element_text(size = 9))
#+ NoLegend()
ggsave(plot = p1/p2 + p3/p4,
       file=glue(dirPLOT,capture,".HLCA.umap.231106.jpeg"), device="jpeg", dpi=600, height=12, width=12)

print("HLCA v2 annotation Done!")


