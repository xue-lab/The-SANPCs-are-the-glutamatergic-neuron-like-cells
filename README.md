# The-SANPCs-are-the-glutamatergic-neuron-like-cells

The sinoatrial node pacemaker cells are the glutamatergic neuron-like cells that reside in the heart

library(Seurat)
library(magrittr)
library(dplyr)
library(cowplot)

SANPC.data <- read.table("/expression_count.txt", header=T, row.names = 1)
CC.data <- read.table("/expression_count.txt", header=T, row.names = 1)

SANPC <- CreateSeuratObject(raw.data = SANPC.data, project = "SC_SANPC", min.cells = 3, min.genes = 1000)
SANPC@meta.data$Group <- "SANPC"
SANPC <- FilterCells(SANPC, subset.names = "nGene", low.thresholds = 1000, high.thresholds = Inf)
SANPC <- NormalizeData(SANPC)
SANPC <- ScaleData(SANPC, display.progress = F)

CC <- CreateSeuratObject(raw.data = CC.data, project = "SC_CC", min.cells = 3, min.genes = 1000)
CC@meta.data$Group <- "CC"
CC <- FilterCells(CC, subset.names = "nGene", low.thresholds = 1000, high.thresholds = 12500)
CC <- NormalizeData(CC)
CC <- ScaleData(CC, display.progress = F)

SANPC <- FindVariableGenes(SANPC, do.plot = F)
CC <- FindVariableGenes(CC, do.plot = F)
g.1 <- head(rownames(SANPC@hvg.info), 2000)
g.2 <- head(rownames(CC@hvg.info), 2000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(SANPC@scale.data))
genes.use <- intersect(genes.use, rownames(CC@scale.data))

SANPC_CC.combined <- RunCCA(SANPC, CC, genes.use = genes.use, num.cc = 30)

SANPC_CC.combined <- AlignSubspace(SANPC_CC.combined, reduction.type = "cca", grouping.var = "Group", 
                                  dims.align = 1:15)

SANPC_CC.combined <- RunTSNE(SANPC_CC.combined, reduction.use = "cca.aligned", dims.use = 1:15, 
                            do.fast = T)
SANPC_CC.combined <- FindClusters(SANPC_CC.combined, reduction.type = "cca.aligned", 
                                 resolution = 0.6, dims.use = 1:15)

