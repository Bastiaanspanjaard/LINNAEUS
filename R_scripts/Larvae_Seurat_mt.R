# Description ####
# Use Seurat to analyze the larvae dataset

# Dependencies ####
library(Seurat)
library(dplyr)

# Load data and create object ####
L1.data <- Read10X("../2017_10X_1/Z1_wt_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
L1 <- CreateSeuratObject(L1.data, project = "L1")

L21.data <- Read10X("../Z2_1_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
L21 <- CreateSeuratObject(L21.data, project = "L21")
Larvae <- MergeSeurat(object1 = L1, object2 = L21, add.cell.id1 = "L1",
                      add.cell.id2 = "L21", project = "Larvae")

L22.data <- Read10X("../Z2_2_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
L22 <- CreateSeuratObject(L22.data, project = "L22")
Larvae <- MergeSeurat(object1 = Larvae, object2 = L22, 
                      add.cell.id2 = "L22", project = "Larvae")

L3.data <- Read10X("../Z3_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
L3 <- CreateSeuratObject(L3.data, project = "L3")
Larvae <- MergeSeurat(object1 = Larvae, object2 = L3,
                      add.cell.id2 = "L3", project = "Larvae")

L4.data <- Read10X("../Z4_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
L4 <- CreateSeuratObject(L4.data, project = "L4")
Larvae <- MergeSeurat(object1 = Larvae, object2 = L4,
                      add.cell.id2 = "L4", project = "Larvae")

L5.data <- Read10X("../Z5_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
L5 <- CreateSeuratObject(L5.data, project = "L5")
Larvae <- MergeSeurat(object1 = Larvae, object2 = L5,
                      add.cell.id2 = "L5", project = "Larvae")

F11.data <- Read10X("../F1_1_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
F11 <- CreateSeuratObject(F11.data, project = "F11")
Larvae <- MergeSeurat(object1 = Larvae, object2 = F11,
                      add.cell.id2 = "F11", project = "Larvae")

F12.data <- Read10X("../F1_2_mt/outs/filtered_gene_bc_matrices/danrer10Zb_mt/")
F12 <- CreateSeuratObject(F12.data, project = "F12")
Larvae <- MergeSeurat(object1 = Larvae, object2 = F12,
                      add.cell.id2 = "F12", project = "Larvae",
                      min.cells = 3, min.genes = 200)

table(Larvae@meta.data$orig.ident)
# str(Larvae@raw.data)
# table(sapply(Larvae@cell.names, function(x) unlist(strsplit(x, "_"))[1]))
rm(L1, L1.data, L21, L21.data, L22, L22.data, L3, L3.data, L4, L4.data, L5, L5.data, F11, F11.data, F12, F12.data)

# Quality control ####
mito.genes <- grep(pattern = "^mt-", x = rownames(x = Larvae@data), value = TRUE)
percent.mito <- Matrix::colSums(Larvae@raw.data[mito.genes, ])/Matrix::colSums(Larvae@raw.data)
Larvae <- AddMetaData(object = Larvae, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = Larvae, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
Larvae <- FilterCells(object = Larvae, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.075))
length(x = Larvae@cell.names)

# cellnames <- data.frame(Larvae@cell.names, stringsAsFactors = F)
# cellnames$Barcode <- sapply(cellnames$Larvae.cell.names,
#                             function(x) unlist(strsplit(x, "_"))[2])
# cellnames$Library <- sapply(cellnames$Larvae.cell.names,
#                             function(x) unlist(strsplit(x, "_"))[1])
# write.csv(cellnames, "./Larvae_cell_names.csv", quote = F, row.names = F)

# Normalization ####
Larvae <- NormalizeData(object = Larvae, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Variable genes ####
Larvae <- FindVariableGenes(object = Larvae, mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Larvae@var.genes) # 2779

# Scaling data and regressing out ####
Larvae <- ScaleData(object = Larvae, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))

# PCA ####
Larvae <- RunPCA(object = Larvae, pc.genes = Larvae@var.genes, pcs.compute = 100,
                 do.print = F)

# TSNE and plot ####
Larvae <- RunTSNE(object = Larvae, dims.use = 1:30, do.fast = TRUE)
# png("Larvae_PCA_60_clustered.png")
TSNEPlot(object = Larvae, pt.size = 0.1)
# dev.off()
# 60, 50, 40, 70 - there are no differences so 60 seems fine as we had this before.
# Extract tSNE + cell type data
tSNE.data <- data.frame(Larvae@dr$tsne@cell.embeddings)
tSNE.data$Cell <- rownames(tSNE.data)
cell.type.data <- data.frame(Cluster = Larvae@ident)
cell.type.data$Cell <- rownames(cell.type.data)

cell.data <- merge(tSNE.data, cell.type.data)
cell.data$Library <-
  sapply(cell.data$Cell,
         function(x) unlist(strsplit(x, "_"))[1])
cell.data$Barcode <-
  sapply(cell.data$Cell,
         function(x) unlist(strsplit(x, "_"))[2])

cell.types <- read.table("/local/Bo/Larvae_seurat/Larvae.celltypes.final.csv",
                       stringsAsFactors = F, sep = ";")[1:3, -1]
cell.types <- data.frame(t(cell.types))
colnames(cell.types) <- c("Cell.type", "Category", "Cluster")

cell.data <- merge(cell.data, cell.types)
cell.data.counts <- data.frame(table(cell.data$Library, cell.data$Cell.type))
# write.csv(cell.data, file = "./Larvae_Seurat_batch_r_out_cells_3.csv",
#           quote = F, row.names = F)

# Cluster ####
Larvae <- FindClusters(object = Larvae, reduction.type = "pca", dims.use = 1:30, 
                     resolution = 1.8, print.output = 0, save.SNN = TRUE, temp.file.location = "./", force.recalc = T)
# save(Larvae, file = "./Larvae_Seurat_batch_r_out_2.Robj")
# save(Larvae, file = "./Larvae_Seurat_batch_r_out_2_30PC.Robj")
load(file = "./Larvae_Seurat_batch_r_out_2.Robj")

# Find differentially expressed genes ####
negbinom.markers <-
  FindAllMarkers(Larvae, min.pct = 0.25, only.pos = T,
                 test.use = "negbinom")
# Cluster 72 has no differentially expressed genes since it has < 3 cells.

