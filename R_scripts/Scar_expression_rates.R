# Description ####
# Calculate average fraction of cells with scars per cell type, and average
# number of scars for non-zero cells per cell type.

# Larvae ####
tsne.coord.L <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells_2.csv",
                         stringsAsFactors = F)
cell.types <- unique(tsne.coord.L[, c("Cluster", "Cell.type")])
tsne.coord.L$Cell.type <-
  factor(tsne.coord.L$Cell.type, levels = cell.types$Cell.type[order(cell.types$Cluster)])
scars.Z1 <- read.csv("./Data/2017_10X_1/Z1_scars_compared.csv", stringsAsFactors = F)
scars.Z1$Cell <- paste("L1", scars.Z1$Barcode, sep = "_")
scars.Z2 <- read.csv("./Data/2017_10X_2/Z2_scars_compared.csv", stringsAsFactors = F)
scars.Z3 <- read.csv("./Data/2017_10X_2/Z3_scars_compared.csv", stringsAsFactors = F)
scars.Z3$Cell <- paste("L3", scars.Z3$Barcode, sep = "_")
scars.Z4 <- read.csv("./Data/2017_10X_10_CR/Z4_scars_compared.csv", stringsAsFactors = F)
scars.Z4$Cell <- paste("L4", scars.Z4$Barcode, sep = "_")
scars.Z5 <- read.csv("./Data/2017_10X_10_CR/Z5_scars_compared.csv", stringsAsFactors = F)
scars.Z5$Cell <- paste("L5", scars.Z5$Barcode, sep = "_")
scars.F11 <- read.csv("./Data/2017_10X_10_CR/F1_1_used_scars_7Larvae.csv",
                      stringsAsFactors = F)
scars.F11$Cell <- paste("F11", scars.F11$Barcode, sep = "_")
scars.F12 <- read.csv("./Data/2017_10X_10_CR/F1_2_used_scars_7Larvae.csv",
                      stringsAsFactors = F)
scars.F12$Cell <- paste("F12", scars.F12$Barcode, sep = "_")

scars.L <- 
  rbind(scars.Z1[, c("Cell", "Scar", "Sequence")], 
        scars.Z2[, c("Cell", "Scar", "Sequence")],
        scars.Z3[, c("Cell", "Scar", "Sequence")], 
        scars.Z4[, c("Cell", "Scar", "Sequence")],
        scars.Z5[, c("Cell", "Scar", "Sequence")])

scars.F11.2 <- merge(scars.F11, unique(scars.L[, c("Scar", "Sequence")]), all.x = T)
scars.F11.3 <- scars.F11.2[is.na(scars.F11.2$Scar), ]
scars.F11.4 <- aggregate(scars.F11.3$Barcode,
                         by = list(Sequence = scars.F11.3$Sequence,
                                   CIGAR = scars.F11.3$CIGAR),
                         length)
colnames(scars.F11.4)[3] <- "Freq"
scars.F11.4 <- scars.F11.4[order(-scars.F11.4$Freq), ]
scars.F11.4$Order <- 1:nrow(scars.F11.4)
scars.F11.4$Scar <- paste("F11_", scars.F11.4$Order, ":", scars.F11.4$CIGAR, sep = "")
scars.F11.3 <- merge(scars.F11.3[, -8], scars.F11.4[, c(1, 2, 5)])
scars.F11.1 <- rbind(scars.F11.2[!is.na(scars.F11.2$Scar), ],
                     scars.F11.3)
scars.F12.2 <- merge(scars.F12, unique(scars.L[, c("Scar", "Sequence")]), all.x = T)
scars.F12.3 <- scars.F12.2[is.na(scars.F12.2$Scar), ]
scars.F12.4 <- aggregate(scars.F12.3$Barcode,
                         by = list(Sequence = scars.F12.3$Sequence,
                                   CIGAR = scars.F12.3$CIGAR),
                         length)
colnames(scars.F12.4)[3] <- "Freq"
scars.F12.4 <- scars.F12.4[order(-scars.F12.4$Freq), ]
scars.F12.4$Order <- 1:nrow(scars.F12.4)
scars.F12.4$Scar <- paste("F12_", scars.F12.4$Order, ":", scars.F12.4$CIGAR, sep = "")
scars.F12.3 <- merge(scars.F12.3[, -8], scars.F12.4[, c(1, 2, 5)])
scars.F12.1 <- rbind(scars.F12.2[!is.na(scars.F12.2$Scar), ],
                     scars.F12.3)

scars.L <- rbind(scars.L, scars.F11.1[, c("Cell", "Scar", "Sequence")], 
                 scars.F12.1[, c("Cell", "Scar", "Sequence")])
scars.L <- merge(scars.L, tsne.coord.L[, c("Cell", "Cell.type")])
cells.w.scars.L <- unique(scars.L[, c("Cell", "Cell.type")])

cell.type.counts.L <- data.frame(table(tsne.coord.L$Cell.type))
colnames(cell.type.counts.L) <- c("Cell.type", "Total")
cells.w.scars.counts.L <- data.frame(table(cells.w.scars.L$Cell.type))
colnames(cells.w.scars.counts.L) <- c("Cell.type", "Cells.w.scars")
celltype.scar.counts.L <- data.frame(table(scars.L$Cell.type))
colnames(celltype.scar.counts.L) <- c("Cell.type", "Scars")

stats.L <- merge(cell.type.counts.L, merge(cells.w.scars.counts.L, celltype.scar.counts.L))
stats.L$Perc.at.least.one <- stats.L$Cells.w.scars/stats.L$Total
stats.L$Av.num.scar <- stats.L$Scars/stats.L$Cells.w.scars
stats.L <- stats.L[order(stats.L$Cell.type), ]
# write.csv(stats.L, file = "./Data/Larvae_data/Scar_detection_stats.csv",
#           quote = F, row.names = F)

# Adults ####
tsne.coord.in.1 <- 
  read.csv("./Data/Adult_data/Adults567_brain_cells.csv",
           stringsAsFactors = F)
tsne.coord.in.1$Origin <- "brain"
  # paste(tsne.coord.in.1$Cell.type, "brain")
tsne.coord.in.2 <- 
  read.csv("./Data/Adult_data/Adults567_heart_cells.csv",
           stringsAsFactors = F)
tsne.coord.in.2$Origin <- "heart"
  # paste(tsne.coord.in.2$Cell.type, "heart")
tsne.coord.in.3 <- 
  read.csv("./Data/Adult_data/Adults567_pancreas_cells.csv",
           stringsAsFactors = F)
tsne.coord.in.3$Origin <- "pancreas"
  # paste(tsne.coord.in.3$Cell.type, "pancreas")
tsne.coord.A <- rbind(tsne.coord.in.1, tsne.coord.in.2, tsne.coord.in.3)

# cell.types <- unique(tsne.coord.L[, c("Cluster", "Cell.type")])
# tsne.coord.L$Cell.type <-
#   factor(tsne.coord.L$Cell.type, levels = cell.types$Cell.type[order(cell.types$Cluster)])
scars.A5 <- read.csv("./Data/2017_10X_7/A5_scars_compared.csv", stringsAsFactors = F)
scars.A6 <- read.csv("./Data/2017_10X_6/A6_scars_compared.csv", stringsAsFactors = F)
scars.A7 <- read.csv("./Data/2018_10X_1/A7_scars_compared.csv", stringsAsFactors = F)

scars.A <- 
  rbind(scars.A5[, c("Cell", "Scar", "Sequence")], 
        scars.A6[, c("Cell", "Scar", "Sequence")],
        scars.A7[, c("Cell", "Scar", "Sequence")])

scars.A <- merge(scars.A, tsne.coord.A[, c("Cell", "Cell.type", "Origin")])
cells.w.scars.A <- unique(scars.A[, c("Cell", "Cell.type", "Origin")])

cell.type.counts.A <- data.frame(table(tsne.coord.A$Cell.type, tsne.coord.A$Origin))
colnames(cell.type.counts.A) <- c("Cell.type", "Origin", "Total")
cell.type.counts.A <- cell.type.counts.A[cell.type.counts.A$Total > 0, ]

cells.w.scars.counts.A <- data.frame(table(cells.w.scars.A$Cell.type, cells.w.scars.A$Origin))
colnames(cells.w.scars.counts.A) <- c("Cell.type", "Origin", "Cells.w.scars")
celltype.scar.counts.A <- data.frame(table(scars.A$Cell.type, scars.A$Origin))
colnames(celltype.scar.counts.A) <- c("Cell.type", "Origin", "Scars")

stats.A <- merge(cell.type.counts.A, merge(cells.w.scars.counts.A, celltype.scar.counts.A))
stats.A$Perc.at.least.one <- stats.A$Cells.w.scars/stats.A$Total
stats.A$Av.num.scar <- stats.A$Scars/stats.A$Cells.w.scars
# stats.L <- stats.L[order(stats.L$Cell.type), ]
# write.csv(stats.A, file = "./Data/Adult_data/Scar_detection_stats.csv",
#           quote = F, row.names = F)

