# Description ####
# Determine cell-specific scar detection rates

# Dependencies ####
source("./scar_helper_functions.R")

# Load data ####
# Filtered scars
F11.scars <- read.csv("./Data/2017_10X_10_CR/F1_1_used_scars_7Larvae.csv",
                      stringsAsFactors = F)
F11.scars$Library <- "F11"
F12.scars <- read.csv("./Data/2017_10X_10_CR/F1_2_used_scars_7Larvae.csv",
                      stringsAsFactors = F)
F12.scars$Library <- "F12"

# Cell types
cell.types <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells_2.csv",
                       stringsAsFactors = F)

# Determine scar content of F1 datasets ####
# Filter data for duplicates and number of reads
F11.scars.nd <- F11.scars[!duplicated(F11.scars[, c("Barcode", "Sequence", "Reads")]), ]
F12.scars.nd <- F12.scars[!duplicated(F12.scars[, c("Barcode", "Sequence", "Reads")]), ]

# Scar frequencies
F11.scar.freqs <- data.frame(table(F11.scars.nd$Sequence, F11.scars.nd$CIGAR))
colnames(F11.scar.freqs) <- c("Sequence", "CIGAR", "Freq.1")
F11.scar.freqs <- F11.scar.freqs[F11.scar.freqs$Freq > 0, ]
F12.scar.freqs <- data.frame(table(F12.scars.nd$Sequence, F12.scars.nd$CIGAR))
colnames(F12.scar.freqs) <- c("Sequence", "CIGAR", "Freq.2")
F12.scar.freqs <- F12.scar.freqs[F12.scar.freqs$Freq > 0, ]
scar.freqs <- merge(F11.scar.freqs, F12.scar.freqs, all = T)
scar.freqs[is.na(scar.freqs)] <- 0

scar.freqs$Rate.1 <- scar.freqs$Freq.1/sum(scar.freqs$Freq.1, na.rm = T)
scar.freqs$Rate.2 <- scar.freqs$Freq.2/sum(scar.freqs$Freq.2, na.rm = T)
scar.freqs$Max.rate <- 
    apply(scar.freqs[, c("Rate.1", "Rate.2")], 1, max)
scar.freqs <- scar.freqs[order(-scar.freqs$Max.rate), ]
scar.freqs$Scar <- paste(1:nrow(scar.freqs), scar.freqs$CIGAR, sep = ":")
scar.freqs$Scar.1 <-
  ifelse(scar.freqs$Max.rate >= 0.01,
         scar.freqs$Scar, "Other")
scar.freqs$Scar.1[scar.freqs$Scar.1 == "Other"] <- 
  paste(sum(scar.freqs$Scar.1 != "Other") + 1, ":Other", sep = "")
scar.freqs.plot.1 <- aggregate(scar.freqs$Rate.1,
                          by = list(Scar = scar.freqs$Scar.1),
                          sum)
colnames(scar.freqs.plot.1)[2] <- "Rate"
scar.freqs.plot.1$Larva <- "1"
scar.freqs.plot.2 <- aggregate(scar.freqs$Rate.2,
                          by = list(Scar = scar.freqs$Scar.1),
                          sum)
colnames(scar.freqs.plot.2)[2] <- "Rate"
scar.freqs.plot.2$Larva <- "2"
scar.freqs.plot <- rbind(scar.freqs.plot.1, scar.freqs.plot.2)
scar.freqs.plot$Scar.order <- sapply(scar.freqs.plot$Scar, function(x) unlist(strsplit(x, ":"))[1])
scar.freqs.plot$Scar.CIGAR <- sapply(scar.freqs.plot$Scar, function(x) unlist(strsplit(x, ":"))[2])
scar.freqs.plot$Scar.CIGAR <- 
  factor(scar.freqs.plot$Scar.CIGAR, 
         levels = unique(scar.freqs.plot$Scar.CIGAR[order(scar.freqs.plot$Scar.order)]))

# pdf("./Images/2017_10X_10/F1_scar_piechart_2.pdf",
#         width = 8, height = 3)
ggplot(scar.freqs.plot) +
  geom_bar(stat = "identity", aes(x = "", y = Rate, fill = Scar.CIGAR)) +
  facet_wrap(~Larva) +
  coord_polar("y", start = 0, direction = -1) +
  labs(x = "", y = "", fill = "") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        text = element_text(size = 16))
# dev.off()

# Determine scar-specific detection rates ####
cells.real.scars.1 <- 
  merge(cell.types[, c("Library", "Barcode", "Cell.type")], 
        F11.scars.nd[, c("Library", "Barcode", "Sequence")])
cells.real.scars.2 <- 
  merge(cell.types[, c("Library", "Barcode", "Cell.type")], 
        F12.scars.nd[, c("Library", "Barcode", "Sequence")])
cells.real.scars <- rbind(cells.real.scars.1, cells.real.scars.2)
cells.real.scars <- merge(cells.real.scars, scar.freqs[, c("Sequence", "Scar.1")])

celltype.counts <- data.frame(table(cell.types$Library, cell.types$Cell.type))
colnames(celltype.counts) <- c("Library", "Cell.type", "Cells")
celltype.all.counts <- data.frame(table(cell.types$Cell.type))
colnames(celltype.all.counts) <- c("Cell.type", "Cells")
celltype.scar.detected <- 
  data.frame(table(cells.real.scars$Library, cells.real.scars$Cell.type,
                                           cells.real.scars$Scar.1))
colnames(celltype.scar.detected) <- c("Library", "Cell.type", "Scar", "Count")
celltype.scar.detected <- merge(celltype.scar.detected, celltype.counts)
celltype.scar.detected$Detection <- 
  celltype.scar.detected$Count/celltype.scar.detected$Cells

celltype.scar.detected$Scar.order <- 
  sapply(as.character(celltype.scar.detected$Scar), function(x) unlist(strsplit(x, ":"))[1])
celltype.scar.detected$Scar.CIGAR <- 
  sapply(as.character(celltype.scar.detected$Scar), function(x) unlist(strsplit(x, ":"))[2])
celltype.scar.detected$Scar.CIGAR <- 
  factor(celltype.scar.detected$Scar.CIGAR, 
         levels = unique(celltype.scar.detected$Scar.CIGAR[order(celltype.scar.detected$Scar.order)]))


celltype.scar.detected.F11 <- 
  celltype.scar.detected[celltype.scar.detected$Library == "F11", ]
celltype.scar.detected.F11 <-
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar %in% 
                               c("1:47M6D28M", "2:47M3D28M", "3:45M30S",
                                 "6:46M3D29M", "7:48M1I1M9D25M", "8:49M1I25M",
                                 "9:Other"), ]
celltype.all.counts <- celltype.all.counts[order(celltype.all.counts$Cells), ]
celltype.scar.detected.F11$Cell.type <-
  factor(celltype.scar.detected.F11$Cell.type, levels = celltype.all.counts$Cell.type)

# pdf("./Images/2017_10X_10/F1_1_detection_rates_3.pdf",
#         width = 5, height = 8)
ggplot(celltype.scar.detected.F11) +
  geom_tile(aes(x = Scar.CIGAR, y = Cell.type, fill = Detection)) +
  labs(title = "1", x = "", y = "", fill = "Detection rate") +
  scale_fill_gradient(low = "gray70", high = "red3", limits = c(0, 1),
                      na.value = "gray70") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 16),
        axis.text.y = element_text(size = 8))
# dev.off()

celltype.scar.detected.F12 <- 
  celltype.scar.detected[celltype.scar.detected$Library == "F12", ]
celltype.scar.detected.F12 <-
  celltype.scar.detected.F12[celltype.scar.detected.F12$Scar %in% 
                               c("1:47M6D28M", "2:47M3D28M", "4:49M4I2M4D20M",
                                 "5:48M4D27M", "8:49M1I25M", "9:Other"), ]
celltype.scar.detected.F12$Cell.type <-
  factor(celltype.scar.detected.F12$Cell.type, levels = celltype.all.counts$Cell.type)

# pdf("./Images/2017_10X_10/F1_2_detection_rates_3.pdf",
#     width = 5.5, height = 8)
ggplot(celltype.scar.detected.F12) +
  geom_tile(aes(x = Scar.CIGAR, y = Cell.type, fill = Detection)) +
  labs(title = "2", x = "", y = "", fill = "Detection rate") +
  scale_fill_gradient(low = "gray70", high = "red3", limits = c(0, 1), 
                      na.value = "gray70") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 16),
        axis.text.y = element_text(size = 8))
# dev.off()

celltype.scar.detected.2 <- celltype.scar.detected
celltype.scar.detected.2$Library <- ifelse(celltype.scar.detected.2$Library == "F11", 1, 2)

# pdf("./Images/2017_10X_10/F1_detection_rates_3.pdf",
#     width = 8.5, height = 8)
ggplot(celltype.scar.detected.2) +
  geom_tile(aes(x = Scar.CIGAR, y = Cell.type, fill = Detection)) +
  labs(x = "", y = "", fill = "Detection rate") +
  facet_wrap(~Library) +
  scale_fill_gradient(low = "gray70", high = "red3", limits = c(0, 1), 
                      na.value = "gray70") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text = element_text(size = 16),
        axis.text.y = element_text(size = 8))
# dev.off()


F12.detection.examples <- 
  celltype.scar.detected.F12[celltype.scar.detected.F12$F1.ident %in% 
                               c("Fibroblasts A", "Skeletal muscle A", "Brain neurons"), ]
# pdf("./Images/2017_10X_10/F1_2_detection_rates_selection.pdf",
#         width = 12, height = 9)
ggplot(F12.detection.examples,
       aes(x = F1.ident, fill = Scar)) +
  geom_bar(stat = "identity", position = "dodge", aes(y = Detection)) +
  geom_errorbar(aes(ymin = Det.min.95, ymax = Det.max.95), position = "dodge") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25)) +
  labs(x = "")
# dev.off()

# Realistic detection rates and abundances ####
F11.high <- celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "2:47M3D28M" &
                                         celltype.scar.detected.F11$Detection > 0.5, ]
F11.medium <- celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "2:47M3D28M" &
                                         celltype.scar.detected.F11$Detection <= 0.5 &
                                           celltype.scar.detected.F11$Detection > 0.2, ]
F11.low <- celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "2:47M3D28M" &
                                           celltype.scar.detected.F11$Detection < 0.2, ]
sum(F11.high$Cells, na.rm =T) # 1726 -> 15%
sum(F11.medium$Cells, na.rm = T) # 2708 > 23%
sum(F11.low$Cells, na.rm = T) # 7283 -> 62%
weighted.mean(F11.high$Detection, F11.high$Cells, na.rm = T) # ~ 70%
weighted.mean(F11.medium$Detection, F11.medium$Cells, na.rm = T) # ~ 30%
weighted.mean(F11.low$Detection, F11.low$Cells, na.rm = T) # ~ 15%

F11.high.8 <- 
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F11$Cell.type %in% F11.high$Cell.type, ]
F11.medium.8 <- 
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F11$Cell.type %in% F11.medium$Cell.type, ]
F11.low.8 <-   
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F11$Cell.type %in% F11.low$Cell.type, ]
weighted.mean(F11.high.8$Detection, F11.high.8$Cells) # ~ 10%
weighted.mean(F11.medium.8$Detection, F11.medium.8$Cells) # ~ 1.5%
weighted.mean(F11.low.8$Detection, F11.low.8$Cells) # ~ 0.4%

F12.high <- celltype.scar.detected.F12[celltype.scar.detected.F12$Scar == "4:49M4I2M4D20M" &
                                         celltype.scar.detected.F12$Detection > 0.5, ]
F12.medium <- celltype.scar.detected.F12[celltype.scar.detected.F12$Scar == "4:49M4I2M4D20M" &
                                           celltype.scar.detected.F12$Detection <= 0.5 &
                                           celltype.scar.detected.F12$Detection > 0.2, ]
F12.low <- celltype.scar.detected.F12[celltype.scar.detected.F12$Scar == "4:49M4I2M4D20M" &
                                        celltype.scar.detected.F12$Detection < 0.2, ]
sum(F12.high$Cells, na.rm =T) # 1381 -> 10%
sum(F12.medium$Cells, na.rm = T) # 2365 > 18%
sum(F12.low$Cells, na.rm = T) # 9656 -> 72%
weighted.mean(F12.high$Detection, F12.high$Cells, na.rm = T) # ~ 79%
weighted.mean(F12.medium$Detection, F12.medium$Cells, na.rm = T) # ~ 29%
weighted.mean(F12.low$Detection, F12.low$Cells, na.rm = T) # ~ 10%

F12.high.8 <- 
  celltype.scar.detected.F12[celltype.scar.detected.F12$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F12$Cell.type %in% F12.high$Cell.type, ]
F12.medium.8 <- 
  celltype.scar.detected.F12[celltype.scar.detected.F12$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F12$Cell.type %in% F12.medium$Cell.type, ]
F12.low.8 <-   
  celltype.scar.detected.F12[celltype.scar.detected.F12$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F12$Cell.type %in% F12.low$Cell.type, ]
weighted.mean(F12.high.8$Detection, F12.high.8$Cells) # ~ 6.3%
weighted.mean(F12.medium.8$Detection, F12.medium.8$Cells) # ~ 0.34%
weighted.mean(F12.low.8$Detection, F12.low.8$Cells) # ~ 0.1%

# Based on these, we take conservative abundances 10, 20, 70. These have
# detection rates 75, 30, 10 for high sites. A low site has 5% of that:
# 3.75%, 1.5% and 0.5% for a cumulative detection of ~1%.