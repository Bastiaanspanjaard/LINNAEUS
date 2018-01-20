# Description ####
# Determine cell-specific scar detection rates

# Dependencies ####
source("./Scripts/linnaeus-scripts/scar_helper_functions.R")

# Load data ####
# Filtered scars
F11.scars <- read.csv("./Data/2017_10X_10_CR/F1_1_used_scars_7Larvae.csv",
                      stringsAsFactors = F)
F11.scars$Library <- "F11"
F12.scars <- read.csv("./Data/2017_10X_10_CR/F1_2_used_scars_7Larvae.csv",
                      stringsAsFactors = F)
F12.scars$Library <- "F12"
# Scar probabilities?
# Cell types
cell.types <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells.csv",
                       stringsAsFactors = F)
clusters <- read.csv("./Data/Larvae_data/Preliminary_celltypes.csv",
                     stringsAsFactors = F)
colnames(clusters)[1] <- "Cluster"
cell.types <- merge(cell.types, clusters)

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

# N.F11 <- sum(cell.types$Library == "F11")
# N.F12 <- sum(cell.types$Library == "F12")
# 
# scar.freqs$Detection.rate.1 <- scar.freqs$Freq.1/N.F11
# scar.freqs$Detection.rate.2 <- scar.freqs$Freq.2/N.F12
# scar.freqs$Max.det.rate <- 
#   apply(scar.freqs[, c("Detection.rate.1", "Detection.rate.2")], 1, max)
# scar.freqs$Scar <- 
#   apply(scar.freqs[, c("CIGAR", "Max.det.rate")], 1,
#         function(x){
#           if(as.numeric(x[2]) >= 0.01){
#             return(x[1])
#           }else{
#             return("Other")
#           }
#         }
#   )
# scar.freqs.F11 <- aggregate(scar.freqs$Detection.rate.1,
#                             by = list(Scar = scar.freqs$Scar),
#                             sum)
# colnames(scar.freqs.F11)[2] <- "Detection"
# scar.freqs.F11$Larva <- "F1_1"
# scar.freqs.F12 <- aggregate(scar.freqs$Detection.rate.2,
#                             by = list(Scar = scar.freqs$Scar),
#                             sum)
# colnames(scar.freqs.F12)[2] <- "Detection"
# scar.freqs.F12$Larva <- "F1_2"
# scar.freqs.plot <- rbind(scar.freqs.F11, scar.freqs.F12)
# 
# ggplot(scar.freqs.plot) +
#   geom_col(aes(x = Scar, y = Detection)) +
#   facet_wrap(~Larva)
# 
# scar.freqs.s <- 
#   melt(scar.freqs[scar.freqs$Max.det.rate >= 0.01, 
#                   c("CIGAR", "Detection.rate.1", "Detection.rate.2")])
# colnames(scar.freqs.s) <- c("Scar", "Larva", "Detection.rate")


# scar.freqs$Perc.1 <- 100 * scar.freqs$Freq.1/sum(scar.freqs$Freq.1, na.rm = T)
# scar.freqs$Perc.2 <- 100 * scar.freqs$Freq.2/sum(scar.freqs$Freq.2, na.rm = T)
scar.freqs$Rate.1 <- scar.freqs$Freq.1/sum(scar.freqs$Freq.1, na.rm = T)
scar.freqs$Rate.2 <- scar.freqs$Freq.2/sum(scar.freqs$Freq.2, na.rm = T)
# scar.freqs$Mean.perc <- rowMeans(scar.freqs[, c("Perc.1", "Perc.2")])
scar.freqs$Max.rate <- 
    apply(scar.freqs[, c("Rate.1", "Rate.2")], 1, max)
# scar.freqs <- scar.freqs[order(-scar.freqs$Mean.perc), ]
scar.freqs <- scar.freqs[order(-scar.freqs$Max.rate), ]
scar.freqs$Scar <- paste(1:nrow(scar.freqs), scar.freqs$CIGAR, sep = ":")
# Make cutoff and piechart of 'real' scars
# scar.freqs$Scar.1 <-
#   ifelse(apply(scar.freqs[, c("Perc.1", "Perc.2")], 1, max) >= 0.5,
#          scar.freqs$Scar, "Other")
scar.freqs$Scar.1 <-
  ifelse(scar.freqs$Max.det.rate >= 0.01,
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
# celltype.counts.F11 <- celltype.counts[celltype.counts$Library == "F11", ]
# celltype.counts.F11 <- celltype.counts.F11[order(celltype.counts.F11$Cells), ]
celltype.all.counts <- celltype.all.counts[order(celltype.all.counts$Cells), ]
celltype.scar.detected.F11$Cell.type <-
  factor(celltype.scar.detected.F11$Cell.type, levels = celltype.all.counts$Cell.type)
# celltype.scar.detected.F11$Det.min.95 <-
#   apply(celltype.scar.detected.F11[, 4:5], 1,
#         function(x){
#           xt <- prop.test(as.integer(x[1]), as.integer(x[2]))
#           return(xt$conf.int[1])
#         }
#   )
# celltype.scar.detected.F11$Det.max.95 <-
#   apply(celltype.scar.detected.F11[, 4:5], 1,
#         function(x){
#           xt <- prop.test(as.integer(x[1]), as.integer(x[2]))
#           return(xt$conf.int[2])
#         }
#   )

# pdf("./Images/2017_10X_10/F1_1_detection_rates_2.pdf",
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
# ggplot(celltype.scar.detected.F11) +
#   geom_bar(stat = "identity", position = "dodge",
#            aes(x = F1.ident, fill = Scar, y = Detection)) +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0))
# F11.detection.examples <- 
#   celltype.scar.detected.F11[celltype.scar.detected.F11$F1.ident %in% 
#                                c("Fibroblasts A", "Skeletal muscle A", "Brain neurons"), ]
# pdf("./Images/2017_10X_10/F1_1_detection_rates_selection.pdf",
#         width = 12, height = 9)
# ggplot(F11.detection.examples,
#        aes(x = F1.ident, fill = Scar)) +
#   geom_bar(stat = "identity", position = "dodge", aes(y = Detection)) +
#   geom_errorbar(aes(ymin = Det.min.95, ymax = Det.max.95), position = "dodge") +
#   theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.25)) +
#   labs(x = "")
# dev.off()

celltype.scar.detected.F12 <- 
  celltype.scar.detected[celltype.scar.detected$Library == "F12", ]
celltype.scar.detected.F12 <-
  celltype.scar.detected.F12[celltype.scar.detected.F12$Scar %in% 
                               c("1:47M6D28M", "2:47M3D28M", "4:49M4I2M4D20M",
                                 "5:48M4D27M", "8:49M1I25M", "9:Other"), ]
# celltype.counts.F12 <- celltype.counts[celltype.counts$Library == "F12", ]
# celltype.counts.F12 <- celltype.counts.F12[order(celltype.counts.F12$Cells), ]
celltype.scar.detected.F12$Cell.type <-
  factor(celltype.scar.detected.F12$Cell.type, levels = celltype.all.counts$Cell.type)
# celltype.scar.detected.F12$Det.min.95 <-
#   apply(celltype.scar.detected.F12[, 4:5], 1,
#         function(x){
#           xt <- prop.test(as.integer(x[1]), as.integer(x[2]))
#           return(xt$conf.int[1])
#         }
#   )
# celltype.scar.detected.F12$Det.max.95 <-
#   apply(celltype.scar.detected.F12[, 4:5], 1,
#         function(x){
#           xt <- prop.test(as.integer(x[1]), as.integer(x[2]))
#           return(xt$conf.int[2])
#         }
#   )

# pdf("./Images/2017_10X_10/F1_2_detection_rates_2.pdf",
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

# pdf("./Images/2017_10X_10/F1_detection_rates_2.pdf",
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
sum(F11.high$Cells) # 2000 -> 15%
sum(F11.medium$Cells) # 3000 > 25%
sum(F11.low$Cells) # 7000 -> 60%
weighted.mean(F11.high$Detection, F11.high$Cells) # ~ 70%
weighted.mean(F11.medium$Detection, F11.medium$Cells) # ~ 30%
weighted.mean(F11.low$Detection, F11.low$Cells) # ~ 15%

F11.high.8 <- 
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F11$F1.ident %in% F11.high$F1.ident, ]
F11.medium.8 <- 
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F11$F1.ident %in% F11.medium$F1.ident, ]
F11.low.8 <-   
  celltype.scar.detected.F11[celltype.scar.detected.F11$Scar == "8:49M1I25M" &
                               celltype.scar.detected.F11$F1.ident %in% F11.low$F1.ident, ]
weighted.mean(F11.high.8$Detection, F11.high.8$Cells) # ~ 11%
weighted.mean(F11.medium.8$Detection, F11.medium.8$Cells) # ~ 1.5%
weighted.mean(F11.low.8$Detection, F11.low.8$Cells) # ~ 0.4%
