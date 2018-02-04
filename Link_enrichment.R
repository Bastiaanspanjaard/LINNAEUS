# Description ####
# Calculate enrichment of scar links between cell types and cluster cell types
# based on that.

# output.fish <- "Z2"

# Dependencies ####
source("./Scripts/linnaeus-scripts/scar_helper_functions.R")

# Parameters ####
# Maximum scar probability to include scar in tree building
max.scar.p <- 0.001
# Maximum number of embryos a scar can be present in to include in tree building
max.larvae <- 1

# Load data ####
tsne.coord <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells_2.csv",
                       stringsAsFactors = F)

# Colors
larvae.colors <- read.csv("./Data/color_table_larvae.csv",
                          stringsAsFactors = F)
colnames(larvae.colors)[2] <- "Cell.type"
colors.use <- larvae.colors$color
names(colors.use) <- larvae.colors$Cell.type
tsne.coord$Cell.type <- factor(tsne.coord$Cell.type,
                               levels = larvae.colors$Cell.type[order(larvae.colors$cluster)])

# Prepare cell type labels
# cell.types <- data.frame(Cluster = unique(tsne.coord$Cluster),
#                          Cell.type = unique(tsne.coord$Cluster))
label.positions <-
  aggregate(tsne.coord[, c("tSNE_1", "tSNE_2")],
            by = list(Cell.type = tsne.coord$Cell.type),
            median)
larvae.colors <- merge(larvae.colors, label.positions, sort = F)

# Prep scar data ####
scar.input <- 
  # read.csv("./Data/2017_10X_1/Z1_scars_compared.csv", stringsAsFactors = F)
  # scar.input$Cell <- paste("L1", scar.input$Barcode, sep = "_")
  # read.csv("./Data/2017_10X_2/Z2_scars_compared.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_2/Z3_scars_compared.csv", stringsAsFactors = F)
  # scar.input$Cell <- paste("L3", scar.input$Barcode, sep = "_")
  read.csv("./Data/2017_10X_10_CR/Z4_scars_compared.csv", stringsAsFactors = F)
  scar.input$Cell <- paste("L4", scar.input$Barcode, sep = "_")
  # read.csv("./Data/2017_10X_10_CR/Z5_scars_compared.csv", stringsAsFactors = F)
# scar.input$Cell <- paste("L5", scar.input$Barcode, sep = "_")

# Scar data preparation ####
scar.input <- merge(scar.input[, c("Cell", "Scar", "Presence", "p")],
                    tsne.coord)
scars <- 
      scar.input[scar.input$p <= max.scar.p & scar.input$Presence <= max.larvae, 
                 c("Cell", "Scar", "Cell.type")]
scars <- scars[!duplicated(scars), ]

# Calculate and visualize connections ####
# Calculate adjacency matrix of cells with p-value of connections given by
# their shared scars.
connections.all <- merge(scars, scars, by.x = "Scar", by.y = "Scar")
connections.all <- unique(connections.all[, -1])
colnames(connections.all) <- c("Cell.1", "Cell.type.1", "Cell.2", "Cell.type.2")
connections <- data.frame(t(combn(unique(scars$Cell), 2)))
colnames(connections) <- c("Cell.1", "Cell.2")
connections <- 
  merge(connections, connections.all)
connections <- 
  merge(connections, 
        tsne.coord[, c("Cell", "tSNE_1", "tSNE_2")], 
        by.x = "Cell.1", by.y = "Cell")
colnames(connections)[5:6] <- c("Cell.1.tSNE_1", "Cell.1.tSNE_2")
connections <- 
  merge(connections, 
        tsne.coord[, c("Cell", "tSNE_1", "tSNE_2")], 
        by.x = "Cell.2", by.y = "Cell")
colnames(connections)[7:8] <- c("Cell.2.tSNE_1", "Cell.2.tSNE_2")

# Plot with cell type labels in figure is commented out since it takes a lot of 
# time.
# png("./Images/2017_10X_2/Larvae_cell_connections_Z2_paper.png",
#     width = 8, height = 6, units = "in", res = 300)
# ggplot() +
#   geom_segment(data = connections,
#                aes(x = Cell.1.tSNE_1, xend = Cell.2.tSNE_1,
#                    y = Cell.1.tSNE_2, yend = Cell.2.tSNE_2), alpha = 0.005) +
#   geom_point(data = tsne.coord, size = 0.01,
#              aes(x = tSNE_1, y = tSNE_2, color = Cell.type)) +
#   scale_color_manual(values = colors.use) +
#   labs(color = "", x = "", y = "") +
#   guides(color = F) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         text = element_text(size = 8),
#         axis.ticks = element_blank(),
#         axis.text = element_blank())
# dev.off()

# Calculate link enrichment between cell types ####
# Connections per cell type
cell.type.link <- 
  data.frame(table(connections$Cell.type.1) + table(connections$Cell.type.2))
colnames(cell.type.link) <- c("Cell.type", "Total.connections")
# We are counting the endpoints here, so the total number of connections is c.total/2
c.total <- sum(cell.type.link$Total.connections)
cell.type.link$Connection.fraction <- cell.type.link$Total.connections/c.total

# Connections per cell type combination
cell.type.pair.link <-
  data.frame(table(connections$Cell.type.1, connections$Cell.type.2))
cell.type.pair.link <-
  acast(cell.type.pair.link, Var1 ~ Var2, value.var = "Freq")
cell.type.pair.link[lower.tri(cell.type.pair.link)] <-
  cell.type.pair.link[lower.tri(cell.type.pair.link)] +
  t(cell.type.pair.link)[lower.tri(t(cell.type.pair.link))]
cell.type.pair.link[upper.tri(cell.type.pair.link)] <- -1
cell.type.pair.link <-
  melt(cell.type.pair.link)
colnames(cell.type.pair.link) <- c("Cell.type.1", "Cell.type.2", "Connections")
cell.type.pair.link <- cell.type.pair.link[cell.type.pair.link$Connections > -1, ]

# Calculate expected number of links between cell types -
# the chance of making a connection between cell types A and B is
# P(A-B) = 2 * C(A)*C(B)/(C(tot))^2 for A!=B, and
# P(A-A) = c(A) * C(A)/(C(tot))^2.
cell.type.pair.link <- merge(cell.type.pair.link, cell.type.link,
                             by.x = "Cell.type.1", by.y = "Cell.type")
colnames(cell.type.pair.link)[4:5] <- 
  paste(colnames(cell.type.pair.link)[4:5], "1", sep = ".")
cell.type.pair.link <- merge(cell.type.pair.link, cell.type.link,
                             by.x = "Cell.type.2", by.y = "Cell.type")
colnames(cell.type.pair.link)[6:7] <- 
  paste(colnames(cell.type.pair.link)[6:7], "2", sep = ".")
cell.type.pair.link$P <- NA
cell.type.pair.link$P <-
  2 * cell.type.pair.link$Connection.fraction.1 * cell.type.pair.link$Connection.fraction.2
cell.type.pair.link$P[cell.type.pair.link$Cell.type.1 == cell.type.pair.link$Cell.type.2] <-
  1/2 * cell.type.pair.link$P[cell.type.pair.link$Cell.type.1 == cell.type.pair.link$Cell.type.2]
cell.type.pair.link$E.connections <- cell.type.pair.link$P * c.total/2 # Connections is c.total/2
cell.type.pair.link$SD <- 
  sqrt(0.5 * c.total * cell.type.pair.link$P * (1 - cell.type.pair.link$P))
cell.type.pair.link$z <- 
  (cell.type.pair.link$Connections - cell.type.pair.link$E.connections)/cell.type.pair.link$SD
cell.type.pair.link$z[is.na(cell.type.pair.link$z)] <- 0
cell.type.pair.link$binom.p <-
  apply(as.matrix(cell.type.pair.link[, c("Connections", "P")]), 1,
        function(x){
          bt <- binom.test(x[1], c.total/2, x[2], alternative = "two.sided")
          return(bt$p.value)
        })
cell.type.pair.link$p.adj <- p.adjust(cell.type.pair.link$binom.p, method = "fdr")

plot.type.link <- cell.type.pair.link[cell.type.pair.link$p.adj < 0.01 & cell.type.pair.link$z > 0,
                                      c("Cell.type.1", "Cell.type.2")]
plot.type.link <- plot.type.link[plot.type.link$Cell.type.1 != plot.type.link$Cell.type.2, ]

plot.type.link <- merge(plot.type.link, larvae.colors[, c("Cell.type", "tSNE_1", "tSNE_2")],
                        by.x = "Cell.type.1", by.y = "Cell.type")
colnames(plot.type.link)[3:4] <- paste(colnames(plot.type.link)[3:4], "type.1", sep = "_")
plot.type.link <- merge(plot.type.link, larvae.colors[, c("Cell.type", "tSNE_1", "tSNE_2")],
                        by.x = "Cell.type.2", by.y = "Cell.type")
colnames(plot.type.link)[5:6] <- paste(colnames(plot.type.link)[5:6], "type.2", sep = "_")

# Plot on tSNE
# png("./Images/2017_10X_2/Larvae_enriched_connections_Z2_paper.png",
#     width = 8, height = 6, units = "in", res = 300)
# ggplot() +
#   geom_point(data = tsne.coord, size = 0.01,
#              aes(x = tSNE_1, y = tSNE_2, color = Cell.type)) +
#   geom_segment(data = plot.type.link, color = "darkgrey", size = 1, alpha = 0.75,
#                aes(x = tSNE_1_type.1, y = tSNE_2_type.1, 
#                    xend = tSNE_1_type.2, yend = tSNE_2_type.2)) +
#   scale_color_manual(values = colors.use) +
#   labs(color = "", x = "", y = "") +
#   guides(color = F) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         text = element_text(size = 8),
#         axis.ticks = element_blank(),
#         axis.text = element_blank())
# dev.off()

# Calculate and plot heatmap
type.connection.heat <-
  acast(cell.type.pair.link, Cell.type.1 ~ Cell.type.2, value.var = "z")
type.connection.heat[upper.tri(type.connection.heat)] <-
  t(type.connection.heat)[upper.tri(t(type.connection.heat))]
full.zeroes <- 
  apply(type.connection.heat, 1,
        function(x) sum(x != 0) == 0)
type.connection.heat <- 
  type.connection.heat[!full.zeroes, !full.zeroes]

type.connection.dist <- type.connection.heat - min(type.connection.heat)
type.connection.dist <- 1 - type.connection.dist/max(type.connection.dist)
type.connection.dist <- as.dist(type.connection.dist)
enr.hclust <- hclust(type.connection.dist, method = "average") 
# Has highest bootstrap stability
enr.hclust.dg <- as.dendrogram(enr.hclust)
# 
# require(fpc)
# clusterboot(type.connection.dist, B = 100,
#                  clustermethod = disthclustCBI, k = 4, cut = "number",
#                  method = "average")

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 49)
my_palette <- c(rep("Blue", 125), my_palette, rep("Red", 125))

# pdf("./Images/2017_10X_10/Z4_scar_enrichment_heatmap_binom_av_hclust_figure_resub.pdf",
#            width = 8, height = 5)
gplots::heatmap.2(type.connection.heat, dendrogram = "column", #keysize = 0.3, 
                  Rowv = enr.hclust.dg, Colv = enr.hclust.dg,
                  trace = "none", labCol = "",#key = T, 
                  col = my_palette, margins = c(1, 12),
                  cexRow = 0.5)
# dev.off()
