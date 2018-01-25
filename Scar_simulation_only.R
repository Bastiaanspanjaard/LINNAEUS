# Description #### 
# Simulate a developmental tree and gradual scar acquisition of cells, and a
# further expansion after scarring.
# Plot the developmental tree with the scarring events, and the resulting
# scar tree.
# Simulate the by definition incomplete sampling of the cells and scars by 
# sequencing.

# Dependencies ####
source("./Scripts/linnaeus-scripts/scar_helper_functions.R")

# Parameters ####
generations <- 3 # Number of developmental generations while scarring takes place
# Note that the first generation consists of one cell, meaning the nth generation
# consists of 2^(n-1) cells.
generations.post <- 9 # Number of cell divisions a single cell undergoes after scarring
# takes place
expansion <- 2^generations.post # Expansion size of a single cell after scarring due to
# [generations.post] divisions.
sites <- 10
scar.speed <- 0.15 # Speed is average chance of transforming a site per cell cycle.
scar.probabilities <- read.csv("./Data/scar_Probs.csv",
                               stringsAsFactors = F)[, 1:2] # NB In the scar creation
# section there is currently a line that ensures uniqueness of all scars created.
scar.probabilities <- scar.probabilities[1:99, ] # Use this to select only a few scars with lower numbers
scar.probabilities$p <- 1/nrow(scar.probabilities) # Use this to set all probabilities equal
# scar.probabilities$p <- 
#   c(0.5, rep(0.5/(nrow(scar.probabilities) - 1), nrow(scar.probabilities) - 1))
# Use this to create one high-p scar and lots of unique ones.
scar.probabilities <- scar.probabilities[order(-scar.probabilities$p), ]
scar.probabilities$Scar <- 1:nrow(scar.probabilities)
scar.probabilities$Cum.p <- cumsum(scar.probabilities$p)
# colors.hm <- rev(colorRampPalette(c("Red", "Yellow", "Blue"))(100))
cell.types <- data.frame(Cell.type = c("Medium"),
                         Detection.rate = 0.3,
                         Abundance = 1)
# cell.types <- data.frame(Cell.type = c("High", "Medium", "Low"),
#                          Detection.rate = c(0.7, 0.3, 0.15),
#                          Abundance = c(0.15, 0.25, 0.6))
# cell.types$Cum.abundance <- cumsum(cell.types$Abundance)
# integration.types <- data.frame(Integration = c("Strong", "Weak"),
# Detection.rate = c(1, 0.05),
# Abundance = c(0.85, 0.15))
integration.types <- data.frame(Integration = c("Strong"),
                                Detection.rate = c(1),
                                Abundance = c(1))

# Create tree ####
# Parameters: generations
cells <- data.frame(Cell = 1:(2^generations - 1))
cells$Observed <- !(cells$Cell < 2^(generations - 1))
cells$Cell.bin <- sapply(cells$Cell,
                         function(x) paste(binary(x), collapse = ""))
cells$Parent.bin <- sapply(cells$Cell.bin,
                           function(x) substr(x, 1, nchar(x) - 1))
cells$Parent <- sapply(cells$Parent.bin,
                       function(x) decimal(x))
cells$phylo.number <- 
  ifelse(cells$Cell < 2^(generations - 1), 
         2^(generations - 1) + cells$Cell,
         cells$Cell - 2^(generations - 1) + 1)
cells$parent.phylo <-
  ifelse(cells$Parent < 2^(generations - 1), 
         2^(generations - 1) + cells$Parent,
         cells$Parent - 2^(generations - 1) + 1) 
# Note that this incorrectly calculates the first tip to be the parent of the root
cells$parent.phylo[1] <- NA

edges <- as.matrix(cells[cells$Parent > 0, c("parent.phylo", "phylo.number")])
leaves <- cells$Cell[cells$Cell > 2^(generations - 1) - 1]
tree<-list(
  edge=edges,
  tip.label = as.character(cells$phylo.number[cells$Observed]),
  node.label = as.character(cells$phylo.number[!cells$Observed]),
  Nnode= sum(!cells$Observed),
  edge.length = rep(1, nrow(edges)))
class(tree)<-"phylo"

# Simulate scar acquisitions ####
# Parameters: integration sites, scar probabilities, scarring speed
scar.cells <- data.frame(matrix(0, nrow = 2^(generations) - 1, ncol = sites + 3))
colnames(scar.cells) <- c("Cell", "Parent", "Generation", 
                          paste(rep("Site", sites), 1:sites, sep = "."))
scar.cells$Cell <- cells$Cell
scar.cells$Parent <- cells$Parent
scar.cells$Generation <- rep(1:generations, 2^(0:(generations - 1)))

set.seed(2)

cells$Scar.acquisition <- ""
all.scars.created <- integer()
for(cell.in.question in 1:nrow(scar.cells)){
  # Per cell
  print(paste("Cell", cell.in.question))
  # ptm <- proc.time()
  cell.scarring <- data.frame(Site = 1:sites,
                              Scar = as.numeric(scar.cells[cell.in.question, -(1:3)]))
  sites.empty <- sum(!cell.scarring$Scar)
  number.scars.acquired <- sum(runif(sites.empty) > 1 - scar.speed)
  # scar.assign.time <- proc.time() - ptm
  
  # ptm <- proc.time()
  if(number.scars.acquired > 0){
    scar.rolls <- runif(number.scars.acquired)
    scars <- 
      sapply(scar.rolls, 
             function(x) scar.probabilities$Scar[scar.probabilities$Cum.p == 
                                                   min(scar.probabilities$Cum.p[scar.probabilities$Cum.p > x])])
    cells$Scar.acquisition[cell.in.question] <- paste(scars, collapse = ",")
    scar.indices <- which(!cell.scarring$Scar)[1:number.scars.acquired]
    cell.scarring$Scar[scar.indices] <- scars
    all.scars.created <- c(all.scars.created, scars)
  }
  scar.cells[cell.in.question, -(1:3)] <- cell.scarring$Scar
  # scarring.time <- proc.time() - ptm
  
  # ptm <- proc.time()
  if(cell.in.question < nrow(scar.cells)){
    # new implementation
    cell.number <- scar.cells$Cell[cell.in.question]
    cell.generation <- scar.cells$Generation[cell.in.question]
    remaining.generations <- generations - cell.generation
    # remaining.generations <- 4
    cell.progeny <- vector()
    for(g in 1:remaining.generations){
      # progeny <- cell.number * 2^g + 0:(2^g-1)
      cell.progeny <- c(cell.progeny, cell.number * 2^g + 0:(2^g-1))
    }
    scar.cells[scar.cells$Cell %in% cell.progeny, -(1:3)] <-
      scar.cells[cell.in.question, -(1:3)]
  }
  # Below line ensures uniqueness of all scars created by removing all scars,
  # once created, from the pool of scars to choose from.
  # scar.probabilities <- scar.probabilities[!(scar.probabilities$Scar %in% scars), ]
  # scar.propagation.time.2 <- proc.time() - ptm

    #
    # for(c in (cell.in.question + 1):nrow(scar.cells)){
    #   current.parent <- scar.cells$Parent[c]
    #   scar.cells[c, -(1:2)] <- scar.cells[current.parent, -(1:2)]
    # }
  # }
  # scar.propagation.time <- proc.time() - ptm
}  

# Determine scar complexity and percentage wt over generations
wt.dynamics <- data.frame(Generation = 0:generations,
                          WT.perc = 100,
                          Scars = 0)
for(g in 1:generations){
  g.cells <- scar.cells[scar.cells$Generation == g, ]
  all.scars.seen <- melt(g.cells[, -(1:3)], id.vars = NULL)
  wt.dynamics$WT.perc[wt.dynamics$Generation == g] <-
    100 * sum(all.scars.seen$value == 0)/nrow(all.scars.seen)
  wt.dynamics$Scars[wt.dynamics$Generation == g] <-
    length(unique(all.scars.seen$value[all.scars.seen$value > 0]))
}
time.step <- 0.25
wt.dynamics$Time <- wt.dynamics$Generation * time.step
extend.time.to <- 2
if(cell.generation < extend.time.to/time.step){
  total.generations <- extend.time.to/time.step
  wt.dynamics.extended <- 
    data.frame(Generation = (max(wt.dynamics$Generation) + 1):total.generations,
               WT.perc = wt.dynamics$WT.perc[generations + 1],
               Scars = wt.dynamics$Scars[generations + 1])
  wt.dynamics.extended$Time <- wt.dynamics.extended$Generation * time.step
  wt.dynamics <- rbind(wt.dynamics, wt.dynamics.extended)
}
extend.time.to <- max(wt.dynamics$Time)

ggplot(wt.dynamics) +
  geom_point(aes(x = Time, y = WT.perc)) +
  # scale_x_continuous(breaks = 2* 0:ceiling(generations/2)) +
  scale_y_continuous(limits = c(0, 100))
wt.rates <- read.csv("./Data/Dynamics_wt_rates.csv",
stringsAsFactors = F)
fit.function <- function(x,l, a) { a * exp(l*x) + 100 - a}
fit.all <- nls( Percentage ~ fit.function(Timepoint,l, a),
                data = wt.rates, 
                start = list(l=-0.3, a = 90))
fit.all.parameters <- as.list(summary(fit.all)[10][[1]][,1])

# postscript("./Images/Simulations/Tree_B_Simulated_dynamics_with_observed_fit_figure_version.eps",
#            width = 2.5, height = 1.75)
ggplot(wt.dynamics) +
  geom_point(aes(x = Time, y = WT.perc)) +
  stat_function(fun = fit.function, args=fit.all.parameters, color="black",
                xlim = c(0.1, extend.time.to)) +
  scale_x_continuous(limits = c(0, extend.time.to))+
  scale_y_continuous(limits = c(0, 100)) +
  labs(y = "Percentage wildtype", x = "Time (h)") +
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 6))
# dev.off()

# Calculate scars generated per division
scars.created <- wt.dynamics
scars.created$Scars.lastgen <- c(0, scars.created$Scars[-nrow(scars.created)])
scars.created$Scars.this.gen <- scars.created$Scars - scars.created$Scars.lastgen
# postscript("./Images/Scars_created_per_division_figure_version.eps",
#     width = 2.5, height = 1.75)
ggplot(scars.created[scars.created$Time > 0, ]) +
  geom_point(aes(x = Time, y = Scars.this.gen), size = 1) +
  # scale_y_log10(limits = c(1, NA),
  #               breaks = c(1, 10, 100, 1000)) +
  # geom_smooth(aes(x = Time, y = Scars.this.gen), se = F, color = "black",
  #             method = "glm") +
  labs(x = "Time (h)",
       y = "New scars created") +
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 8))
# dev.off()


# Plot developmental and scar tree ####
# Developmental tree
edges.scars <- merge(edges, cells[, c("phylo.number", "Scar.acquisition")])
edges.scars <- edges.scars[, c("parent.phylo", "phylo.number", "Scar.acquisition")]
edges.scars <- edges.scars[order(edges.scars$parent.phylo, edges.scars$phylo.number), ]



# write.csv(edges.scars,
#           "./Data/Simulations/tree_A_dev_tree.csv", quote = F, row.names = F)

# pdf("Images/Simulations/Input_scartree_C.pdf",
#            width = 24, height = 10)
# plot(tree, show.node.label = F, show.tip.label = F, edge.width = 3, no.margin = T)
# edgelabels(edges.scars$Scar.acquisition, frame = "none", adj = c(0.5, -0.2))
# title(main = cells$Scar.acquisition[is.na(cells$parent.phylo)])
# dev.off()

dev.tree.edges <- cells[, c("parent.phylo", "phylo.number", "Scar.acquisition")]
colnames(dev.tree.edges)[1:2] <- c("Parent", "Child")
dev.tree.edges$Parent[is.na(dev.tree.edges$Parent)] <- 0

dev.tree.edges$fill <- "black"
dev.tree.edges$size <- 1
dev.tree.edges$Cell.type <- NA
dev_tree <- generate_tree(dev.tree.edges)
dev_tree_wg <- collapsibleTree(dev_tree, root = dev_tree$scar, collapsed = F,
                               fontSize = 8, width = 300, height = 200,
                               fill = "fill", nodeSize = "size")
dev_tree_wg
# htmlwidgets::saveWidget(dev_tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_C2_dev_tree.html")


# write.table(dev.tree.edges,
#           "./Data/Simulations/tree_B2_dev_tree.txt", quote = F, row.names = F,
#           sep = " ")

# Scar tree
# Start with edges.scars. First collapse all nodes and tips that did not get a
# scar. Then collapse all (resulting) singles. Finally, plot the tree.
edges.scars.collapse <- dev.tree.edges # edges.scars
# Loop over dataframe to find all nodes and tips that did not get a scar
index <- 1
while(index <= nrow(edges.scars.collapse)){
  if(edges.scars.collapse$Scar.acquisition[index] == ""){
    cell.wo.scar.name <- edges.scars.collapse$Child[index]
    upstream.name <- edges.scars.collapse$Parent[index]
    edges.scars.collapse$Parent[edges.scars.collapse$Parent == 
                                        cell.wo.scar.name] <-
      upstream.name
    edges.scars.collapse <- edges.scars.collapse[-index, ]
  }else{index <- index + 1}
}

# Collapse singles
edges.scars.collapse.2 <- edges.scars.collapse
# Loop over dataframe to find all singles.
index <- 1
while(index <= nrow(edges.scars.collapse.2)){
  non.singles <-
    c(unique(edges.scars.collapse.2$Parent[duplicated(edges.scars.collapse.2$Parent)]),
      0)
  if(!(edges.scars.collapse.2$Parent[index] %in% non.singles)){
    # Identify the single, its child, and the edge where the single is the
    # child.
    single.name <- edges.scars.collapse.2$Parent[index]
    child.name <- edges.scars.collapse.2$Child[index]
    single.index <- which(edges.scars.collapse.2$Child == single.name)

    # Collapse the single and its child: replace the name of the single with
    # the name of its child; add the child scar to the single's scars; remove
    # the edge to the child.
    edges.scars.collapse.2$Child[single.index] <- child.name
    edges.scars.collapse.2$Scar.acquisition[single.index] <-
      paste(edges.scars.collapse.2$Scar.acquisition[single.index],
            edges.scars.collapse.2$Scar.acquisition[index], sep = ",")
    edges.scars.collapse.2 <- edges.scars.collapse.2[-index, ]
    
  }else{
    index <- index + 1
  }
}

# tip.converter <- data.frame(Old.name = setdiff(edges.scars.collapse.2$Child, edges.scars.collapse.2$Parent))
# tip.converter$New.name <- 1:nrow(tip.converter)
# node.converter <- data.frame(Old.name = sort(unique(edges.scars.collapse.2$Parent)))
# node.converter$New.name <- (nrow(tip.converter) + 1):(nrow(tip.converter) + nrow(node.converter))
# converter <- rbind(tip.converter, node.converter)
# 
# edges.scars.collapse.2 <- merge(edges.scars.collapse.2, converter,
#                                 by.x = "parent.phylo", by.y = "Old.name")
# colnames(edges.scars.collapse.2)[4] <- "Parent"
# edges.scars.collapse.2 <- merge(edges.scars.collapse.2, converter,
#                                 by.x = "phylo.number", by.y = "Old.name")
# colnames(edges.scars.collapse.2)[5] <- "Child"
# edges.scars.collapse.2 <- 
#   edges.scars.collapse.2[order(edges.scars.collapse.2$Parent,
#                                edges.scars.collapse.2$Child), ]
# write.table(edges.scars.collapse.2[, c("Parent", "Child", "Scar.acquisition")],
#           "./Data/Simulations/tree_B2_scar_tree.csv", quote = F, row.names = F,
#           sep = " ")

# scar.tree <- list(
#   edge = as.matrix(edges.scars.collapse.2[, c("Parent", "Child")]),
#   tip.label = as.character(tip.converter$New.name),
#   root.edge = 1,
#   Nnode = length(unique(edges.scars.collapse.2$V1)),
#   edge.length = rep(1, nrow(edges.scars.collapse.2))
# )
# class(scar.tree) <- "phylo"
# checkValidPhylo(scar.tree)
# pdf("Images/Simulations/Collapsed_scartree_C.pdf",
#            width = 24, height = 14)
# plot(scar.tree, edge.width = 3, no.margin = T, root.edge = T)
# edgelabels(edges.scars.collapse.2$Scar.acquisition, frame = "none", adj = c(0.5, -0.2))
# dev.off()

# Expand cells after scarring ####
# Make [expansion] copies of every leaf cell of the scarred generations.
# Naming convention for cells is [xx]_[yy], with [xx] the number of the
# leaf cell of the scarred generations, and [yy] the identifier of the
# cell within the expansion of [xx]. Note that cell numbers [xx] start with
# 1 and end with 2^(generations - 1), in contrast with the cell numbering
# in [cells] and [scar.cells], but in accordance with phylogenetic methods.
scar.cells.final <-
  data.frame(matrix(ncol = 1 + sites,
                    nrow = 2^(generations - 1) * expansion))
colnames(scar.cells.final) <-
  c("Cell", paste("Site", 1:sites, sep = "."))
scar.cells.final$Cell <-
  paste(rep(1:2^(generations - 1), each = expansion),
        rep(1:expansion, times = 2^(generations - 1)), sep = "_")
scar.cells.final[, -1] <-
  scar.cells[rep(2^(generations - 1):(2^generations - 1), each = expansion),
             -(1:3)]

# Create tree for expanded cells
# The cells that have to be added, together with the cells that they descend
# from.
final.cell.names <- 
  data.frame(Cell = scar.cells.final$Cell,
             Parent = sapply(scar.cells.final$Cell,
                             function(x) unlist(strsplit(x, "_"))[1]))
final.cell.names$Cell  <- as.character(final.cell.names$Cell)
final.cell.names$Parent <- as.character(final.cell.names$Parent)

# We will add the same tree on each tip of the original tree, a tree with
# depth [generations.post] that is a simple binary expansion. We generate
# the edge list for such a tree and then use that same edge list to add all
# trees to the original tree.
exp.edges <- 
  matrix(
    data = c(
      rep((2^generations.post + 1):(2^(generations.post + 1) - 1), each = 2),
      c((2^generations.post + 2):(2^(generations.post + 1) - 1), 1:2^generations.post)),
    nrow = 2^(generations.post + 1) - 2,
    ncol = 2)

expanded.tree <- tree
for(c in 1:sum(cells$Observed)){
  cell.to.expand <- (cells$phylo.number[cells$Observed])[c]
  exp.tip.labels <- final.cell.names$Cell[final.cell.names$Parent == cell.to.expand]
  
  exp.cell.tree <- 
    list(
      edge = exp.edges,
      tip.label = exp.tip.labels,
      Nnode= expansion - 1,
      edge.length = rep(1, 2 * expansion - 2))
  class(exp.cell.tree)<-"phylo"
  
  expanded.tree <- expand.tip(expanded.tree, exp.cell.tree, cell.to.expand)
}
# png("./Images/Simulations/10gen_5gen_tree.png",
#     width = 1600, height = 800, units = "px", type = "cairo")
# plot(expanded.tree, direction = "downward")
# dev.off()

# Assign cell types and scar types ####
set.seed(1)
scar.cells.final$Cell.type<- 
  sample(cell.types$Cell.type, nrow(scar.cells.final), replace = T, 
         prob = cell.types$Abundance)
scar.cells.final <- merge(scar.cells.final, 
                          cell.types[, c("Cell.type", "Detection.rate")])
integration.sites <- 
  data.frame(Site = paste("Site.", 1:sites, sep = ""),
             Type = sample(integration.types$Integration, sites, replace = T, 
                           prob = integration.types$Abundance))
integration.sites <- merge(integration.sites, integration.types[, c("Integration", "Detection.rate")],
                           by.x = "Type", by.y = "Integration")

# Sample cells and scars ####
# Parameters 
cells.sampled <- 125
# detection.rate <- 0.5
set.seed(2)

# Sample cells and scars
readout.wt <- get.readout(scar.cells.final, cells.sampled, integration.sites,
                          doublet.rate = 0)
cells.in.tree <- readout.wt[readout.wt$Scar != 0, ]

length(unique(readout.wt$Cell)) # 1665 cells, including those with only wt.
length(unique(cells.in.tree$Cell)) # 1203 cells with more than wt.

# Write output ####
# write.csv(cells.in.tree, "./Data/Simulations/Tree_B2_2000cellsout_d005.csv",
# quote = F, row.names = F)
# Write output for PHYLIP
cells.in.tree.phylip <- cells.in.tree
cells.in.tree.phylip$Presence <- 1
phylip.out.array <- acast(cells.in.tree.phylip, Cell ~ Scar, value.var = "Presence")
phylip.out.array[is.na(phylip.out.array)] <- "?"
phylip.out <- data.frame(Cell = rownames(phylip.out.array),
                         Scar.vector = apply(phylip.out.array, 1,
                                             function(x) paste(x, collapse = "")),
                         stringsAsFactors = F)
# x <- as.character(phylip.out$Cell[1])
# y <- paste(x, paste(rep("x", 10 - nchar(x)), collapse = ""), sep = "")
if(sum(grepl(";", phylip.out$Cell)) > 0){
  # Doublet formatting does not work for PHYLIP: too long and ";" is forbidden.
  doublets <- data.frame(Doublet = phylip.out$Cell[grepl(";", phylip.out$Cell)])
  doublets$New.name <- paste("d_", 1:nrow(doublets), sep = "")
  phylip.out$Cell <- 
    mapvalues(phylip.out$Cell, from = doublets$Doublet, to = doublets$New.name)
}
phylip.out$Cell <-
    sapply(as.character(phylip.out$Cell), function(x){
      name.length = nchar(x)
      if(name.length < 10){
        y <- paste(x, paste(rep("x", 10 - nchar(x)), collapse = ""), sep = "")
      }
      return(y)
    })

colnames(phylip.out) <- c(nrow(phylip.out), ncol(phylip.out.array))
phylip.scar.conversion <- data.frame(Scar.order = 1:ncol(phylip.out.array),
                                     Scar = colnames(phylip.out.array))
# write.table(phylip.out, "./Data/Simulations/Tree_B2_2000cellsout_d005_phylip_?.txt",
# sep = " ", row.names = F, quote = F)
# write.csv(phylip.scar.conversion, "./Data/Simulations/Tree_B2_scar_conversion_d005.csv",
#           row.names = F, quote = F)
# write.table(doublets, "./Data/Simulations/Tree_B2_2000cellsout_d005_doublets.csv",
#           row.name = F, quote = F, sep = ",")
