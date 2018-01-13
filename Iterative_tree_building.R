# Description ####
# Iterative tree building based on the comparison of degrees and detection rates
# in the scar graph.
# The test we use is an extension of the naive "the-first-scar-should-have-the-
# highest-degree" scheme one would use if there are no dropouts. The extension
# is that we now look for the scar that would have had the highest degree if
# there wouldn't have been dropouts. Once this scar has been determined, it is
# removed from the dataset, in which we again determine the first scar, and so
# on.
# To find the first scar, we simply try out each scar. Under the assumption that
# a scar is the first created scar, we can calculate its detection rate as the
# amount of times we see it coincide with another scar, divided by the total
# amount of times we see cells with other scars. Given this detection rate and
# the amount of times every other scar is observed, we calculate the expected
# number of connections the first scar should have and compare it to the
# observed number of connections. The underlying distribution for the number
# of connections is the Poisson binomial poisson distribution, so we use this
# to calculate the p-value of the observed value under the assumption of the
# expected value. After calculating this p-value for all scars, we select the
# scar with the highest p-value as the first scar.

# Dependencies ####
source("./Scripts/linnaeus-scripts/scar_helper_functions.R")

# Parameters ####
# Fraction of doublets expected; number of connections has to be higher than
# the expected number of doublets + 2sigma under the assumption that the
# number of doublets is binomially distributed.
doublet.rate <- 0.1 # Default is 0.1, set to 0 to turn off.
# The minimum detection rate for a scar to be considered as top scar.
min.detection.rate <- 0.01 # Default value is 0.1
# Minimum cell number ratio between branches.
branch.size.ratio <- 0 # Default 0.25, set to 0 to turn off
# Maximum scar probability to include scar in tree building
max.scar.p <- 0.001
# Maximum number of embryos a scar can be present in to include in tree building
max.larvae <- 1

# For testing purposes: how many scars to include in tree building (takes the
# most frequent scars, set to NA to include all)
number.scars <- 10

# Load data ####
# mRNA
tsne.coord.in <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells.csv")
# Count total number of cells present even without scars
# For Z2
tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L2", c("Barcode", "Cluster")]
# For Z4
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L4", c("Barcode", "Cluster")]
# For Z5
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L5", c("Barcode", "Cluster")]
# For A5
# N <- sum(grepl("B5|H5|P5", tsne.coord$Cell))
# For (simulated) tree B
# N <- 3000
N <- nrow(tsne.coord)

# Scars
scar.input <- 
  # read.csv("./Data/Simulations/Tree_B_3k_cells_3celltypes_2sites.csv")
  # read.csv("./Data/2017_10X_7/A5_used_scars_2.csv", stringsAsFactors = F)
  read.csv("./Data/2017_10X_2/Z2_scars_compared.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_10_CR/Z4_scars_compared.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_10_CR/Z5_scars_compared.csv", stringsAsFactors = F)
scar.input <- merge(scar.input[, c("Barcode", "Scar", "Presence", "p")], 
                    tsne.coord)
colnames(scar.input)[which(colnames(scar.input) == "Cluster")] <-
  "Cell.type"
colnames(scar.input)[which(colnames(scar.input) == "Barcode")] <-
  "Cell"

if(!("Cell.type" %in% names(scar.input))){
  scar.input$Cell.type <- "Type.O.Negative"
}
if("p" %in% names(scar.input)){
  cells.in.tree <- scar.input[scar.input$p <= max.scar.p, 
                              c("Cell", "Scar", "Cell.type")]
  if("Presence" %in% names(scar.input)){
    cells.in.tree <- 
      scar.input[scar.input$p <= max.scar.p & scar.input$Presence <= max.larvae, 
                 c("Cell", "Scar", "Cell.type")]
  }
}else{
  cells.in.tree <- scar.input
}
cells.in.tree <- cells.in.tree[!duplicated(cells.in.tree), ]

scar.freqs <- data.frame(table(cells.in.tree$Scar))
colnames(scar.freqs)[1] <- "Scar"
scar.freqs <- scar.freqs[order(-scar.freqs$Freq), ]
set.seed(1)
if(is.na(number.scars)){
  include.scars <- scar.freqs$Scar
}else{
  include.scars <- scar.freqs$Scar[1:number.scars]
}
cells.in.tree <- cells.in.tree[cells.in.tree$Scar %in% include.scars, ]

# Filter out low-frequency scar connections ####
# Count how often every scar-scar connection is seen
scar.connections <- connections.for.graph(cells.in.tree)
only.once.connections <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
colnames(only.once.connections) <- c("Scar.A", "Scar.B")
only.once.connections <- 
  merge(only.once.connections, 
        scar.connections)

# Calculate how many doublets we'd expect given a doublet rate, and (under
# assumption that the number of doublets is binomially distributed) calculate
# an expected value + 2*sigma threshold
only.once.connections$AB.doublets <- 
  2 * doublet.rate * only.once.connections$x_A * only.once.connections$x_B/N
only.once.connections$AB.doublet.rate <- only.once.connections$AB.doublets/N
only.once.connections$AB.doublet.sd <- 
  sqrt(N * only.once.connections$AB.doublet.rate * 
         (1 - only.once.connections$AB.doublet.rate))
only.once.connections$AB.doublet.threshold <- 
  only.once.connections$AB.doublets + 2 * only.once.connections$AB.doublet.sd

ooc.cutoff <- 
  only.once.connections[only.once.connections$x_AB > 
                          only.once.connections$AB.doublet.threshold, ]

ooc.cutoff.graph <-
  graph_from_data_frame(ooc.cutoff[, c("Scar.A", "Scar.B")],
                        directed = F, vertices = union(ooc.cutoff$Scar.A, 
                                                       ooc.cutoff$Scar.B))
# pdf("Images/Simulations/Z2_network_10scars_p01_coincr1_detratio01_branchratio0.pdf",
#            width = 20, height = 10)
plot(ooc.cutoff.graph)
# dev.off()

# Identify incorrect connections
incorrect.connections <- 
  only.once.connections[only.once.connections$x_AB > 0 & 
                          only.once.connections$x_AB <= 
                          only.once.connections$AB.doublet.threshold, ]

# Identify and remove cells with incorrect connections
# ic <- 1
inc.cells <- character()
for(ic in 1:nrow(incorrect.connections)){
  inc.scar.A <- incorrect.connections$Scar.A[ic]
  inc.scar.B <- incorrect.connections$Scar.B[ic]
  cells.scar.A <- cells.in.tree$Cell[cells.in.tree$Scar == inc.scar.A]
  cells.scar.B <- cells.in.tree$Cell[cells.in.tree$Scar == inc.scar.B]
  
  inc.cells <- unique(c(inc.cells, intersect(cells.scar.A, cells.scar.B)))
}
cells.in.tree.f <- cells.in.tree[!(cells.in.tree$Cell %in% inc.cells), ]

# Iteration start conditions ####
scar.amount <- length(unique(cells.in.tree.f$Scar))
it.tree.building <- vector("list", scar.amount)
tree.summary <- 
  initialize.branches(cells.in.tree.f, scar.remove = "Root",
                      size.ratio = branch.size.ratio)

# Iterative tree building ####
scar.index <- 1
while(scar.index <= scar.amount){
  print(paste("Iterative tree building, identifying scar",
              scar.index, "of", scar.amount))
  
  # Select the topmost incomplete line in the tree summary to define the scar 
  # graph branch for which to determine the first scar created. Deduce the other 
  # scars removed for this branch and which component corresponds to the branch.
  top.incomplete.edge.index <- min(which(is.na(tree.summary$Node.2)))
  top.incomplete.edge <- tree.summary[top.incomplete.edge.index, ]
  
  # Determine which scars to remove
  starting.scar <- as.character(top.incomplete.edge$Node.1)
  scars.to.remove <- find.scars.to.remove(starting.scar, tree.summary)
  
  # Finding the correct component is not fully straightforward: it may be the
  # second connected component of the third connected component of the seventh
  # connected component, for example, when scars are removed starting from the 
  # top. In tree.summary we only record the component number of the last 
  # branching event, and we record that mostly to ensure we are covering all
  # components. While this means the 'component history' of any connected
  # component at any level can be reconstructed (simply follow the scars back up
  # and read out the component numbers), this is not necessary.
  # To select the correct component, we remove scars in two steps: we first 
  # remove all scars but the bottommost one (the [starting.scar]) and create the
  # graph with those scars removed. As a second step, we select the component
  # that includes the [starting.scar], remove that scar, and find the correct
  # component in the resulting disconnected graph.
  preceding.scars <- setdiff(scars.to.remove, starting.scar)
  current.cs <- cells.in.tree.f[!(cells.in.tree.f$Scar %in% preceding.scars), ]
  first.decomposition <- graph.and.decompose(current.cs)
  for(comp in 1:length(first.decomposition)){
    if(starting.scar %in% V(first.decomposition[[comp]])$name){
      break
    }
  }
  if(starting.scar == "Root"){
    second.decomposition <- 
      first.decomposition
  }else{
    second.decomposition <- 
      decompose(delete_vertices(first.decomposition[[comp]], starting.scar))
  }
  current.component <- top.incomplete.edge$Component
  current.graph <- second.decomposition[[current.component]]
  current.cs.component <- current.cs[current.cs$Scar %in% V(current.graph)$name, ]
  
  # cs <- current.cs.component
  # graph <- current.graph
  
  if(length(V(current.graph)) == 1){
    # Condition for last scar in branch - this will have a graph without any
    # connections, and with one vertex.
    scar.remove <- V(current.graph)$name
    it.tree.element <- list(Scar = scar.remove)
  }else{
    # Test congruence of degree and detection efficiency of all scars, assuming
    # they are the topmost scar of the current data.
    scar.lls <- create.degree.lls(current.cs.component, 
                                  current.graph)
    
    # Calculate the weighted average detection rate of all scars.
    average.det.rates <- 
      ddply(scar.lls, .(Scar),
            function(x) data.frame(Mean.p_A = weighted.mean(x$p_A, x$Total.other)))
    
    # Select only scars whose average detection rate is higher than cutoff
    scar.lls <- merge(scar.lls, average.det.rates)
    scar.lls <- scar.lls[order(-scar.lls$Degree.p,
                               -scar.lls$Degree,
                               -scar.lls$Mean.p_A), ]
    scar.lls.select <- scar.lls[scar.lls$Mean.p_A > min.detection.rate, ]
    scar.lls.unique <- 
      unique(scar.lls[, c("Scar", "Degree", "Scar.count", "Expected.degree", 
                          "Degree.p", "Mean.p_A")])
    scar.lls.select.unique <- 
      scar.lls.unique[scar.lls.unique$Mean.p_A > min.detection.rate, ]
    
    if(nrow(scar.lls.select) > 0){
      scar.remove <- scar.lls.select.unique$Scar[1]
    }else{
      print("No scar above minimum detection rate. Taking best scar under minimum detection rate")
      scar.remove <- scar.lls.unique$Scar[1]
    }
    
    it.tree.element <- list(Scar = scar.remove,
                            LLS = scar.lls,
                            LLS.select = scar.lls.select,
                            LLS.unique = scar.lls.unique,
                            LLS.select.unique = scar.lls.select.unique)
  }
  it.tree.building[[scar.index]] <- it.tree.element
  
  # Add scar to tree.summary as Node.2
  tree.summary$Node.2[top.incomplete.edge.index] <- scar.remove
  
  # Add newfound scar and components to tree.summary.
  remaining.cs <- current.cs.component[current.cs.component$Scar != scar.remove, ]
  if(nrow(remaining.cs) > 0){
    tree.summary.add <- 
      initialize.branches(remaining.cs, scar.remove = scar.remove, 
                          size.ratio = branch.size.ratio)
    tree.summary.add$Depth <- tree.summary$Depth[top.incomplete.edge.index] + 1
    # If current branch was not a main branch, set non-main flags for all
    # consecutive scars as well.
    if(!tree.summary$Main[top.incomplete.edge.index]){
      tree.summary.add$Main <- F
    }
    tree.summary <- rbind(tree.summary, tree.summary.add)
  }
  
  scar.index <- scar.index + 1
}

# Collapse tree ####
# Remove 'singles' (i.e. scars in tree.summary$Node.1 that only occur once)
# from the tree.summary by collapsing them with their successors while changing
# the name to "[scar 1], [scar 2]".
tree.summary.collapse <- tree.summary[tree.summary$Main, ]
tree.summary.old <- tree.summary
# Loop over dataframe to find all singles.
index <- 1
while(index <= nrow(tree.summary.collapse)){
  non.singles <-
    unique(tree.summary.collapse$Node.1[duplicated(tree.summary.collapse$Node.1)])
  if(!(tree.summary.collapse$Node.1[index] %in% non.singles)){
    single.name <- tree.summary.collapse$Node.1[index]
    downstream.name <- tree.summary.collapse$Node.2[index]
    collapsed.name <- paste(single.name, downstream.name, sep = ",")
    tree.summary.collapse$Node.1[tree.summary.collapse$Node.1 %in% 
                                   c(single.name, downstream.name)] <- 
      collapsed.name
    tree.summary.collapse$Node.2[tree.summary.collapse$Node.2 %in% 
                                   c(single.name, downstream.name)] <- 
      collapsed.name
    tree.summary.collapse <- tree.summary.collapse[-index, ]
  }else{index <- index + 1}
  if(nrow(tree.summary.collapse) == 1){break}
}
tree.summary <- tree.summary.collapse

# Create phylogenetic tree ####
# Create a list that includes edges, tips, nodes and possibly labels,
# then turn that into an object of class "phylo", then plot.
tips <- data.frame(Name = setdiff(tree.summary$Node.2, tree.summary$Node.1),
                   stringsAsFactors = F)
tips$Index <- 1:nrow(tips)
nodes <- data.frame(
  Name = setdiff(c(tree.summary$Node.1, tree.summary$Node.2), tips$Name),
  stringsAsFactors = F)
nodes$Index <- NA
nodes$Index[grep("Root", nodes$Name)] <- nrow(tips) + 1
nodes <- nodes[order(nodes$Index), ]
nodes$Index[-1] <- (nrow(tips) + 2):(nrow(tips) + nrow(nodes))
nodes$Index <- as.integer(nodes$Index)
nodestips <- rbind(nodes, tips)
phylo.edges <- merge(tree.summary[, c("Node.1", "Node.2")], nodestips,
                     by.x = "Node.1", by.y = "Name")
colnames(phylo.edges)[3] <- "V1"
phylo.edges <- merge(phylo.edges, nodestips,
                     by.x = "Node.2", by.y = "Name")
colnames(phylo.edges)[4] <- "V2"

nodes.2 <- nodes
nodes.2$Name[grep("Root", nodes.2$Name)] <-
  sub("Root,", "", nodes.2$Name[grep("Root", nodes.2$Name)])

scar.phylo <-
  list(
    edge = as.matrix(phylo.edges[, c("V1", "V2")]),
    tip.label = tips$Name,
    edge.length = rep(1, nrow(phylo.edges)),
    Nnode = nrow(nodes.2),
    node.label = nodes.2$Name,
    root.edge = 1)
class(scar.phylo) <- "phylo"

# Plot tree ####
# pdf("Images/Simulations/Z2_15scars_L1_p001_doubletrate01_detratio01_branchratio0.pdf",
# width = 20, height = 10)
plot(scar.phylo, show.node.label = F, show.tip.label = F, root.edge = T,
     edge.width = 3, no.margin = T, direction = "leftward")
# title(main = sub("Root,", "", nodes$Name[grep("Root", nodes$Name)]))
edgelabels(phylo.edges$Node.2, frame = "none", adj = c(0.5, 0), cex = 2,
           col = "red")
# dev.off()

# Investigate tree building ####
View(tree.summary.old)
View(tree.summary)
# View(it.tree.building[[1]]$LLS)
# View(it.tree.building[[1]]$LLS.select.unique)
View(it.tree.building[[5]]$LLS.unique)
detection.rate.progression <-
  data.frame(Scar = character(),
             Step = integer(),
             Detection.rate = numeric())
for(n.scar in 1:length(it.tree.building)){
  detection.rate.add <- it.tree.building[[n.scar]]$LLS.unique[, c("Scar", "Mean.p_A")]
  if(length(detection.rate.add) >= 1){
    colnames(detection.rate.add)[2] <- "Detection.rate"
    detection.rate.add$Step = n.scar
    
    detection.rate.progression <- rbind(detection.rate.progression,
                                        detection.rate.add)
  }
}
print(ggplot(detection.rate.progression) +
        geom_tile(aes(x = Step, y = Scar, fill = Detection.rate)) +
        scale_fill_gradient(low = "grey", high = "red")
)
# ggplot(detection.rate.progression) +
#   geom_line(aes(x = Step, y = Detection.rate, color = Scar)) +
#   scale_color_manual(values = rep("black", 9))

# Investigations ####
scar.1 <- "515:39M7D2M10D34M"
scar.2 <- "622:48M1D27M"
cells.with.1 <- cells.in.tree.f$Cell[cells.in.tree.f$Scar == scar.1]
cells.with.2 <- cells.in.tree.f$Cell[cells.in.tree.f$Scar == scar.2]
cells.with.12 <- intersect(cells.with.1, cells.with.2)
View(cells.in.tree.f[cells.in.tree.f$Cell %in% cells.with.12, ])
cs.with.1 <- cells.in.tree.f[cells.in.tree.f$Cell %in% cells.with.1, ]
cs.with.2 <- cells.in.tree.f[cells.in.tree.f$Cell %in% cells.with.2, ]
View(cells.in.tree.f[cells.in.tree.f$Cell %in% cells.with.1, ])
View(scars.in.1[scars.in.1$Barcode %in% cells.with.12, ])
View(scars.in.2[scars.in.2$Barcode %in% cells.with.12, ])
