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
source("./Scripts/scar_helper_functions.R")

# Parameters ####
# Minimum ratio between scar coincidence rate and constituent scar occurrence 
# rates - p_AB/(p_Ap_B) for scars A and B.
min.coinc.occurence.ratio <- 0 # Default 1, set to 0 to include all connections.
# Minimal detection rate ratio for scar to be considered as top scar in an 
# iteration.
min.detection.rate.ratio <- 0.1 # Default 0.1, set to 0 to turn off.
# Minimum cell number ratio between branches.
branch.size.ratio <- 0 # Default 0.25, set to 0 to turn off
# Maximum scar probability to include scar in tree building
max.scar.p <- 0.01


# Load data ####
scar.input <- read.csv("./Data/Simulations/Tree_B_3k_cells_3celltypes_2sites.csv")
  # read.csv("./Data/2017_10X_7/A5_used_scars_2.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_2/Z2_used_scars_2.csv", stringsAsFactors = F)
if(!("Cell.type" %in% names(scar.input))){
  scar.input$Cell.type <- "Type.O.Negative"
}
if("p" %in% names(scar.input)){
  cells.in.tree <- scar.input[scar.input$p <= max.scar.p, 
                              c("Cell", "Scar", "Cell.type")]
}else{
  cells.in.tree <- scar.input
}
cells.in.tree <- cells.in.tree[!duplicated(cells.in.tree), ]

# Count total number of cells present even without scars
tsne.coord <- read.csv("./Data/2017_10X_2/X10_final_all_tsne_Seurat_Bo_ID.csv")
N <- sum(grepl("Z2", tsne.coord$Cell)) 

scar.freqs <- data.frame(table(cells.in.tree$Scar))
colnames(scar.freqs)[1] <- "Scar"
scar.freqs <- scar.freqs[order(-scar.freqs$Freq), ]
include.scars <- scar.freqs$Scar #[1:10]
cells.in.tree <- cells.in.tree[cells.in.tree$Scar %in% include.scars, ]

# Filter out low-frequency scar connections ####
# Count how often every scar-scar connection is seen
scar.connections <- connections.for.graph(cells.in.tree)
only.once.connections <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
colnames(only.once.connections) <- c("Scar.A", "Scar.B")
only.once.connections <- 
  merge(only.once.connections, 
        scar.connections)

# Compare the connection frequency to the constituent scars frequencies
# and remove connections that aren't seen often enough.
only.once.connections$AB.ratio.obs <- only.once.connections$x_AB/N
only.once.connections$AB.ratio.exp <- 
  only.once.connections$x_A * only.once.connections$x_B/N^2
only.once.connections$AB.obs.exp <- 
  only.once.connections$AB.ratio.obs/only.once.connections$AB.ratio.exp
ooc.cutoff <- only.once.connections[only.once.connections$AB.obs.exp >= 
                                      min.coinc.occurence.ratio, ]
ooc.cutoff.graph <-
  graph_from_data_frame(ooc.cutoff[, c("Scar.A", "Scar.B")],
                        directed = F, vertices = union(ooc.cutoff$Scar.A, 
                                                       ooc.cutoff$Scar.B))
# plot(ooc.cutoff.graph)

# Identify incorrect connections
incorrect.connections <- 
  only.once.connections[only.once.connections$AB.obs.exp > 0 & 
                          only.once.connections$AB.obs.exp < 
                          min.coinc.occurence.ratio, ]

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
      first.decomposition #[[top.incomplete.edge$Component]]
  }else{
    second.decomposition <- 
      decompose(delete_vertices(first.decomposition[[comp]], starting.scar))
  }
  current.component <- top.incomplete.edge$Component
  current.graph <- second.decomposition[[current.component]]
  current.cs.component <- current.cs[current.cs$Scar %in% V(current.graph)$name, ]
  
  # tree.summary$Size[top.incomplete.edge.index] <- 
  #   length(unique(current.cs.component$Cell))
  
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
    
    # Only test for scars that have a high average detection rate; set to
    # on average 1/10 of the maximum scar detection rate per cell type. To
    # determine this average detection rate ratio, first compute the max
    # detection rate per cell type.
    max.p.celltypes <- aggregate(scar.lls$p_A,
                                 by = list(Cell.type = scar.lls$Cell.type),
                                 max)
    colnames(max.p.celltypes)[2] <- "Max.p"
    # Then list all combinations of scars and cell types, and add their
    # detection rates and the maximum detection rate per cell type.
    possible.top.scars <- expand.grid(unique(scar.lls$Scar), unique(scar.lls$Cell.type))
    colnames(possible.top.scars) <- c("Scar", "Cell.type")
    possible.top.scars <- 
      merge(merge(possible.top.scars, scar.lls, all.x = T), max.p.celltypes)
    possible.top.scars[is.na(possible.top.scars)] <- 0
    # Compute the ratio of the detection rate of a scar in a cell type to the
    # maximum detection rate in that cell type.
    possible.top.scars$p_A.to.max <- 
      possible.top.scars$p_A/possible.top.scars$Max.p
    # Remove scar/cell type combinations that have no maximum detection rate
    # (i.e. none of the scars under consideration are in this cell type)
    possible.top.scars <- possible.top.scars[complete.cases(possible.top.scars), ]
    # Calculate the mean detection rate to max detection rate ratio per scar
    top.scars <- aggregate(possible.top.scars$p_A.to.max,
                           by = list(Scar = possible.top.scars$Scar),
                           mean)
    colnames(top.scars)[2] <- "Mean.to.max.rate"
    # Determine valid top scar candidates and subset the likelihood computations
    valid.top.scars <- 
      top.scars$Scar[top.scars$Mean.to.max.rate >= min.detection.rate.ratio]
    scar.lls <- scar.lls[scar.lls$Scar %in% valid.top.scars, ]

    # possible.top.scars <- merge(scar.lls, max.p.celltypes)
    # possible.top.scars$p_A.to.max <- 
    #   possible.top.scars$p_A/possible.top.scars$Max.p
    # top.scars <- aggregate(possible.top.scars$p_A.to.max,
    #                        by = list(Scar = possible.top.scars$Scar),
    #                        mean)
    # 
    # possible.top.scars <-
    #   possible.top.scars[possible.top.scars$p_A >= possible.top.scars$Max.p/10, ]
    # scar.lls <- scar.lls[scar.lls$Scar %in% unique(possible.top.scars$Scar), ]
    
    # Old filtering version (worked with simulations but real data has more
    # variation)
    # possible.top.scars <- 
    #   scar.lls$Scar[scar.lls$p_A >= max(scar.lls$p_A)/10]
    # scar.lls <- scar.lls[scar.lls$Scar %in% unique(possible.top.scars$Scar), ]
    
    scar.lls.unique <- 
      unique(scar.lls[, c("Scar", "Degree", "Scar.count", "Expected.degree", "Degree.p")])
    # scar.lls.unique <- scar.lls.unique[scar.lls.unique$Scar %in% possible.top.scars, ]
        
    scar.remove <- scar.lls.unique$Scar[1]
    it.tree.element <- list(Scar = scar.remove,
                            LLS = scar.lls)
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
# pdf("Images/Simulations/Tree_B_3types_2sites_coincr0_detratio01_branchratio0.pdf",
#            width = 12, height = 7)
plot(scar.phylo, show.node.label = F, show.tip.label = F, root.edge = T,
     edge.width = 3, no.margin = T, direction = "leftward")
# title(main = sub("Root,", "", nodes$Name[grep("Root", nodes$Name)]))
edgelabels(phylo.edges$Node.2, frame = "none", adj = c(0.5, -0.5))
# dev.off()

# Investigate tree building ####
View(tree.summary.old)
# View(it.tree.building[[1]]$LLS)
