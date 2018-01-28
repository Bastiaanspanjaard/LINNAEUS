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
# Fraction of doublets expected.
doublet.rate <- 0.1 # Default is 0.1, set to 0 to turn off.
# The minimum detection rate for a scar to be considered as top scar.
min.detection.rate <- 0.01 # Default value is 0.01
# Minimum cell number ratio between branches.
branch.size.ratio <- 0.25 # Default 0.25, set to 0 to turn off
# Maximum scar probability to include scar in tree building
max.scar.p <- 0.001
# Maximum number of embryos a scar can be present in to include in tree building
max.larvae <- 1

# For testing purposes: how many scars to include in tree building (takes the
# most frequent scars, set to NA to include all)
number.scars <- NA

# Load data ####
print("Loading data")
# mRNA
# tsne.coord.in <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells_2.csv")
# Count total number of cells present even without scars
# For Z2
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library %in% c("L21", "L22"), 
#                             c("Cell", "Cluster", "Cell.type")]
# For Z4
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L4", c("Cell", "Cluster", "Cell.type")]
# For Z5
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L5", c("Barcode", "Cluster")]
# For A5
# N <- sum(grepl("B5|H5|P5", tsne.coord$Cell))
# For (simulated) tree B
N <- 3000 #125 #
# N <- nrow(tsne.coord)

# Scars
scar.input <- 
  # read.csv("./Data/Simulations/Tree_C2_100cellsout_detection03.csv")
  # read.csv("./Data/Simulations/Tree_B2_2000cellsout.csv")
  # read.csv("./Data/Simulations/Tree_B2_2000cellsout_d005.csv")
  # read.csv("./Data/Simulations/Tree_B2_2000cellsout_d0_wweakint.csv")
  read.csv("./Data/Simulations/Tree_B2_2000cellsout_d005_wweakint.csv")
  # read.csv("./Data/Simulations/Tree_B3_2000cellsout_d0_wweakint.csv")
  # read.csv("./Data/2017_10X_7/A5_used_scars_2.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_2/Z2_scars_compared.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_10_CR/Z4_scars_compared.csv", stringsAsFactors = F)
# scar.input$Cell <- paste("L4", scar.input$Barcode, sep = "_")
  # read.csv("./Data/2017_10X_10_CR/Z5_scars_compared.csv", stringsAsFactors = F)
# scar.input <- merge(scar.input[, c("Cell", "Scar", "Presence", "p")],
#                     tsne.coord)
# colnames(scar.input)[which(colnames(scar.input) == "Cluster")] <-
#   "Cell.type"
# colnames(scar.input)[which(colnames(scar.input) == "Barcode")] <-
#   "Cell"

if(!("Cell.type" %in% names(scar.input))){
  scar.input$Cell.type <- "Type.O.Negative"
}
if("p" %in% names(scar.input)){
  cells.in.tree.pre.f <- scar.input[scar.input$p <= max.scar.p, 
                              c("Cell", "Scar", "Cell.type")]
  if("Presence" %in% names(scar.input)){
    cells.in.tree.pre.f <- 
      scar.input[scar.input$p <= max.scar.p & scar.input$Presence <= max.larvae, 
                 c("Cell", "Scar", "Cell.type")]
  }
}else{
  cells.in.tree.pre.f <- scar.input
}
cells.in.tree.pre.f <- cells.in.tree.pre.f[!duplicated(cells.in.tree.pre.f), ]

# cells.in.tree <- cells.in.tree[!grepl(";", cells.in.tree$Cell), ]

scar.freqs <- data.frame(table(cells.in.tree.pre.f$Scar))
colnames(scar.freqs)[1] <- "Scar"
scar.freqs <- scar.freqs[order(-scar.freqs$Freq), ]
set.seed(1)
if(is.na(number.scars)){
  include.scars <- scar.freqs$Scar[scar.freqs$Freq > 1]
}else{
  include.scars <- scar.freqs$Scar[1:number.scars]
}
cells.in.tree.pre.f <- cells.in.tree.pre.f[cells.in.tree.pre.f$Scar %in% include.scars, ]


# Start tree building ####
# Try to build a tree and detect any weak scars until a tree has been built
# without finding weak scars.
weak.scars <- character()
repeat{
  print(paste("Tree building without scars", paste(weak.scars, collapse = ", ")))
  weak.scar.found <- F
  # Remove weak scars
  # weak.scars <- NULL
  # c(934, 206) weak scars for B2 with doublets.
  cells.in.tree <- 
    cells.in.tree.pre.f[!(cells.in.tree.pre.f$Scar %in% weak.scars), ]
  
  
  # Filter out low-frequency scar connections ####
  print("Filtering doublets")
  dfilter.start <- Sys.time()
  
  # dataset.graph <- graph.and.decompose(cells.in.tree)
  # dataset.degrees <-
  #   data.frame(lapply(dataset.graph, function(x) degree(x, mode = "all", loops = F)))
  # colnames(dataset.degrees) <- "Degree"
  # dataset.degrees$Scar <- as.character(rownames(dataset.degrees))
  # dataset.counts <- data.frame(table(cells.in.tree$Scar))
  # colnames(dataset.counts) <- c("Scar", "Count")
  # dataset.dc <- merge(dataset.counts, dataset.degrees)
  # dataset.dc$scardeg.ratio <- dataset.dc$Count/dataset.dc$Degree
  # ggplot(dataset.dc) +
  #   geom_histogram(aes(x = scardeg.ratio))
  # dataset.dc$degscar.ratio <- dataset.dc$Degree/dataset.dc$Count
  # ggplot(dataset.dc) +
  #   geom_histogram(aes(x = log(degscar.ratio)))
  # dataset.dc$log.dsr <- log(dataset.dc$degscar.ratio)
  # 
  # 
  # cd.lm <- lm(Degree ~ Count + 0, dataset.dc)
  # dataset.dc$cd.lm.ratio <- cd.lm$residuals/dataset.dc$Degree
  # dataset.dc$cd.lm <- cd.lm$residuals
  # 
  # dataset.dc$Exp.d.max <- dataset.dc$Count/20
  # dataset.dc$exp.d.max.diff <- dataset.dc$Exp.d.max - dataset.dc$Degree
  
  # 
  # scars.high.ratio <- dataset.dc$Scar[dataset.dc$scardeg.ratio > 3]
  # 
  # cells.in.tree <- cells.in.tree[cells.in.tree$Scar %in% scars.high.ratio, ]
  
  # dc.lm <- lm(Count ~ Degree + 0, dataset.dc)
  # dataset.dc$LM.res <- dc.lm$residuals/dataset.dc$Count
  # dataset.dc$Scaled.res <- dataset.dc$LM.res/dataset.dc$Count
  
  # Count how often every scar-scar connection is seen
  scar.connections <- connections.for.graph(cells.in.tree)
  only.once.connections <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
  colnames(only.once.connections) <- c("Scar.A", "Scar.B")
  only.once.connections <- 
    merge(only.once.connections, 
          scar.connections)
  
  # NEW
  only.once.connections <- only.once.connections[only.once.connections$x_AB > 0, ]
  
  # Investigate difference between dataset with doublets and without
  # cells.in.tree.no.d <- cells.in.tree[cells.in.tree$Cell.type != "Doublet", ]
  # scar.connections.no.d <- connections.for.graph(cells.in.tree.no.d)
  # only.once.connections.no.d <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
  # colnames(only.once.connections.no.d) <- c("Scar.A", "Scar.B")
  # only.once.connections.no.d <- 
  #   merge(only.once.connections.no.d, 
  #         scar.connections.no.d)
  # colnames(only.once.connections.no.d)[3:5] <- c("xnd_A", "xnd_B", "xnd_AB")
  
  
  # Calculate how many doublets we'd expect given a general doublet rate, and 
  # calculate the doublet rate for every scar connection.
  only.once.connections$AB.doublets <- 
    2 * doublet.rate * only.once.connections$x_A * only.once.connections$x_B/N
  only.once.connections$AB.doublet.rate <- only.once.connections$AB.doublets/N
  
  # Calculate the probability that the connections we see can all be explained by
  # doublets.
  only.once.connections$Doublet.p <-
    apply(only.once.connections[, c("x_AB", "AB.doublet.rate")], 1,
          function(x){
            x_d <- as.integer(x[1])
            p_d <- as.numeric(x[2])
            binom.test(x_d, N, p_d, alternative = "greater")$p.value
          }
    )
  only.once.connections$Doublet.padj <-
    p.adjust(only.once.connections$Doublet.p, "fdr")
  
  ooc.cutoff <- 
    only.once.connections[only.once.connections$Doublet.padj > 0.01, ]
  
  # Identify incorrect connections
  incorrect.connections <- 
    only.once.connections[only.once.connections$Doublet.padj >= 0.01, ]
  
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
  dfilter.end <- Sys.time()
  dfilter.time <- dfilter.end - dfilter.start
  print(dfilter.time)
  
  # Iteration start conditions ####
  iterative.sc.start <- Sys.time()
  scar.amount <- length(unique(cells.in.tree.f$Scar))
  it.tree.building <- vector("list", scar.amount)
  tree.summary <- 
    initialize.branches(cells.in.tree.f, scar.remove = "Root",
                        size.ratio = branch.size.ratio)
  iterative.sc.end <- Sys.time()
  iterative.sc.time <- iterative.sc.end - iterative.sc.start
  print(iterative.sc.time)
  
  # Iterative tree building ####
  scar.index <- 1
  iterative.start <- Sys.time()
  while(scar.index <= scar.amount & !weak.scar.found){
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
    
    cs <- current.cs.component
    graph <- current.graph
    
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
        difficult.scars <-
          scar.lls$Scar[scar.lls$Degree.p >
                          scar.lls$Degree.p[scar.lls$Scar == scar.remove]]
        if(length(difficult.scars) > 0){
          print(paste("Weak scar", difficult.scars[1], "found"))
          weak.scar.found <- T
          weak.scar <- difficult.scars[1]
          # print(cat("Scars", difficult.scars, "will be difficult to place"))
        }
      }else{
        # How to relate this to weak scars?
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
  iterative.end <- Sys.time()
  iterative.time <- iterative.end - iterative.start
  print(iterative.time)

  weak.scars <- c(weak.scars, weak.scar)
  weak.scar <- character()
  if(!weak.scar.found){break}
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
# pdf("Images/Simulations/tree_B_wdoublets_doubletrate009_detratio01_branchratio025.pdf",
# width = 20, height = 10)
plot(scar.phylo, show.node.label = F, show.tip.label = F, root.edge = T,
     edge.width = 3, no.margin = T, direction = "leftward")
edgelabels(phylo.edges$Node.2, frame = "none", adj = c(0.5, 0), cex = 2,
           col = "red")
# dev.off()

# View(it.tree.building[[1]]$LLS.unique)

# Place cells in tree ####
print("Placing cells")
# Name nodes
tree.summary.old.pc <- tree.summary.old[tree.summary.old$Node.1 == "Root", ]
tree.summary.old.pc$Node <-
  paste("0", tree.summary.old.pc$Component, sep = "_")
for(d in 1:max(tree.summary.old$Depth)){
  tree.summary.pc.add <- tree.summary.old[tree.summary.old$Depth == d, ]
  tree.summary.pc.add <- merge(tree.summary.pc.add, tree.summary.old.pc[, c("Node.2", "Node")],
                               by.x = "Node.1", by.y = "Node.2")
  tree.summary.pc.add$Node <- 
    paste(tree.summary.pc.add$Node, tree.summary.pc.add$Component, sep = "_")
  tree.summary.old.pc <- rbind(tree.summary.old.pc, tree.summary.pc.add)
}

# Place correct cells in the lowest possible place
correct.cell.placement.positions <- 
  merge(cells.in.tree.f, tree.summary.old[, c("Node.2", "Depth", "Main")],
        by.x = "Scar", by.y = "Node.2")
correct.cell.depths <-
  aggregate(correct.cell.placement.positions$Depth,
            by = list(Cell = correct.cell.placement.positions$Cell),
            max)
colnames(correct.cell.depths)[2] <- "Depth"                      
correct.cell.placement <- merge(correct.cell.placement.positions, correct.cell.depths)
correct.cell.placement <- 
  merge(correct.cell.placement, tree.summary.old.pc[, c("Node.2", "Node")],
        by.x = "Scar", by.y = "Node.2")

# Determine which doublet-flagged cells can be placed: 
# unplaceable cells will not have any scars in the actual tree; of the remainder,
#   incorrect cells will have conflicting scar placements; 
#   correct cells are the remaining cells.
if(length(inc.cells) > 0){
  cells.in.tree.flagged <- cells.in.tree[cells.in.tree$Cell %in% inc.cells, ]
  cells.in.tree.flagged <- 
    merge(cells.in.tree.flagged, tree.summary.old[, c("Node.2", "Depth", "Main")],
          by.x = "Scar", by.y = "Node.2", all.x = T)
  unplaceable.cells <- 
    unique(cells.in.tree.flagged$Cell[is.na(cells.in.tree.flagged$Depth)])
  placeable.cells <- cells.in.tree.flagged[!is.na(cells.in.tree.flagged$Depth), ]
  unplaceable.cells <- setdiff(unplaceable.cells, placeable.cells$Cell)
  placeable.cells <- merge(placeable.cells, tree.summary.old.pc[, c("Node.2", "Node")],
                           by.x = "Scar", by.y = "Node.2")
  correct.conflicting <- unique(placeable.cells[, c("Cell", "Cell.type")])
  correct.conflicting$Node <- NA                               
  correct.conflicting$Node <-
    sapply(correct.conflicting$Cell,
           function(x){
             this.cell <- placeable.cells[placeable.cells$Cell == x, ]
             if(max(table(this.cell$Depth)) > 1){
               return("-1")
             }else{
               this.cell <- this.cell[order(this.cell$Depth), ]
               conflict <- F
               for(i in 2:nrow(this.cell)){
                 if(!grepl(this.cell$Node[i-1], this.cell$Node[i])){
                   conflict <- T
                   break
                 }
               }
               if(conflict){
                 return("-1")
               }else{
                 return(this.cell$Node[nrow(this.cell)])
               }
             }
           }
    )
  really.conflicting <- correct.conflicting[correct.conflicting$Node == "-1", ]
  actually.not.conflicting <- correct.conflicting[correct.conflicting$Node != "-1", ]
  actually.not.conflicting <-
    merge(actually.not.conflicting, tree.summary.old.pc[, c("Node.2", "Depth", "Node", "Main")])
  colnames(actually.not.conflicting)[which(colnames(actually.not.conflicting) == "Node.2")] <-
    "Scar"
  
  # Place placeable doublet-flagged cells
  correct.cell.placement <- rbind(correct.cell.placement, actually.not.conflicting)
}else{
  really.conflicting <- cells.in.tree.f[0, ]
  unplaceable.cells <- character()
  actually.not.conflicting <- cells.in.tree.f[0, ]
}

# Calculate tree statistics
tree.statistics <- data.frame(Cells = length(unique(cells.in.tree$Cell)),
                              Doublets = nrow(really.conflicting),
                              Unplaceable = length(unplaceable.cells),
                              Placeable.main = sum(correct.cell.placement$Main),
                              Placeable.off.main = sum(!correct.cell.placement$Main),
                              Recovered.suspected.doublets = nrow(actually.not.conflicting))

tree.summary.out.2 <- tree.summary.old.pc[, c("Node.2", "Depth", "Main", "Node", "Size")]
colnames(tree.summary.out.2)[1] <- "Scar"

# write.csv(correct.cell.placement,
#           "./Data/2017_10X_2/Z2_tree_sc_positions_dr09_det01_bs025_scarp001_larvae1.csv",
#           row.names = F, quote = F)

# Create tree summary and make node piecharts ####
print("Making tree summary and node piecharts")
# Aggregate cells (stratified by cell type) to nodes to make a cumulative node
# count.
cumulative.node.count <- 
  expand.grid(Node = unique(tree.summary.old.pc$Node),
              Cell.type = unique(correct.cell.placement$Cell.type),
              stringsAsFactors = F)

node.count <- 
  data.frame(table(correct.cell.placement$Node, correct.cell.placement$Cell.type))
colnames(node.count)[1:2] <- c("Node", "Cell.type")
node.count <- merge(node.count, tree.summary.old.pc[, c("Node", "Main")])

cumulative.node.count$Cumulative.count.main <- NA
cumulative.node.count$Cumulative.count.all <- NA
for(i in 1:nrow(cumulative.node.count)){
  c.node <- cumulative.node.count$Node[i]
  c.node.pattern <- paste(c.node, "(_|$)", sep = "")
  c.type <- cumulative.node.count$Cell.type[i]
  
  nodes.under.and.including <-
    node.count[node.count$Cell.type == c.type &
                 grepl(c.node.pattern, node.count$Node), ]
  
  cumulative.node.count$Cumulative.count.main[i] <-
    sum(nodes.under.and.including$Freq[nodes.under.and.including$Main])
  cumulative.node.count$Cumulative.count.all[i] <-
    sum(nodes.under.and.including$Freq)
}

# Calculate total node sizes, make tree summary and calculate cell type ratios 
# per node
tree.summary.out.1 <- aggregate(cumulative.node.count$Cumulative.count.main,
                          by = list(Node = cumulative.node.count$Node),
                          sum)
colnames(tree.summary.out.1)[2] <- "Total.main"
tree.summary.out.2 <- aggregate(cumulative.node.count$Cumulative.count.all,
                                by = list(Node = cumulative.node.count$Node),
                                sum)
colnames(tree.summary.out.2)[2] <- "Total.all"
tree.summary.out <- merge(tree.summary.out.1, tree.summary.out.2)

tree.summary.out <- merge(tree.summary.out, 
                          tree.summary.old.pc[, c("Node", "Node.2", "Depth", "Main")])
colnames(tree.summary.out)[which(colnames(tree.summary.out) == "Node.2")] <- "Scar"
tree.summary.out <- tree.summary.out[, c("Scar", "Node", "Depth", "Total.main", "Total.all", "Main")]

cumulative.node.count <- merge(cumulative.node.count, tree.summary.out)
cumulative.node.count$Ratio.main <-
  cumulative.node.count$Cumulative.count.main/cumulative.node.count$Total.main
cumulative.node.count$Ratio.all <-
  cumulative.node.count$Cumulative.count.all/cumulative.node.count$Total.all

# Output tree summary
# write.csv(tree.summary.out,
#           "./Data/Simulations/Tree_A_1kcellsout_det03_reconstructed_tree.csv",
#           row.names = F, quote = F)

# Plot pie charts main only
# ggplot(cumulative.node.count[cumulative.node.count$Main, ]) +
#   geom_bar(aes(x = "", y = Ratio.main, fill = as.factor(Cell.type)), stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   coord_polar("y", start = 0) +
#   facet_wrap(~ Node) +
#   labs(x = "", y = "") +
#   theme(axis.ticks.y = element_blank(),
#         axis.text.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank())

# Plot pie charts all
# ggplot(cumulative.node.count) +
#   geom_bar(aes(x = "", y = Ratio.all, fill = as.factor(Cell.type)), stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   coord_polar("y", start = 0) +
#   facet_wrap(~ Node) +
#   labs(x = "", y = "") +
#   theme(axis.ticks.y = element_blank(),
#         axis.text.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank())

# Create edgelist of scars and cells ####
# Create edgelist of scars
tree.summary.edge <- tree.summary.out[, c("Node", "Scar")]
colnames(tree.summary.edge) <- c("Child", "Scar.acquisition")
tree.summary.edge$Parent <-
  sapply(tree.summary.edge$Child,
         function(x){
           y <- unlist(strsplit(x, "_"))
           z <- paste(y[-length(y)], collapse = "_")
         }
  )
cell.edge <- correct.cell.placement[, c("Node", "Cell")]
colnames(cell.edge) <- c("Parent", "Child")
cell.edge$Scar.acquisition <- ""

tree.edgelist <- rbind(tree.summary.edge, cell.edge)[, c("Parent", "Child", "Scar.acquisition")]
# write.csv(tree.edgelist, "./Data/Simulations/Tree_B2_2000cell_LINNAEUS_tree.csv",
#           row.names = F, quote = F)

# Do same for collapsed tree
tree.summary.c <- tree.summary[, c("Node.1", "Node.2")]
# Create root node
root.scars <- unique(tree.summary$Node.1[grepl("Root", tree.summary$Node.1)])
root.scars <- paste(unlist(strsplit(root.scars, ","))[-1], collapse = ",")
tree.summary.c$Node.1[grepl("Root", tree.summary$Node.1)] <- root.scars
tree.summary.c <- 
  rbind(data.frame(Node.1 = 0, Node.2 = root.scars, stringsAsFactors = F),
        tree.summary.c)
# Rename nodes and create scar.acquisition column
known.nodes <- 
  data.frame(Old.node = unique(c(tree.summary.c$Node.1, tree.summary.c$Node.2)))
known.nodes$New.node <- 0:(nrow(known.nodes) - 1)
tree.summary.c <- merge(tree.summary.c, known.nodes,
                        by.x = "Node.1", by.y = "Old.node")
colnames(tree.summary.c)[3] <- "Parent"
tree.summary.c <- merge(tree.summary.c, known.nodes,
                        by.x = "Node.2", by.y = "Old.node")
colnames(tree.summary.c)[4] <- "Child"
tree.summary.c <- tree.summary.c[, c("Parent", "Child", "Node.2")]
colnames(tree.summary.c)[3] <- "Scar.acquisition"
# Place cells in tree
cells.to.place <- correct.cell.placement[correct.cell.placement$Main, c("Cell", "Scar")]
# x <- 70
cells.to.place$Parent <-
  sapply(cells.to.place$Scar,
         function(x){
           scar.pattern <- 
             paste("^", x, ",|^", x, "$|,", x, ",|,", x, "$", sep = "")
           parent <- tree.summary.c$Child[grepl(scar.pattern, tree.summary.c$Scar.acquisition)]
           return(parent)
         }
  )
colnames(cells.to.place)[1] <- "Child"
cells.to.place <- cells.to.place[, c("Parent", "Child")]
cells.to.place$Scar.acquisition <- ""

collapsed.tree <- rbind(tree.summary.c, cells.to.place)

if(exists("tsne.coord")){
collapsed.tree <- merge(collapsed.tree, tsne.coord[, c("Cell", "Cell.type")],
                        by.x = "Child", by.y = "Cell", all.x = T)
}else{
  collapsed.tree$Cell.type <- 
    sapply(collapsed.tree$Child,
           function(x) {
             if(grepl("_", x)){
               type.output <- "Cell"
             }else{
               type.output <- NA
             }
             return(type.output)
           }
    )
}

# write.table(collapsed.tree, "./Data/Simulations/Tree_B2_2000cell_LINNAEUS_tree.csv",
#           row.names = F, quote = F, sep = " ")

# Visualize tree ####
# Without cells
tree.summary.c.plot <- tree.summary.c
tree.summary.c.plot$fill <- "black"
tree.summary.c.plot$size <- 1
tree.summary.c.plot$Cell.type <- NA
# Order in scar_tree order
# scar_tree <- read.table("./Data/Simulations/tree_B2_scar_tree.csv",
#                         header = T, fill = T, stringsAsFactors = F)
# tentative.reorder <- data.frame(Scar.acquisition = tree.summary.c.plot$Scar.acquisition)
# tentative.reorder$Order <-
#   sapply(tentative.reorder$Scar.acquisition,
#          function(x) {
#            if(x %in% scar_tree$Scar.acquisition){
#              which(scar_tree$Scar.acquisition == x)
#            }else{0}
#          }
#   )
# tentative.reorder$Order[tentative.reorder$Order == 0] <-
#   c(11, 22, 13)
# tree.summary.c.plot <- merge(tree.summary.c.plot, tentative.reorder)
# tree.summary.c.plot <- tree.summary.c.plot[order(tree.summary.c.plot$Order),
#                                    c("Parent", "Child", "Scar.acquisition",
#                                      "fill", "size", "Cell.type")]
tree.summary.c.plot <- tree.summary.c.plot[order(tree.summary.c.plot$Parent), ]
# save(tree.summary.c.plot, file = "./Data/Simulations/B2_wweak_dlet0_edges.Robj")
LINNAEUS.tree <- generate_tree(tree.summary.c.plot)
# save(LINNAEUS.tree, file = "./Data/Simulations/B2_wweak_dlet0_tree.Robj")
LINNAEUS.tree_wg <-
  collapsibleTree(LINNAEUS.tree, root = LINNAEUS.tree$scar, pieNode = F,
                  pieSummary = F, fill = "fill", nodeSize = "size",
                  collapsed = F, ctypes=unique(LINNAEUS.tree$Get("Cell.types")))
LINNAEUS.tree_wg

# LINNAEUS.tree_wg <- 
#   collapsibleTree(LINNAEUS.tree, root = LINNAEUS.tree$scar, collapsed = F,
#                   fontSize = 8, width = 300, height = 400, fill = "fill",
#                   nodeSize = "size", pieSummary = F)
# save(LINNAEUS.tree_wg, file = "./Data/Simulations/B2_wweak_dlet0_widget.Robj")
# LINNAEUS.tree_wg
# htmlwidgets::saveWidget(LINNAEUS.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_d005_LINNAEUS_tree.html")

# With cells
tree.cells.c.plot <- collapsed.tree
tree.cells.c.plot$fill <-  
  sapply(tree.cells.c.plot$Child,
         function(x){
           if(grepl("_", x)){
             return("lightgrey")
           }else{
             return("black")
           }
         })
tree.cells.c.plot$size <- 
  sapply(tree.cells.c.plot$Child,
         function(x){
           if(grepl("_", x)){
             return(0.5)
           }else{
             return(1)
           }
         })
tree.cells.c.plot <- tree.cells.c.plot[order(tree.cells.c.plot$Parent), ]
LINNAEUS.cell.tree <- generate_tree(tree.cells.c.plot)
# save(LINNAEUS.cell.tree, 
#      file = "./Scripts/linnaeus-scripts/collapsibleTree/sand/C2_correct_tree_wcells.Robj")
# load(file = "./Scripts/linnaeus-scripts/collapsibleTree/sand/C2_correct_tree_wcells.Robj")
collapsibleTree(LINNAEUS.cell.tree, collapsed = F, pieSummary=F, pieNode=F, 
                nodeSize='size', ctypes=unique(LINNAEUS.cell.tree$Get("Cell.type")), 
                fill='fill')
# load("./Scripts/linnaeus-scripts/collapsibleTree/sand/C2_phylip_0_wcells.Robj");
collapsibleTree(phylip.tree, collapsed = F, pieSummary=F, pieNode=F, 
                nodeSize='size', ctypes=unique(phylip.tree$Get("Cell.type")), 
                fill='fill')



LINNAEUS.cell.tree_wg <- 
  collapsibleTree(LINNAEUS.cell.tree, root = LINNAEUS.cell.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 600)
# , fill = "fill",
#                   nodeSize = "size")
LINNAEUS.cell.tree_wg
LINNAEUS.cell.tree.pie <-
  collapsibleTree(LINNAEUS.cell.tree, root = LINNAEUS.cell.tree$scar,
                  collapsed = F, pieNode = T)
LINNAEUS.cell.tree.pie
# htmlwidgets::saveWidget(LINNAEUS.cell.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_C2_03det_iterative.html")
# save(tree.cells.c.plot, file = "./Data/2017_10X_2/Z2_tree_2.Robj")
# save(LINNAEUS.cell.tree, file = "./Data/2017_10X_2/Z2_Ltree_2.Robj")
# save(LINNAEUS.cell.tree.pie, file = "./Data/2017_10X_2/Z2_Ltree_pie.Robj")

# Extract subtree ####
get.node.comp <- function(node, tree.edges){
  # Return the cell type composition of a node, including its subnodes (recursive function)
  node.children <- extracted.tree$Child[is.na(extracted.tree$Cell.type) & 
                                          extracted.tree$Parent == node]
  cell.children <- extracted.tree$Child[!is.na(extracted.tree$Cell.type) & 
                                          extracted.tree$Parent == node]
  comp.this.node <-
    data.frame(table(extracted.tree$Cell.type[extracted.tree$Child %in% 
                                                cell.children]))
  colnames(comp.this.node) <- c("Cell.type", "Count")
  if(length(node.children) > 0){
    for(n in 1:length(node.children)){
      comp.below <- get.node.comp(node.children[n], tree.edges)
      colnames(comp.below)[2] <- "Count.1"
      comp.this.node <- merge(comp.this.node, comp.below)
      comp.this.node$Count <- comp.this.node$Count + comp.this.node$Count.1
      comp.this.node <- comp.this.node[, c("Cell.type", "Count")]
    }
  }
  
  return(comp.this.node)
}

celltypes.to.extract <- c("Chondrocytes A", "Retinal cells A", "Erythrocytes B",
                          "Epidermal cells B", "Fibroblasts B (Fin)", "Hepatocytes A")
extracted.tree <- tree.cells.c.plot
extracted.tree <- extracted.tree[is.na(extracted.tree$Cell.type) | 
                                   extracted.tree$Cell.type %in% celltypes.to.extract, ]
extracted.tree$Keep <-
  sapply(extracted.tree$Child,
         function(x){
           if(grepl("_", x)){
             return(T)
           }else{
             node.comp <- get.node.comp(x, extracted.tree)
             node.comp <- node.comp[node.comp$Cell.type %in% celltypes.to.extract, ]
             
             return(sum(node.comp$Count) > 0)
           }
         }
  )

extracted.tree <- extracted.tree[extracted.tree$Keep, ]

# Tree vis tests ####
tree.summary.c.plot <- tree.summary.c
tree.summary.c.plot$fill <- "black"
tree.summary.c.plot$size <- 1
tree.summary.c.plot$Cell.type <- NA
tree.summary.c.plot <- tree.summary.c.plot[order(tree.summary.c.plot$Parent), ]

LINNAEUS.tree <- generate_tree(tree.summary.c.plot)
# save(LINNAEUS.tree, file = "./Scripts/linnaeus-scripts/collapsibleTree/sand/B2_correct_tree.Robj")
collapsibleTree(LINNAEUS.tree, pieNode=T, pieSummary=F, collapsed=F, 
                width=500, height=300, ctypes=c('internal'))
collapsibleTree(LINNAEUS.tree, pieNode=T, pieSummary=F, collapsed=F, 
                width=1.5e3, height=1e3, nodeLabel_sc=40, ctypes=c('internal'))
collapsibleTree(LINNAEUS.tree, root = LINNAEUS.tree$scar, collapsed = F,  
                pieNode=T, ctypes=unique(LINNAEUS.tree$Get("Cell.type")))

LINNAEUS.tree_wg <-
  collapsibleTree(LINNAEUS.tree, root = LINNAEUS.tree$scar, pieNode = F, 
                  collapsed = F, ctypes=unique(LINNAEUS.tree$Get("Cell.types")))
LINNAEUS.tree_wg

# load(file = "./Scripts/linnaeus-scripts/collapsibleTree/sand/B2_correct_tree.Robj")
# collapsibleTree(LINNAEUS.tree, collapsed = F, pieSummary=F, pieNode=F, 
#                 nodeSize='size', ctypes=unique(LINNAEUS.tree$Get("Cell.type")), 
#                 fill='fill')

# Data structure tests ####
basic.structure <- it.tree.building[[1]]$LLS
ggplot(basic.structure) +
  geom_point(aes(x = Scar.count, y = Degree))
basic.structure$Scars.per.deg <- basic.structure$Scar.count/basic.structure$Degree
ggplot(basic.structure) +
  geom_histogram(aes(x = Scars.per.deg))
