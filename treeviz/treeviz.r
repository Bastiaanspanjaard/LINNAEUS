# Dependencies ####
source("./Scripts/linnaeus-scripts/scar_helper_functions.R")
# require(data.tree)
# require(treemap)
# source("./Scripts/linnaeus-scripts/collapsibleTree.R")
# source("./Scripts/linnaeus-scripts/collapsibleTree.data.tree.R")
# 
# generate_tree = function(df){
#   # columns.include <- c("Parent", "Child", "Scar.acquisition")
#   # if(!is.null(fill.col)){
#   #   columns.include <- c(columns.include, fill.col)
#   # }
#   # if(!is.null(size.col)){ columns.include <- c(columns.include, size.col)}
#   for(i in 1:nrow(df)){
#     parent = paste0('nd', as.character(df$Parent[i]))
#     child = paste0('nd', as.character(df$Child[i]))
#     scar <- df$Scar.acquisition[i]
# 
#     if(!exists(child)){
#       eval_txt = sprintf('%s <<- Node$new("%s", name="%s", scar = "%s")',
#                          child, child, child, scar)
#       eval(parse(text=eval_txt))
#       if("fill" %in% colnames(df)){
#         eval_txt <- paste(child, "$fill <- \"", df$fill[i], "\"", sep = "")
#         eval(parse(text = eval_txt))
#       }
#       if("size" %in% colnames(df)){
#         eval_txt <- paste(child, "$size <- ", df$size[i], sep = "")
#         eval(parse(text = eval_txt))
#       }
#     }
#     if(exists(parent)){
#       add_txt = sprintf('%s$AddChildNode(%s)', parent, child)
#       eval(parse(text=add_txt))
#     }
#   }
# 
#   return_tree <- eval(parse(text=sprintf('%s$root', ls(envir=globalenv(), pattern='^nd')[1])))
#   rm(list=ls(envir=globalenv(), pattern='^nd'), envir=globalenv())
# 
#   return(return_tree)
# }
    
# Create developmental tree C ####
dev_tree <- read.table("./Data/Simulations/tree_C_dev_tree_2.txt", 
                       header = T, fill = T, stringsAsFactors = F)
dev_tree$fill <- "black"
dev_tree$size <- 1.5
dev_tree_C <- generate_tree(dev_tree)
dev_tree_wg <- collapsibleTree(dev_tree_C, root = dev_tree_C$scar, collapsed = F,
                               fontSize = 12, width = 300, height = 200,
                               fill = "fill", nodeSize = "size")
dev_tree_wg
# htmlwidgets::saveWidget(dev_tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_C2_dev_tree.html")

# Create CS 0 tree C ####
phylip.edges <- read.table("./Data/Simulations/Tree_C2_100cells_03det_0_PHYLIP_tree1",
                                     stringsAsFactors = F)
phylip.scars <- read.csv("./Data/Simulations/Tree_C2_scar_conversion.csv",
                            stringsAsFactors = F)

# Add scar acquisition
phylip.edges$Scars <- do.call(paste0, phylip.edges[-(1:3)])
phylip.edges$Scar.acquisition <- 
  sapply(phylip.edges$Scars,
         function(x) paste(phylip.scars$Scar[unlist(strsplit(x, "")) == "1"], 
                           collapse = ","))

# Change cell names (remove x's)
phylip.edges$Child <- 
  sapply(phylip.edges$V2,
         function(x){
           if(grepl("_", x)){
             return(unlist(strsplit(x, "x"))[1])
           }else{
             return(x)
           }
         }
  )

# Get correctly-formatted edge list with scars
phylip.edges.f <- phylip.edges[, c("V1", "Child", "Scar.acquisition")]
colnames(phylip.edges.f)[1] <- "Parent"

# Collapse empty (i.e. scar-less) nodes, unless that node is a cell.
phylip.edges.collapse <- phylip.edges.f
# Loop over dataframe to find all nodes and tips that did not get a scar
index <- 1
while(index <= nrow(phylip.edges.collapse)){
  if(phylip.edges.collapse$Scar.acquisition[index] == "" &
     !grepl("_", phylip.edges.collapse$Child[index])){
    cell.wo.scar.name <- phylip.edges.collapse$Child[index]
    upstream.name <- phylip.edges.collapse$Parent[index]
    phylip.edges.collapse$Parent[phylip.edges.collapse$Parent == 
                                        cell.wo.scar.name] <-
      upstream.name
    phylip.edges.collapse <- phylip.edges.collapse[-index, ]
  }else{index <- index + 1}
}

# Create nodes before cells if a scar is created in a cell
e <- 1
while(e <= nrow(phylip.edges.collapse)){
  # For every edge, test if a scar is created in a cell. If so, create a new
  # node to create the scar in that has the cell as daughter.
  if(phylip.edges.collapse$Scar.acquisition[e] != "" & 
     grepl("_", phylip.edges.collapse$Child[e])){
    # Determine the cell
    cell.child <- phylip.edges.collapse$Child[e]
    
    # Determine the new node number
    numbered.parent.nodes <- 
      as.integer(phylip.edges.collapse$Parent[which(phylip.edges.collapse$Parent != "root" & 
                                         !grepl("_", phylip.edges.collapse$Parent))])
    numbered.child.nodes <- 
      as.integer(phylip.edges.collapse$Child[which(phylip.edges.collapse$Child != "root" & 
                                         !grepl("_", phylip.edges.collapse$Child))])
    new.node <- max(c(numbered.child.nodes, numbered.parent.nodes)) + 1
    
    # Determine scar
    node.scar <- phylip.edges.collapse$Scar.acquisition[e]
    
    # Create and add new edge
    new.edge <- data.frame(Parent = new.node,
                           Child = cell.child,
                           Scar.acquisition = "")
    phylip.edges.collapse <- rbind(phylip.edges.collapse, new.edge)
    
    # Change old edge
    phylip.edges.collapse$Child[e] <- new.node
  }
  e <- e + 1
}

# phylip.edges.collapse <- 
#   rbind(data.frame(Parent = 0, Child = "root", Scar.acquisition = ""),
#         phylip.edges.collapse)
phylip.edges.collapse$fill <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return("lightgrey")
           }else{
             return("black")
           }
  })
phylip.edges.collapse$size <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return(0.5)
           }else{
             return(1)
           }
         })
phylip.edges.collapse$Cell.type <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return("Cell")
           }else{
             return(NA)
           }
         })
# Add entries for nodesize and fill; get fields into tree (fill has to be named
# "fill", nodesize can be anything but its name has to be supplied in the 
# collapsibleTree functioncall.
phylip.tree <- generate_tree(phylip.edges.collapse)
# save(phylip.tree,
     # file = "./Scripts/linnaeus-scripts/collapsibleTree/sand/C2_phylip_0_wcells.Robj")

phylip.tree_wg <- 
  collapsibleTree(phylip.tree, root = phylip.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 800, fill = "fill",
                  nodeSize = "size")
phylip.tree_wg
# htmlwidgets::saveWidget(phylip.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_C2_03det_CamSok0.html")


# Create CS ? tree C ####
phylip.edges <- read.table("./Data/Simulations/Tree_C2_100cells_03det_phylip_q_tree1",
                           stringsAsFactors = F)
phylip.scars <- read.csv("./Data/Simulations/Tree_C2_scar_conversion.csv",
                         stringsAsFactors = F)

# Add scar acquisition
phylip.edges$Scars <- do.call(paste0, phylip.edges[-(1:3)])
phylip.edges$Scar.acquisition <- 
  sapply(phylip.edges$Scars,
         function(x) paste(phylip.scars$Scar[unlist(strsplit(x, "")) == "1"], 
                           collapse = ","))
phylip.edges$Scar.acquisition[grepl("\\?", phylip.edges$Scars)] <- "?"

# Change cell names (remove x's)
phylip.edges$Child <- 
  sapply(phylip.edges$V2,
         function(x){
           if(grepl("_", x)){
             return(unlist(strsplit(x, "x"))[1])
           }else{
             return(x)
           }
         }
  )

# Get correctly-formatted edge list with scars
phylip.edges.f <- phylip.edges[, c("V1", "Child", "Scar.acquisition")]
colnames(phylip.edges.f)[1] <- "Parent"

# Collapse empty (i.e. scar-less) nodes, unless that node is a cell.
phylip.edges.collapse <- phylip.edges.f
# Loop over dataframe to find all nodes and tips that did not get a scar
index <- 1
while(index <= nrow(phylip.edges.collapse)){
  if(phylip.edges.collapse$Scar.acquisition[index] == "" &
     !grepl("_", phylip.edges.collapse$Child[index])){
    cell.wo.scar.name <- phylip.edges.collapse$Child[index]
    upstream.name <- phylip.edges.collapse$Parent[index]
    phylip.edges.collapse$Parent[phylip.edges.collapse$Parent == 
                                   cell.wo.scar.name] <-
      upstream.name
    phylip.edges.collapse <- phylip.edges.collapse[-index, ]
  }else{index <- index + 1}
}

# Create nodes before cells if a scar is created in a cell
e <- 1
while(e <= nrow(phylip.edges.collapse)){
  # For every edge, test if a scar is created in a cell. If so, create a new
  # node to create the scar in that has the cell as daughter.
  if(phylip.edges.collapse$Scar.acquisition[e] != "" & 
     grepl("_", phylip.edges.collapse$Child[e])){
    # Determine the cell
    cell.child <- phylip.edges.collapse$Child[e]
    
    # Determine the new node number
    numbered.parent.nodes <- 
      as.integer(phylip.edges.collapse$Parent[which(phylip.edges.collapse$Parent != "root" & 
                                                      !grepl("_", phylip.edges.collapse$Parent))])
    numbered.child.nodes <- 
      as.integer(phylip.edges.collapse$Child[which(phylip.edges.collapse$Child != "root" & 
                                                     !grepl("_", phylip.edges.collapse$Child))])
    new.node <- max(c(numbered.child.nodes, numbered.parent.nodes)) + 1
    
    # Determine scar
    node.scar <- phylip.edges.collapse$Scar.acquisition[e]
    
    # Create and add new edge
    new.edge <- data.frame(Parent = new.node,
                           Child = cell.child,
                           Scar.acquisition = "")
    phylip.edges.collapse <- rbind(phylip.edges.collapse, new.edge)
    
    # Change old edge
    phylip.edges.collapse$Child[e] <- new.node
  }
  e <- e + 1
}

# phylip.edges.collapse <- 
#   rbind(data.frame(Parent = 0, Child = "root", Scar.acquisition = ""),
#         phylip.edges.collapse)
# Add entries for nodesize and fill; get fields into tree (fill has to be named
# "fill", nodesize can be anything but its name has to be supplied in the 
# collapsibleTree functioncall.
phylip.edges.collapse$fill <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return("lightgrey")
           }else{
             return("black")
           }
         })
phylip.edges.collapse$size <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return(0.5)
           }else{
             return(1)
           }
         })
phylip.edges.collapse$Cell.type <-  
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             "Cell"
           }else{
             NA
           }
         })
phylip.tree <- generate_tree(phylip.edges.collapse)
phylip.tree_wg <- 
  collapsibleTree(phylip.tree, root = phylip.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 800, fill = "fill",
                  nodeSize = "size", pieSummary = F)
phylip.tree_wg
# htmlwidgets::saveWidget(phylip.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_C2_03det_CamSok?.html")



# Create LINNAEUS tree C ####
iterative.edges <- read.csv("./Data/Simulations/Tree_C2_100cell_03det_iterative_tree.csv",
                            stringsAsFactors = F)
iterative.edges$Scar.acquisition[is.na(iterative.edges$Scar.acquisition)] <- ""
iterative.edges$fill <-  
  sapply(iterative.edges$Child,
         function(x){
           if(substring(x, 1, 1) != 0){
             return("lightgrey")
           }else{
             return("black")
           }
         })
iterative.edges$size <- 
  sapply(iterative.edges$Child,
         function(x){
           if(substring(x, 1, 1) != 0){
             return(0.5)
           }else{
             return(1)
           }
         })
iterative.tree <- generate_tree(iterative.edges)
iterative.tree_wg <- 
  collapsibleTree(iterative.tree, root = iterative.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 500, fill = "fill",
                  nodeSize = "size")
iterative.tree_wg
# htmlwidgets::saveWidget(iterative.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_C2_03det_iterative.html")

  




# Create developmental tree B ####
dev_tree <- read.table("./Data/Simulations/tree_B2_dev_tree.txt", 
                       header = T, fill = T, stringsAsFactors = F)
dev_tree$fill <- "black"
dev_tree$size <- 1
dev_tree_B <- generate_tree(dev_tree)
dev_tree_wg <- collapsibleTree(dev_tree_B, root = dev_tree_B$scar, collapsed = F,
                               fontSize = 8, width = 300, height = 400,
                               fill = "fill", nodeSize = "size")
dev_tree_wg
# htmlwidgets::saveWidget(dev_tree_wg,
                        # file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_dev_tree.html")

# Create scar tree B ####
scar_tree <- read.table("./Data/Simulations/tree_B2_scar_tree.csv", 
                       header = T, fill = T, stringsAsFactors = F)
# scar_tree <- scar_tree[order(scar_tree$Parent), ]
# Sort scar tree to conform with the developmental tree as good as possible

scar_tree$fill <- "black"
scar_tree$size <- 1
scar_tree_B <- generate_tree(scar_tree)
scar_tree_wg <- collapsibleTree(scar_tree_B, root = scar_tree_B$scar, collapsed = F,
                               fontSize = 8, width = 300, height = 400,
                               fill = "fill", nodeSize = "size")
scar_tree_wg
# htmlwidgets::saveWidget(scar_tree_wg,
# file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_scar_tree.html")

# Create LINNAEUS tree B ####
iterative.edges <- read.csv("./Data/Simulations/Tree_B2_2000cell_LINNAEUS_tree.csv",
                            stringsAsFactors = F, sep = " ")
iterative.edges$Scar.acquisition[is.na(iterative.edges$Scar.acquisition)] <- ""
iterative.edges <- iterative.edges[!grepl("_", iterative.edges$Child), ]
iterative.edges$fill <-  "black"
iterative.edges$size <- 1

iterative.edge.reorder <- which(iterative.edges$Scar.acquisition == scar_tree$Scar.acquisition)
tentative.reorder <- data.frame(Scar.acquisition = iterative.edges$Scar.acquisition)
tentative.reorder$Order <- 
  sapply(tentative.reorder$Scar.acquisition,
       function(x) {
         if(x %in% scar_tree$Scar.acquisition){
           which(scar_tree$Scar.acquisition == x)
         }else{0}
       }
)
tentative.reorder$Order[tentative.reorder$Order == 0] <-
  c(6, 13, 23, 28, 27)
iterative.edges <- merge(iterative.edges, tentative.reorder)
iterative.edges <- iterative.edges[order(iterative.edges$Order),
                                   c("Parent", "Child", "Scar.acquisition", 
                                     "fill", "size")]
iterative.tree <- generate_tree(iterative.edges)
iterative.tree_wg <- 
  collapsibleTree(iterative.tree, root = iterative.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 400, fill = "fill",
                  nodeSize = "size")
iterative.tree_wg
# htmlwidgets::saveWidget(iterative.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_LINNAEUS_tree.html")







# Create CS 0 tree B ####
phylip.edges <- read.table("./Data/Simulations/Tree_B2_2000cellsout_phylip_0_tree1",
                           stringsAsFactors = F)
phylip.scars <- read.csv("./Data/Simulations/Tree_B2_scar_conversion.csv",
                         stringsAsFactors = F)

# Add scar acquisition
phylip.edges$Scars <- do.call(paste0, phylip.edges[-(1:3)])
phylip.edges$Scar.acquisition <- 
  sapply(phylip.edges$Scars,
         function(x) paste(phylip.scars$Scar[unlist(strsplit(x, "")) == "1"], 
                           collapse = ","))

# Change cell names (remove x's)
phylip.edges$Child <- 
  sapply(phylip.edges$V2,
         function(x){
           if(grepl("_", x)){
             return(unlist(strsplit(x, "x"))[1])
           }else{
             return(x)
           }
         }
  )

# Get correctly-formatted edge list with scars
phylip.edges.f <- phylip.edges[, c("V1", "Child", "Scar.acquisition")]
colnames(phylip.edges.f)[1] <- "Parent"

# Collapse empty (i.e. scar-less) nodes, unless that node is a cell.
phylip.edges.collapse <- phylip.edges.f
# Loop over dataframe to find all nodes and tips that did not get a scar
index <- 1
while(index <= nrow(phylip.edges.collapse)){
  if(phylip.edges.collapse$Scar.acquisition[index] == "" &
     !grepl("_", phylip.edges.collapse$Child[index])){
    cell.wo.scar.name <- phylip.edges.collapse$Child[index]
    upstream.name <- phylip.edges.collapse$Parent[index]
    phylip.edges.collapse$Parent[phylip.edges.collapse$Parent == 
                                   cell.wo.scar.name] <-
      upstream.name
    phylip.edges.collapse <- phylip.edges.collapse[-index, ]
  }else{index <- index + 1}
}

# Create nodes before cells if a scar is created in a cell
e <- 1
while(e <= nrow(phylip.edges.collapse)){
  # For every edge, test if a scar is created in a cell. If so, create a new
  # node to create the scar in that has the cell as daughter.
  if(phylip.edges.collapse$Scar.acquisition[e] != "" & 
     grepl("_", phylip.edges.collapse$Child[e])){
    # Determine the cell
    cell.child <- phylip.edges.collapse$Child[e]
    
    # Determine the new node number
    numbered.parent.nodes <- 
      as.integer(phylip.edges.collapse$Parent[which(phylip.edges.collapse$Parent != "root" & 
                                                      !grepl("_", phylip.edges.collapse$Parent))])
    numbered.child.nodes <- 
      as.integer(phylip.edges.collapse$Child[which(phylip.edges.collapse$Child != "root" & 
                                                     !grepl("_", phylip.edges.collapse$Child))])
    new.node <- max(c(numbered.child.nodes, numbered.parent.nodes)) + 1
    
    # Determine scar
    node.scar <- phylip.edges.collapse$Scar.acquisition[e]
    
    # Create and add new edge
    new.edge <- data.frame(Parent = new.node,
                           Child = cell.child,
                           Scar.acquisition = "")
    phylip.edges.collapse <- rbind(phylip.edges.collapse, new.edge)
    
    # Change old edge
    phylip.edges.collapse$Child[e] <- new.node
  }
  e <- e + 1
}

phylip.edges.collapse <- 
  rbind(data.frame(Parent = 0, Child = "root", Scar.acquisition = ""),
        phylip.edges.collapse)
phylip.edges.collapse$fill <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return("lightgrey")
           }else{
             return("black")
           }
         })
phylip.edges.collapse$size <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return(0.5)
           }else{
             return(1)
           }
         })
# Add entries for nodesize and fill; get fields into tree (fill has to be named
# "fill", nodesize can be anything but its name has to be supplied in the 
# collapsibleTree functioncall.
phylip.tree <- 
  generate_tree(phylip.edges.collapse[!grepl("_", phylip.edges.collapse$Child), ])
phylip.tree_wg <- 
  collapsibleTree(phylip.tree, root = phylip.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 800, fill = "fill",
                  nodeSize = "size")
phylip.tree_wg
# htmlwidgets::saveWidget(phylip.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_CamSok0.html")




# Create doublet CS 0 tree B ####
phylip.edges <- read.table("./Data/Simulations/Tree_B2_2000cellsout_d005_phylip_tree1",
                           stringsAsFactors = F)
phylip.scars <- read.csv("./Data/Simulations/Tree_B2_scar_conversion.csv",
                         stringsAsFactors = F)

# Add scar acquisition
phylip.edges$Scars <- do.call(paste0, phylip.edges[-(1:3)])
phylip.edges$Scar.acquisition <- 
  sapply(phylip.edges$Scars,
         function(x) paste(phylip.scars$Scar[unlist(strsplit(x, "")) == "1"], 
                           collapse = ","))

# Change cell names (remove x's)
phylip.edges$Child <- 
  sapply(phylip.edges$V2,
         function(x){
           if(grepl("_", x)){
             return(unlist(strsplit(x, "x"))[1])
           }else{
             return(x)
           }
         }
  )

# Get correctly-formatted edge list with scars
phylip.edges.f <- phylip.edges[, c("V1", "Child", "Scar.acquisition")]
colnames(phylip.edges.f)[1] <- "Parent"

# Collapse empty (i.e. scar-less) nodes, unless that node is a cell.
phylip.edges.collapse <- phylip.edges.f
# Loop over dataframe to find all nodes and tips that did not get a scar
index <- 1
while(index <= nrow(phylip.edges.collapse)){
  if(phylip.edges.collapse$Scar.acquisition[index] == "" &
     !grepl("_", phylip.edges.collapse$Child[index])){
    cell.wo.scar.name <- phylip.edges.collapse$Child[index]
    upstream.name <- phylip.edges.collapse$Parent[index]
    phylip.edges.collapse$Parent[phylip.edges.collapse$Parent == 
                                   cell.wo.scar.name] <-
      upstream.name
    phylip.edges.collapse <- phylip.edges.collapse[-index, ]
  }else{index <- index + 1}
}

# Create nodes before cells if a scar is created in a cell
e <- 1
while(e <= nrow(phylip.edges.collapse)){
  # For every edge, test if a scar is created in a cell. If so, create a new
  # node to create the scar in that has the cell as daughter.
  if(phylip.edges.collapse$Scar.acquisition[e] != "" & 
     grepl("_", phylip.edges.collapse$Child[e])){
    # Determine the cell
    cell.child <- phylip.edges.collapse$Child[e]
    
    # Determine the new node number
    numbered.parent.nodes <- 
      as.integer(phylip.edges.collapse$Parent[which(phylip.edges.collapse$Parent != "root" & 
                                                      !grepl("_", phylip.edges.collapse$Parent))])
    numbered.child.nodes <- 
      as.integer(phylip.edges.collapse$Child[which(phylip.edges.collapse$Child != "root" & 
                                                     !grepl("_", phylip.edges.collapse$Child))])
    new.node <- max(c(numbered.child.nodes, numbered.parent.nodes)) + 1
    
    # Determine scar
    node.scar <- phylip.edges.collapse$Scar.acquisition[e]
    
    # Create and add new edge
    new.edge <- data.frame(Parent = new.node,
                           Child = cell.child,
                           Scar.acquisition = "")
    phylip.edges.collapse <- rbind(phylip.edges.collapse, new.edge)
    
    # Change old edge
    phylip.edges.collapse$Child[e] <- new.node
  }
  e <- e + 1
}

phylip.edges.collapse <- 
  rbind(data.frame(Parent = 0, Child = "root", Scar.acquisition = ""),
        phylip.edges.collapse)
phylip.edges.collapse$fill <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return("lightgrey")
           }else{
             return("black")
           }
         })
phylip.edges.collapse$size <- 
  sapply(phylip.edges.collapse$Child,
         function(x){
           if(grepl("_", x)){
             return(0.5)
           }else{
             return(1)
           }
         })
# Add entries for nodesize and fill; get fields into tree (fill has to be named
# "fill", nodesize can be anything but its name has to be supplied in the 
# collapsibleTree functioncall.
phylip.tree <- 
  generate_tree(phylip.edges.collapse[!grepl("_", phylip.edges.collapse$Child), ])
phylip.tree_wg <- 
  collapsibleTree(phylip.tree, root = phylip.tree$scar, collapsed = F,
                  fontSize = 8, width = 300, height = 800, fill = "fill",
                  nodeSize = "size")
phylip.tree_wg
# htmlwidgets::saveWidget(phylip.tree_wg,
#                         file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_d005_CamSok0.html")
