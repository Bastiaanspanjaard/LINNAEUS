#' @rdname collapsibleTree
#' @method collapsibleTree Node
#' @export
rename.node <- function(node){
  if("scar" %in% names(node)){
    node$name <- node$scar
  }else{
    node$name <- ""
  }

  if("children" %in% names(node)){
    for(i in 1:length(node$children)){
      node$children[[i]] <- rename.node(node$children[[i]])
    }
  }

  return(node)
}

collapsibleTree.Node <- function(df, hierarchy_attribute = "level",
                                 root = df$name, inputId = NULL, attribute = "leafCount",
                                 aggFun = sum, fill = "lightsteelblue",
                                 linkLength = NULL, fontSize = 10, tooltip = FALSE,
                                 tooltipHtml = NULL,nodeSize = NULL, collapsed = TRUE,
                                 zoomable = TRUE, width = NULL, height = NULL,
				# PO
    				nodeSize_sc = 2, nodeLabel_sc = FALSE,
				colors = NULL, ctypes = NULL, 
				nodeSize_class = c(   10, 15, 20, 35),
				nodeSize_breaks = c( 0, 5, 20, 100, 1e6),
				pieSummary = TRUE,
				pieNode = FALSE,  
				...) {


  # acceptable inherent node attributes
  nodeAttr <- c("leafCount", "count")

  # reject bad inputs
  if(!is(df) %in% "Node") stop("df must be a data tree object")
  if(!is.character(fill)) stop("fill must be a either a color or column name")
  if(!is.null(tooltipHtml)) if(!(tooltipHtml %in% df$fields)) stop("tooltipHtml column name is incorrect")
  if(!is.null(nodeSize)) if(!(nodeSize %in% c(df$fields, nodeAttr))) stop("nodeSize column name is incorrect")

  # calculate the right and left margins in pixels
  leftMargin <- nchar(root)
  rightLabelVector <- df$Get("name", filterFun = function(x) x$level==df$height)
  rightMargin <- max(sapply(rightLabelVector, nchar))

  # Deriving hierarchy variable from data.tree input
  hierarchy <- unique(ToDataFrameTree(df, hierarchy_attribute)[[hierarchy_attribute]])
  if(length(hierarchy) <= 1) stop("hierarchy vector must be greater than length 1")
  # PO tmp code for dealing with ctype colour attribution
  # PO TODO, something still looks weird on the pie charts. Colors are not sorted.
  # ctypes = sort(unique(df$Get("Cell.type")))
  # PO create colors
  if(is.null(colors)){
	colors = c('#7400C2', '#8E14E0', '#A929FF', '#BB56FF', '#CD82FF', '#DFAEFF', '#F1DBFF', '#844420', '#904A23', '#9D5126', '#A85A2D', '#AE673D', '#B5734D', '#BB805D', '#C28C6D', '#C9997E', '#944D72', '#DF74AC', '#F89ACB', '#FABFDE', '#FDE5F2', '#173416', '#214B1F', '#2B6229', '#357933', '#3F913D', '#4AA847', '#57B354', '#66BA64', '#75C173', '#84C782', '#93CE91', '#A2D5A0', '#B1DCAF', '#C0E2BE', '#CFE9CD', '#DEF0DC', '#EDF7EC', '#B75B01', '#ED7600', '#FF9833', '#FFBE7F', '#FFE5CC', '#99991E', '#E6E62D', '#FFFF57', '#FFFF8C', '#FFFFC1', '#9F1213', '#B91516', '#D31819', '#E52729', '#EA4E4F', '#EE7475', '#F39B9B', '#F7C1C1', '#FCE8E8', '#214B6E', '#265780', '#2C6493', '#3171A6', '#377EB8', '#4E8DC0', '#649BC7', '#7BA9CF', '#91B8D7', '#A7C6DF', '#BED5E7', '#D4E3EF', '#EBF2F7')
  }

  if(is.null(ctypes)){
	ctypes = c('Hepatocytes A', 'Intestinal cells A', 'Epithelial cells (Mucosa)', 'Intestinal cells B (Exocrine)', 'Pancreas exocrine cells', 'Hepatocytes B', 'Hepatocytes C', 'Erythrocytes A', 'Erythrocytes B', 'Lymphocytes A', 'Nephric duct cells', 'Endothelial cells', 'Macrophages', 'Neutrophils', 'Kidney cells', 'Lymphocytes B', 'Pigment cells (Xantophores)', 'Pigment cells (Melanocytes)', 'Pigment cells (Iridophores)', 'Glial cells (Peripheral)', 'Neurons E (Peripheral)', 'Neurons A (Brain)', 'Neuronal precursors A', 'Neurons B (Proliferating)', 'Cone cells', 'Retinal cells A', 'Neuronal precursors B', 'Neurons C (Spinal cord)', 'Retinal cells B', 'Rod cells', 'Neurons D (Brain)', 'Radial glia A', 'Retinal pigment epithelial cells', 'Radial glia B', 'Olfactory receptor cells', 'Epithelial cells (Otic vesicle)', 'Oligodendrocytes', 'Inner ear cells', 'Chondrocytes A', 'Chondrocytes B', 'Chondrocytes C', 'Osteoblasts', 'Chondrocytes D', 'Undifferentiated mesodermal cells', 'Fibroblasts A', 'Fibroblasts B (Fin)', 'Fibroblasts C', 'Fibroblasts D', 'Skeletal muscle cells A', 'Muscle cells (Pectoral fin)', 'Skeletal muscle cells B', 'Satellite cells', 'Smooth muscle cells A (Vasculature)', 'Skeletal muscle cells C', 'Skeletal muscle cells D', 'Smooth muscle B', 'Skeletal muscle cells E (slow)', 'Epidermal cells A', 'Epithelial cells (Fin)', 'Epidermal cells B', 'Lens cells A', 'Keratinocytes A', 'Epidermal cells C', 'Lens cells B', 'Epidermal cells D', 'Neuromast cells', 'Epidermal precursors', 'Corneal cells', 'Epithelial cells (other)', 'Keratinocytes B')
  }
  ctypes = ctypes[!is.na(ctypes)]
  ctypes = ctypes[is.na(match(ctypes, "NA"))]



  # create a list that contains the options
  options <- list(
    hierarchy = hierarchy,
    input = inputId,
    attribute = attribute,
    linkLength = linkLength,
    fontSize = fontSize,
    tooltip = tooltip,
    collapsed = collapsed,
    pieNode = pieNode, #PO 
    useColors = !is.null(colors), #PO 
    colors = colors, #PO
    nodeLabel_sc = ifelse(is.null(nodeLabel_sc), TRUE, nodeLabel_sc ),
    zoomable = zoomable,
    margin = list(
      top = 20,
      bottom = 20,
      left = (leftMargin * fontSize/2) + 25,
      right = (rightMargin * fontSize/2) + 25
    )
  )

  # these are the fields that will ultimately end up in the json
  jsonFields <- "scar"

  if(fill %in% df$fields) {
    # fill in node colors based on column name
    df$Do(function(x) x$fill <- x[[fill]])
    jsonFields <- c(jsonFields, "fill")
  } else {
    # default to using fill value as literal color name
    options$fill <- fill
  }
  # PO determine size classes
  df = Clone(df)
  #if(pieSummary){
  #  df = Clone(df)
  #}
  SortNumeric(df, decreasing=T, recursive=T,  attribute = function(x){ifelse(x$Cell.type == "NA", "1e4", as.numeric(match(x$Cell.type, ctypes)))})

  if(pieNode){
    t <- data.tree::Traverse(df, 'level')
    data.tree::Do(t, function(x) {
	x$isScar = !x$isLeaf & !x$isRoot
	if(x$isRoot) {x$isScar = TRUE; x$Cell.type = "_"}
	xpieNode = x$Get("Cell.type")
	x$ct = x$Cell.type
	x$pieNode = table(factor(array(xpieNode), levels=ctypes))
	x$SizeOfNode = nodeSize_class[cut(sum(x$pieNode), breaks=nodeSize_breaks, include.lowest=T, labels=F )]
	if(!x$isScar) {
		x$SizeOfNode = nodeSize_sc
	}else{
		sapply(x$children, function(child){child$parSize = length(x$children)}) 
	}
    })
    jsonFields <- c(jsonFields, "pieNode")
    jsonFields <- c(jsonFields, "parSize") # keeps memory of the size of parent; used to decide wheter to show cell type of single cell
    jsonFields <- c(jsonFields, "ct")
    jsonFields <- c(jsonFields, "SizeOfNode")
    jsonFields <- c(jsonFields, "isScar")
  }

  if(pieSummary){
   # Only after collecting the statistics for the scar nodes we get rid of the scells
   	t <- data.tree::Traverse(df, 'post-order')
   	data.tree::Do(t, function(x) {
		if(x$isLeaf & !x$isRoot & !x$isScar){	
			x$parent$RemoveChild(x$name)
		}
	})
  }

  # only necessary to perform these calculations if there is a tooltip
  if(tooltip & is.null(tooltipHtml)) {
    t <- data.tree::Traverse(df, hierarchy_attribute)
    if(substitute(identity)=="identity") {
      # for identity, leave the tooltips as is
      data.tree::Do(t, function(x) {
        x$WeightOfNode <- x[[attribute]]
      })
    } else {
      # traverse down the tree and compute the weights of each node for the tooltip
      data.tree::Do(t, function(x) {
        x$WeightOfNode <- data.tree::Aggregate(x, attribute, aggFun)
        # make the tooltips look nice
        x$WeightOfNode <- prettyNum(
          x$WeightOfNode, big.mark = ",", digits = 3, scientific = FALSE
        )
      })
    }
    jsonFields <- c(jsonFields, "WeightOfNode")
  }

  # if tooltipHtml is specified, pass it on in the data
  if(tooltip & !is.null(tooltipHtml)) {
    df$Do(function(x) x$tooltip <- x[[tooltipHtml]])
    jsonFields <- c(jsonFields, "tooltip")
  }

  # only necessary to perform these calculations if there is a nodeSize specified
  if(!is.null(nodeSize) & !pieNode) {
    # Scale factor to keep the median leaf size around 10
    scaleFactor <- 10/data.tree::Aggregate(df, nodeSize, stats::median)
    t <- data.tree::Traverse(df, hierarchy_attribute)
    # traverse down the tree and compute the size of each node
    data.tree::Do(t, function(x) {
      # x$SizeOfNode <- data.tree::Aggregate(x, nodeSize, aggFun)
      # scale node growth to area rather than radius and round
      # x$SizeOfNode <- round(sqrt(x$SizeOfNode*scaleFactor)*pi, 2)
      x$SizeOfNode <- eval(parse(text = paste("x$", nodeSize, sep = "")))
      # scale node growth to area rather than radius and round
      x$SizeOfNode <- round(sqrt(x$SizeOfNode)*pi, 2)
    })
    # update left margin based on new root size
    jsonFields <- c(jsonFields, "SizeOfNode")
  }

  # keep only the JSON fields that are necessary
  if(is.null(jsonFields)) jsonFields <- NA
  data <- data.tree::ToListExplicit(df, unname = TRUE, keepOnly = jsonFields)
  data <- rename.node(data)

  # pass the data and options using 'x'
  x <- list(
    data = data,
    options = options
  )

  # create the widget
  htmlwidgets::createWidget(
    "collapsibleTree", x, width = width, height = height,
    htmlwidgets::sizingPolicy(viewer.padding = 0)
  )
}

# PO helper function left for the record
pieProportions <- function(node) {
  return(c(node$Cell.type, sapply(node$children, pieProportions)))
}

# color legend
#plot_color_map <- function(colors, ctypes){
#	pdf('color_map.pdf', width=2.5, height=10)
#	par(mar=c(1.1, 10, 1.1, 1.1))
#		image(y=1:length(colors), x=1, t(as.matrix(1:length(colors))), col= colors, axes=F, ylab='', xlab=''); axis(2, at=1:length(colors), las=2, labels = ctypes, cex.axis=0.5)
#	dev.off()
#}
# helper function to sort children by attribute, casting the given value to numeric
SortNumeric = function (node, attribute, ..., decreasing = FALSE, recursive = TRUE)
{
    if (node$isLeaf)
        return()
    ChildL <- sapply(node$children, function(x) GetAttribute(x,
        attribute, ...))
    names(ChildL) <- names(node$children)
    node$children <- node$children[order(as.numeric(ChildL), decreasing = decreasing,
        na.last = TRUE)]
    if (recursive)
        for (child in node$children) SortNumeric(child, attribute, ...,
            decreasing = decreasing, recursive = recursive)
    invisible(node)
}
