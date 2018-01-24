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
				# PO
				colors = NULL, 
				pieNode = FALSE, nCell.types=ifelse(!is.null(colors, length(colors), 1:60)), 
                                 zoomable = TRUE, width = NULL, height = NULL, ...) {

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

  # PO create colors
  if(is.null(colors)){
	colors = c('red', colorRampPalette(c('#fa9fb5','#f768a1','#ae017e','#4400d9', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8' ))(nCell.types))
  }

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
  
  # PO tmp code for dealing with ctype colour attribution
  ctypes = sort(unique(df$Get("Cell.type")))
  # if next line commented, we are keeping NAs
  #cytpes = cytpes[!is.na(ctypes)]


  # PO adding cell type counts 
  if(pieNode){
    t <- data.tree::Traverse(df, 'level')
    data.tree::Do(t, function(x) {
#message(x$levelName)
        xpieNode =  data.tree::Aggregate(x, 'Cell.type', function(j){
	#	message(collapse = " ", j)
		return(array(match(j, ctypes)))
		})
	
	x$pieNode = table(c(nCell.types, unlist(xpieNode))) -1 
	x$SizeOfNode = sum(x$pieNode)
#browser()
	})
    jsonFields <- c(jsonFields, "pieNode")
    jsonFields <- c(jsonFields, "SizeOfNode")
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
    options$margin$left <- options$margin$left + df$SizeOfNode - 10
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
