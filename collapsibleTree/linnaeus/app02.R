library(shiny)
library(collapsibleTree)
library(data.tree)

#load('/Users/polivar/src/linnaeus-scripts/collapsibleTree/sand/C2_correct_tree_wcells.Robj')
#load('/Users/polivar/Downloads/bs/Z2_Ltree.Robj')
#ttt = Clone(LINNAEUS.cell.tree$nd2)
load('src/linnaeus-scripts/collapsibleTree/sand/Z2_Ltree_pie.Robj')

ctypes = linnaeus.colors_larva$Cell.type
ct_colors = linnaeus.colors_larva$color

#names(ct_colors) = ctypes

orit = Clone(LINNAEUS.pie)
get_pieNode(orit, ctypes=ctypes)
ttt = Clone(orit)

namess = ttt$Get(function(x) if(x$isScar) x$name)
namess = namess[!is.na(namess)]
names(namess) = namess


do_barplot = function(foc_pie, ct_colors, zeros=F){
 par(mar=c(15, 3, 0.1, 0.1))
 #names(ct_colors) = ctypes
 if(!zeros){
  f_n = foc_pie > 0 
 }else{
  f_n = 1:length(foc_pie)
 }
 foc_pie = foc_pie[f_n]
 barplot(foc_pie, col = ct_colors[f_n], las = 2, space = 0.2)
 grid()
}

linkLength = 180
pieSummary = TRUE
first = TRUE
# Define server logic required to draw a collapsible tree diagram
server <- function(input, output) {
   output$plot <- renderCollapsibleTree({

	linkLength <<- ifelse(first, linkLength, input$linkLength)
	pieSummary <<- ifelse(first, pieSummary, input$pieSummary)
	if(first){first <<- FALSE}
	if(input$root != orit$name){
		collapsibleTree(FindNode(ttt, input$root), collapsed=F, inputId = "node", pieNode=T, pieSummary=pieSummary, hide_scars=T, linkLength = linkLength, do_collapse= FALSE)
	} else{
     	#collapsibleTree(ttt, collapsed=F, inputId = "node", pieNode=T, pieSummary=pieSummary, hide_scars=F, linkLength = linkLength)
	# TODO hide_scars must be set to true as current trees lack the scar field, if F then barplots break as there is no name for a node
	print(pieSummary)
     	collapsibleTree(ttt, collapsed=F, inputId = "node", pieNode=T, pieSummary=pieSummary, hide_scars=T, linkLength = linkLength, do_collapse= FALSE)}
   })
	observe(
		if(!is.null(input$node)){

			clicked_node = input$node[length(input$node)]
			if(length(input$node) == 0){
   				foc_pie = ttt$pieNode
   				output$str <- renderPrint(cat(ttt$name))
			}else{
	   			foc_node = FindNode(ttt, clicked_node)
	   			foc_pie = foc_node$pieNode
   				output$str <- renderPrint(cat(foc_node$name))
			}
	   		output$barplot_ctypes <- renderPlot(do_barplot(foc_pie, ct_colors = ct_colors))
			output$treeheight <- renderText({input$treeheight})
	})
}

get_clicked_node_name = function(dt, input){
	clicked_node = input$node[length(input$node)]
				if(length(input$node) == 0){
					return(dt$name)
				}else{
		   			foc_node = FindNode(dt, clicked_node)
					return(foc_node$name)
				}
}


# Define UI for application that draws a collapsible tree
ui <- fluidPage(

   # Application title
   titlePanel("Linnaeus"),

   # Sidebar with a select input for the root node
   sidebarLayout(position ='right',
      sidebarPanel(
         selectInput("root", "Select a node to render a sub-tree", namess),
         tags$p("Scar from the most recently clicked node:"),
         verbatimTextOutput("str"),
	
	checkboxInput(inputId = "pieSummary",
	      label = strong("Hide single cells"),  value = TRUE),

	sliderInput(inputId = "linkLength",
        	label = "Link Length",
	        min = 20, max = 500, value = 100, step = 5),

	sliderInput(inputId = "treeheight",
        	label = "Tree panel height",
	        min = 500, max = 2000, value = 1000, step = 5)

	#conditionalPanel(condition = TRUE,
	),

      # Show a tree diagram with the selected root node
      mainPanel(
        collapsibleTreeOutput("plot", height = '450px'),
  	plotOutput(outputId = "barplot_ctypes", height = "500px")
      )
   )
)

# Run the application
shinyApp(ui = ui, server = server)
