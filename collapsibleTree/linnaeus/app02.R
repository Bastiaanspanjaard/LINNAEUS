#library(shiny)
#library(collapsibleTree)
#library(data.tree)

#load('/Users/polivar/src/linnaeus-scripts/collapsibleTree/sand/C2_correct_tree_wcells.Robj')
#load('/Users/polivar/Downloads/bs/Z2_Ltree.Robj')
#ttt = Clone(LINNAEUS.cell.tree$nd2)
load('src/linnaeus-scripts/collapsibleTree/sand/Z2_Ltree_pie.Robj')
orit = Clone(LINNAEUS.pie)
namess = names(ttt$children)
names(namess) = namess

ctypes = linnaeus.defaults('larva')$Cell.type
ct_colors = linnaeus.defaults('larva')$color

#names(ct_colors) = ctypes

get_pieNode(orit, ctypes=ctypes)
ttt = Clone(orit)

first = TRUE
# Define UI for application that draws a collapsible tree
ui <- fluidPage(

   # Application title
   titlePanel("Linnaeus"),

   # Sidebar with a select input for the root node
   sidebarLayout(position ='right',
      sidebarPanel(
         selectInput("root", "[This doens't work yet] Select a node", namess),
         tags$p("Scar from the most recently clicked node:"),
         verbatimTextOutput("str"),
	
	checkboxInput(inputId = "collapsible",
	      label = strong("Show single cells"),  value = TRUE),

	checkboxInput(inputId = "branch",
	      label = strong("render branch"),  value = TRUE),

	sliderInput(inputId = "linkLength",
        	label = "Link Length",
	        min = 20, max = 200, value = 170, step = 10)

	#conditionalPanel(condition = TRUE,
	),

      # Show a tree diagram with the selected root node
      mainPanel(
        collapsibleTreeOutput("plot", height = '1000px'),
  	plotOutput(outputId = "bbb", height = "500px")
      )
   )
)

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

linkLength = 170
first = TRUE
# Define server logic required to draw a collapsible tree diagram
server <- function(input, output) {
   output$plot <- renderCollapsibleTree({
	if(input$branch & !is.null(input$node) & length(input$node) != 0){
		ttt = FindNode(ttt, input$node[length(input$node)])
 	}
	linkLength = ifelse(first, 170, input$linkLength)
	if(first){first = FALSE}
     collapsibleTree(ttt, collapsed=F, inputId = "node", pieNode=T, pieSummary=input$collapsible, hide_scars=T, linkLength = 170, height=500)
   })
	observe(
		if(!is.null(input$node)){

			clicked_node = input$node[length(input$node)]
			if(length(input$node) == 0){
   				foc_pie = ttt$pieNode
   				output$str <- renderPrint(cat(ttt$scar))
			}else{
	   			foc_node = FindNode(ttt, clicked_node)
	   			foc_pie = foc_node$pieNode
   				output$str <- renderPrint(cat(foc_node$scar))
			}
	   		output$bbb <- renderPlot(do_barplot(foc_pie, ct_colors = ct_colors))
	})
}

# Run the application
shinyApp(ui = ui, server = server)
