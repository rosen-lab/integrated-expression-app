library(shiny)
library(shinyjs)
# shiny::shinyAppDir(getwd(),options=list(port=3880,launch.browser=FALSE))

# Contact info footer
contact <- div(
	class = "text-center",
	hr(),
	"This app was made for ",
	a(href = "http://www.evanrosenlab.net/", "Evan Rosen's lab"),
	"at BIDMC. The source code can be found on ",
	a(href = "https://github.com/rosen-lab/motif_diff.R", "github")
)

# Define app tab
app <- tabPanel("App",
	useShinyjs(),
	
	# Some CSS styling
	tags$head(
		tags$style(HTML("
			.x-large-text {
				font-size: x-large;
			}
			.large-text {
				font-size: large;
			}"
			))
		),
	
	# Main controls panel
	tabsetPanel(type='tabs',
		# "Gene-to-Peak Integration"
		tabPanel("Gene->Peak Integration",
			div(class = "well",
				# Gene data upload
				fluidRow(
					class = "x-large-text",
					column(4,"Gene Targets")
					),
				fluidRow(
					column(4,
						fileInput('gene.targets.file',
							label=NULL,
							accept=c(".txt",".csv",".tsv")
							)
						),
					column(1,
						actionButton(
							"gene.targets.show",
							"Show",
							class="btn btn-primary btn-md btn-block"
							)
						)
					),
				# Gene data filters
				span(class="x-large-text","Gene Target Filters"),
				fluidRow(
					column(
						4,
						selectizeInput("gene.targets.filter.fdr",label="Maximum FDR",
# 								choices=c("Choose (or enter) any number of values:"="",0.05,0.25),
							choices=c("Choose or enter any value [0-1]:"="",0.05,0.25),
							selected=c(0.05),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						),
					column(
						4,
						selectizeInput("gene.targets.filter.fc",label="Minimum Fold-Change (log2)",
# 								choices=c("Choose (or enter) any number of values:"="",0.5,1),
							choices=c("Choose or enter any value [0-):"="",0.5,1),
							selected=c(1),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						)
					),
				# ChIP data upload
				fluidRow(
					class = "x-large-text",
					column(4,"ChIP Peak Table")
					),
				fluidRow(
					column(4,
						fileInput('peak.table.file',
							label=NULL,
							accept=c(".txt",".csv",".tsv")
							)
						),
					column(1,
						actionButton(
							"peak.table.show",
							"Show",
							class="btn btn-primary btn-md btn-block"
							)
						)
					),
				# ChIP data filters
				span(class = "x-large-text", "ChIP Peak Filters"),
				fluidRow(
					column(
						4,
						selectizeInput("peak.table.filter.fdr",label="Maximum FDR",
# 								choices=c("Choose (or enter) any number of values:"="",0.05,0.25),
							choices=c("Choose or enter any value [0-1]:"="",0.05,0.25),
							selected=c(0.05),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						),
					column(
						4,
						selectizeInput("peak.table.filter.fc",label="Minimum Fold-Change (log2)",
# 								choices=c("Choose (or enter) any number of values:"="",0.5,1),
							choices=c("Choose or enter any value [0-):"="",0.5,1),
							selected=c(1),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						),
					column(
						4,
						selectizeInput("peak.table.filter.fdr",label="Minimum CPM (log2)",
# 								choices=c("Choose (or enter) any number of values:"="",10,100,1000),
							choices=c("Choose or enter any value [0-):"="",0.5,1),
							selected=NULL,
							multiple=FALSE,
							options=list(create=TRUE)
							)
						)
					),
				# Integration filters
				span(class = "x-large-text", "Additional Integration Filters"),
				fluidRow(
					column(
						4,
						selectizeInput("g2p.filter.dist",label="Maximum Distance btw. TSS and Peak (kb)",
							choices=c("Choose any number of values:"="",10,100,250,1000,2500),
							selected=c(10),
							multiple=TRUE,
							options=list(create=TRUE)
							)
						)
					),
				# Integration button
				fluidRow(
					column(
						4,
						offset = 4,
						actionButton(
							"g2p.run",
							"Integrate!",
							class="btn btn-primary btn-lg btn-block"
							)
						)
					)
				)
			),
		# "Peak-to-Gene Integration"
		tabPanel("Peak->Gene Integration",
			div(class = "well",
				# Peak data upload
				fluidRow(
					class = "x-large-text",
					column(4,"Peak Targets")
					),
				fluidRow(
					column(4,
						fileInput('peak.targets.file',
							label=NULL,
							accept=c(".txt",".csv",".tsv")
							)
						),
					column(1,
						actionButton(
							"peak.targets.show",
							"Show",
							class="btn btn-primary btn-md btn-block"
							)
						)
					),
				# Peak data filters
				span(class="x-large-text","Gene Target Filters"),
				fluidRow(
					column(
						4,
						selectizeInput("peak.targets.filter.fdr",label="Maximum FDR",
# 								choices=c("Choose (or enter) any number of values:"="",0.05,0.25),
							choices=c("Choose or enter any value [0-1]:"="",0.05,0.25),
							selected=c(0.05),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						),
					column(
						4,
						selectizeInput("peak.targets.filter.fc",label="Minimum Fold-Change (log2)",
# 								choices=c("Choose (or enter) any number of values:"="",0.5,1),
							choices=c("Choose or enter any value [0-):"="",0.5,1),
							selected=c(1),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						)
					),
				# Gene data upload
				fluidRow(
					class = "x-large-text",
					column(4,"Gene Table")
					),
				fluidRow(
					column(4,
						fileInput('gene.table.file',
							label=NULL,
							accept=c(".txt",".csv",".tsv")
							)
						),
					column(1,
						actionButton(
							"gene.table.show",
							"Show",
							class="btn btn-primary btn-md btn-block"
							)
						)
					),
				# Gene data filters
				span(class = "x-large-text", "Gene Filters"),
				fluidRow(
					column(
						4,
						selectizeInput("gene.table.filter.fdr",label="Maximum FDR",
# 								choices=c("Choose (or enter) any number of values:"="",0.05,0.25),
							choices=c("Choose or enter any value [0-1]:"="",0.05,0.25),
							selected=c(0.05),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						),
					column(
						4,
						selectizeInput("gene.table.filter.fc",label="Minimum Fold-Change (log2)",
# 								choices=c("Choose (or enter) any number of values:"="",0.5,1),
							choices=c("Choose or enter any value [0-):"="",0.5,1),
							selected=c(1),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						),
					column(
						4,
						selectizeInput("gene.table.filter.cpm",label="Minimum CPM (log2)",
# 								choices=c("Choose (or enter) any number of values:"="",10,100,1000),
							choices=c("Choose or enter any value [0-):"="",0.5,1),
							selected=NULL,
							multiple=FALSE,
							options=list(create=TRUE)
							)
						)
					),
				# Integration filters
				span(class = "x-large-text", "Additional Integration Filters"),
				fluidRow(
					column(
						4,
						selectizeInput("p2g.filter.dist",label="Maximum Distance btw. TSS and Peak (kb)",
							choices=c("Choose any number of values:"="",10,100,250,1000,2500),
							selected=c(10),
							multiple=FALSE,
							options=list(create=TRUE)
							)
						)
					),
				fluidRow(
					column(
						4,
						checkboxInput("p2g.filter.matches",label="Singleton matching?",
							value=FALSE,
							)
						)
					),
				# Integration button
				fluidRow(
					column(
						4,
						offset = 4,
						actionButton(
							"p2g.run",
							"Integrate!",
							class="btn btn-primary btn-lg btn-block"
							)
						)
					)
				)
			)
		),
	br(),
	downloadButton("download","Download"),
	br(),
	dataTableOutput("table.integrated"),
	
	# Add contact info
	br(),
	contact
	)

about <- tabPanel("About",HTML("<h1>TODO</h1>"))

# Define UI as a navbar with app and about tabs
shinyUI(navbarPage("Integrated Data Motif Enchrichment App", app, about))
