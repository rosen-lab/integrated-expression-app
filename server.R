library(shiny)
library(shinyjs)
library(tibble)
library(dplyr)
library(readr)
# shiny::shinyAppDir(getwd(),options=list(port=3880,launch.browser=FALSE))

options(shiny.maxRequestSize=32*1024^2)

shinyServer(
	function(input,output,session)
	{
		# Gene targets
		gene.targets <- reactive(
			{
				file <- input$gene.targets.file
				# stop if the file is empty
				if (is.null(file))
				{
					return(NULL)
				}
				# fail if the file is not a table
				data <- try(read_tsv(file$datapath,header=TRUE),silent=TRUE)
				if (!is.data.frame(data))
				{
					showNotification(
						paste0("Unable to read table from input file. Please make sure the file is of an appropriate format to be read by the `readr::read_tsv` function"),
						duration=NULL,
						closeButton=TRUE,
						type='error'
						)
					return(NULL)
				}
				# fail if the data does not have the required column (headers)
				names <- colnames(data)
				for (name in c('ID','logFC','FDR','Direction'))
				{
					if (! name %in% names)
					{
						showNotification(
							paste0('Gene target file missing required column: ',name),
							duration=NULL,
							closeButton=TRUE,
							type='error'
							)
						return(NULL)
					}
				}
				# finally, return the table
				return(data)
			})
		disable('gene.targets.show')
		observeEvent(input$gene.targets.file,
			{
				enable('gene.targets.show')
			})
		observeEvent(input$gene.targets.show,
			{
				showModal(modalDialog(
					renderDataTable(gene.targets()),
					title='Gene Targets',
					footer=modalButton('Close'),
					size='l',
					easyClose=FALSE
					))
			})
		
		# Gene table file was selected
		gene.table <- reactive(
			{
				file <- input$gene.table.file
				# stop if the file is empty
				if (is.null(file))
				{
					return(NULL)
				}
				# fail if the file is not a table
				data <- try(read_tsv(file$datapath),silent=TRUE)
				if (!is.data.frame(data))
				{
					showNotification(
						paste0("Unable to read table from input file. Please make sure the file is of an appropriate format to be read by the `readr::read_tsv` function"),
						duration=NULL,
						closeButton=TRUE,
						type='error'
						)
					return(NULL)
				}
				# fail if the data does not have the required column (headers)
				names <- colnames(data)
				for (name in c('ID','logFC','logCPM','PValue','FDR'))
				{
					if (! name %in% names)
					{
						showNotification(
							paste0('Peak table missing required column: ',name),
							duration=NULL,
							closeButton=TRUE,
							type='error'
							)
						return(NULL)
					}
				}
				# enable run button (if this was the last file)
				if (! is.null(peak.targets))
				{
					enable('p2g.run')
				}
				# finally, return the table
				return(data)
			})
		disable('gene.table.show')
		observeEvent(input$gene.table.file,
			{
				enable('gene.table.show')
			})
		observeEvent(input$gene.table.show,
			{
				showModal(modalDialog(
					renderDataTable(gene.table()),
					title='Gene Table',
					footer=modalButton('Close'),
					size='l',
					easyClose=FALSE
					))
			})
		
		# Peak targets
		peak.targets <- reactive(
			{
				file <- input$peak.targets.file
				# stop if the file is empty
				if (is.null(file))
				{
					return(NULL)
				}
				# fail if the file is not a table
				data <- try(read_tsv(file$datapath),silent=TRUE)
				if (!is.data.frame(data))
				{
					showNotification(
						paste0("Unable to read table from input file. Please make sure the file is of an appropriate format to be read by the `readr::read_tsv` function"),
						duration=NULL,
						closeButton=TRUE,
						type='error'
						)
					return(NULL)
				}
				# fail if the data does not have the required column (headers)
				names <- colnames(data)
				for (name in c('ID','logFC','FDR','Direction','Neutral'))
				{
					if (! name %in% names)
					{
						showNotification(
							paste0('Peak target file missing required column: ',name),
							duration=NULL,
							closeButton=TRUE,
							type='error'
							)
						return(NULL)
					}
				}
				# finally, return the table
				return(data)
			})
		disable('peak.targets.show')
		observeEvent(input$peak.targets.file,
			{
				enable('peak.targets.show')
			})
		observeEvent(input$peak.targets.show,
			{
				showModal(modalDialog(
					renderDataTable(peak.targets()),
					title='ChIP Peak Targets',
					footer=modalButton('Close'),
					size='l',
					easyClose=FALSE
					))
			})
		
		# Peak table
		peak.table <- reactive(
			{
				file <- input$peak.table.file
				# stop if the file is empty
				if (is.null(file))
				{
					return(NULL)
				}
				# fail if the file is not a table
				data <- try(read_tsv(file$datapath,header=TRUE),silent=TRUE)
				if (!is.data.frame(data))
				{
					showNotification(
						paste0("Unable to read table from input file. Please make sure the file is of an appropriate format to be read by the `readr::read_tsv` function"),
						duration=NULL,
						closeButton=TRUE,
						type='error'
						)
					return(NULL)
				}
				# fail if the data does not have the required column (headers)
				names <- colnames(data)
				for (name in c('ID','logFC','logCPM','PValue','FDR'))
				{
					if (! name %in% names)
					{
						showNotification(
							paste0('Peak table missing required column: ',name),
							duration=NULL,
							closeButton=TRUE,
							type='error'
							)
						return(NULL)
					}
				}
				# finally, return the table
				return(data)
			})
		disable('peak.table.show')
		observeEvent(input$peak.table.file,
			{
				enable('peak.table.show')
			})
		observeEvent(input$peak.table.show,
			{
				showModal(modalDialog(
					renderDataTable(peak.table()),
					title='Peak Table',
					footer=modalButton('Close'),
					size='l',
					easyClose=FALSE
					))
			})
		
		# Gene-to-Peak integrate button
		observeEvent(input$g2p.run,
			{
				targets <- input$gene.targets.file
				table <- input$peak.table.file
				session <<- showNotification(
					paste0('Ima doosum integrtin of ',targets,'+',table,' now, kk?'),
					duration=10,
					closeButton=TRUE,
					type='message'
					)
			})
		
		# Peak-to-Gene integrate button was clicked
		integrated.table <- reactive({
			peaks2genes(peak.targets(),gene.table(),
				filt.fc.tgt=input$peak.targets.filter.fc,
				filt.fdr.tgt=input$peak.targets.filter.fdr
				)
			})
		observeEvent(input$p2g.run,
			{
				output$table.integrated <- renderDataTable(integrated.table())
				enable('download')
			})
		
		# Download button
		disable('download')
		output$download <- downloadHandler(
			filename = 'file.tsv',
			content = function(file)
				{
					write_tsv(integrated.table(),path=file)
				}
			)
	})

peaks2genes <- function (targets,table,
	filt.fdr.tgt=0.05,filt.fc.tgt=1.0,
	filt.fdr.tab=0.05,filt.fc.tab=1.0,filt.cpm.tab=0,
	filt.dist=10
	)
{
	# filter base targets (peaks)
	peaks <- dplyr::filter(targets,
		FDR <= filt.fdr.tgt,
		abs(logFC) >= filt.fc.tgt,
		Neutral == FALSE
		)
	
	# filter mapping table (genes)
	genes <- dplyr::filter(table,
		FDR <= filt.fdr.tab,
		abs(logFC) >= filt.fc.tab,
		logCPM >= filt.cpm.tab
		)
	genes <- add_column(genes,UCSCID=sub('([^|]*)\\|.*','\\1',genes$ID))
	
	# import local gene start site information
	tsses <- read_tsv('TSS.gtf',
		col_types='c--i----c',
		col_names=c('Chromosome','Position','Attributes')
		)
	tsses <- add_column(tsses,UCSCID=sub('gene_id "([^"]*)".*','\\1',tsses$Attributes))
	
	# match genes up with their transcription start sites
	genes <- dplyr::left_join(genes,tsses,by='UCSCID')
	
	# map genes to the peak targets
	final <- tibble(Peak=character(),
		Peak.Chromosome=character(),Peak.Start=integer(),Peak.Stop=integer(),
		Peak.logFC=double(),Peak.FDR=double(),
		Gene=character(),Gene.Position=integer(),
		Gene.logCPM=double(),Gene.logFC=double(),Gene.FDR=double()
		)
	withProgress(message = 'Integrating Peaks->Genes',value=0,
		{
			n.rows = nrow(peaks)
# 			n.rows = 100
			progress = 0.0;
			for (row in 1:n.rows)
			{
				peak <- peaks[row,]
				id <- peak$ID
				chr <- sub('.*\\|(chr[^:]+):.*','\\1',id)
				start <- as.numeric(sub('.*\\|chr[^:]+:(\\d+)-(\\d+).*','\\1',id))
				stop <- as.numeric(sub('.*\\|chr[^:]+:(\\d+)-(\\d+)','\\2',id))
				dir <- peak$Direction
				matched <- dplyr::filter(genes,
					between(Position,start-(1000 * filt.dist),stop+(1000 * filt.dist)),
					(dir == 'up') == (logFC >= 0)
					)
				if (nrow(matched) == 0)
				{
					next
				}
				final <- add_row(final,Peak=sub('(\\d+)\\|.*','\\1',id),
					Peak.Chromosome=chr,Peak.Start=start,Peak.Stop=stop,
					Peak.logFC=peak$logFC,Peak.FDR=peak$FDR,
					Gene=matched$UCSCID,Gene.Position=matched$Position,
					Gene.logCPM=matched$logCPM,Gene.logFC=matched$logFC,Gene.FDR=matched$FDR
					)
				incProgress((row / n.rows) - progress)
				progress = row / n.rows
			}
		})
	
	return (final)
}
