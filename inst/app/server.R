#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#remotes::install_github("rstudio/bootstraplib")
library(vcfR)
library(plotly)
library(biomaRt)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)

# Define server 
options(shiny.maxRequestSize=500*1024^2)
shinyServer(function(input, output) {
    # set up biomart
    ensembl <- useMart("ensembl", host="grch37.ensembl.org",  dataset="hsapiens_gene_ensembl")
    snp_mart <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_snp")

    # read in vcf
    vcfInput <- reactive({
        req(input$vcf)
        data <- read.vcfR(input$vcf$datapath, verbose = F)
        })
    
    # print how many variants in vcf
    output$success <- renderText({
        n_vars <- nrow(vcfInput())
        print(paste("File with", n_vars, "variants successfully processed!"))
        }) 
    
    # define gene of interest when action button is pressed
    goi <- eventReactive(input$goi_go, {
      input$goi
      })
    
    # get location of gene of interest
    gene_pos <- reactive({
      req(input$goi)
      df <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position'), 
                                               filters = c('hgnc_symbol'),
                                               values= goi(),
                                               mart=ensembl)
      })
    
    # get variants in that region 
    my_vars <- reactive({
      req(gene_pos(), vcfInput())
      subset(vcfInput(), as.numeric(str_extract(vcfInput()@fix[,'CHROM'], '[0-9]+')) == gene_pos()$chromosome_name & 
                         as.numeric(vcfInput()@fix[,'POS']) >= gene_pos()$start_position &
                         as.numeric(vcfInput()@fix[,'POS']) <= gene_pos()$end_position)
      })
    
    # create gene of interest df for plotting
    goi_df <- reactive({
      req(c(gene_pos(), my_vars()))
      data.frame(position = as.numeric(my_vars()@fix[, 'POS']), 
                 variant =my_vars()@fix[,'ALT'],
                 num_alleles = ifelse(extract.gt(my_vars())[,1]== '1/1' | extract.gt(my_vars())[,1]== '1|1', 2, 1),
                 id = as.character(my_vars()@fix[,'ID']), stringsAsFactors = F)
      })
    
    # output gene of interest plot
    output$goi_plot <- renderPlotly({
      goi_plot(goi_df(), gene_pos(), goi())
    })
    
    
    dt <- reactive({
      if(is.na(unique(goi_df()$id)[1])){
        dt <- goi_df()}
      else{
        id <- paste0('<a href="https://www.ncbi.nlm.nih.gov/snp/', goi_df()$id, '">', goi_df()$id, '</a>')
        sift <- getBM(attributes = c('refsnp_id', 'sift_prediction', 'consequence_allele_string'), 
                      filters = 'snp_filter', 
                      values = goi_df()$id, 
                      mart = snp_mart) %>% 
          filter(sift_prediction != "") %>%
          mutate(id = refsnp_id) %>% 
          inner_join(goi_df(), by = 'id') %>%
          filter(str_extract(consequence_allele_string, '([^/]+$)') == variant) %>%
          `rownames<-`(.[,1])
        predicted_consequence <- sift[goi_df()$id, 'sift_prediction']
        dt <- cbind(goi_df()[,c(1:3)], id, predicted_consequence)}
      dt
    })
    
    # output gene of interest df
    output$goi_table <- renderDataTable({
      dt()
    }, escape = c('position', 'variant', 'num_alleles'))
    
    
    # Downloadable csv of rare vars
    output$download_goi <- downloadHandler(
      filename = paste0(input$goi, '_variants.tsv'),
      content = function(file) {
        chr <- rep(gene_pos()$chromosome_name, nrow(dt()))
        goi_out <- cbind(chr, dt())
        goi_out$id <- str_extract(goi_out$id, 'rs[0-9]+')
        write.table(goi_out, file, row.names = FALSE, sep = '\t', quote = F)
      }
    )
    
    # output 1k genomes pca plot
    output$genomes_1k_pca <-  renderPlotly({
      req(input$pca_go)
      genomes_1k_pca(vcfInput())
    })
    
    # compute overlap with ukbb
    overlap_df <- reactive({
      req(input$rare_var_go)
      get_ukbb_overlap(vcfInput())
    })
    
    # define rare var df
    rare_var_df <- reactive({
      subset(overlap_df(), overlap_df()$percent_people < input$p_threshold)
    })
    
    # print rare var df
    output$rare_var_df <- renderDataTable({
      req(input$rare_var_go)
      output_df <- data.frame(rare_var_df(), id = paste0('<a href="https://www.ncbi.nlm.nih.gov/snp/', rare_var_df()$rsid, '">', rare_var_df()$rsid, '</a>'))
      output_df[,c('varid', 'id', 'percent_people', 'consequence', 'consequence_category')]
    }, options = list(pageLength = 10), escape = c(1,3:5))
    
    # Downloadable csv of rare vars
    output$download_rare_var <- downloadHandler(
      filename = 'rare_vars.tsv',
      content = function(file) {
        write.table(rare_var_df(), file, row.names = FALSE, sep = '\t', quote = F)
      }
    )
    
    # print title for rare var plot
    output$plot_title <- renderUI(h5(paste0('This is what sets you apart from ~', prettyNum(round(361194-(input$p_threshold/100*361194)),big.mark = ','), ' people!')))
    
    # print rare_var_plot
    output$rare_var_plot <- renderPlot({
      print_rare_vars(overlap_df(), as.numeric(input$p_threshold), pal =  pals[[as.numeric(input$pal)]])
      })
})

