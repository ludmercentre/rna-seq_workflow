# load required libraries
library(dplyr)
library(gprofiler2)
library(plotly)

run_gProfiler <- function(bg_genes_list_filepath, results_folder_name, comparison_name) {

    result_folder_path = ""

    # import background gene list
    bg_genes_list_df <- read.csv(bg_genes_list_filepath)

    # import significant differentially expressed (DE) gene list
    sg_genes_list_df <- read.csv(results_folder_name, comparison_name)

    # Sort gene list by -1log10(p-value)*sign(logFC)
    ## Generate log10_pvalue_signed column
    sg_genes_list_df$log10_pvalue_signed = -1*log10(sg_genes_list_df[,'P.Value'])*sign(sg_genes_list_df[,'logFC'])

    ## Sort by absolute values of log10_pvalue_signed
    sg_genes_list_df <- sg_genes_list_df[order(abs(sg_genes_list_df$log10_pvalue_signed), decreasing=T),]


    # Keep significant genes only (filter by adjusted p-value)
    sg_genes_list_df <- sg_genes_list_df %>% filter(adj.P.Val <= 0.05)


    # Separate between UP-regulated and DOWN-regulated gene lists (using log of Fold Change values)
    UP_sg_genes_list_df <- sg_genes_list_df %>% filter(logFC > 0)

    DOWN_sg_genes_list_df <- sg_genes_list_df %>% filter(logFC < 0)


    # Vectorize lists
    bg_genes_list_l = as.vector(bg_genes_list_df$ensembl_gene_id)

    UP_sg_genes_list_l = as.vector(UP_sg_genes_list_df$ensembl_gene_id)
    DOWN_sg_genes_list_l = as.vector(DOWN_sg_genes_list_df$ensembl_gene_id)


    # Check if a list is empty and replace it wit a vector containing only "0" if that is the case
    if (length(UP_sg_genes_list_l) == 0) {
        UP_sg_genes_list_l <- c("0")
    }

    if (length(DOWN_sg_genes_list_l) == 0) {
        DOWN_sg_genes_list_l <- c("0")
    }


    # Output path
    output_folder_path = paste0(results_folder_name, comparison_name)


    # If both UP and DOWN
    ## Generate g:GOSt object
    multi_gostres1 <- gost(query = list("UP" = UP_sg_genes_list_l,
                                        "DOWN" = DOWN_sg_genes_list_l),
                            evcodes = TRUE, 
                            multi_query = FALSE,
                            organism = "mmusculus", 
                            ordered_query = TRUE, 
                            domain_scope = "custom", 
                            custom_bg = bg_genes_list_l)


    # Plot interactive manhattan like plot of UP and DOWN regulated genes (gostplot())
    p <- gostplot(multi_gostres1, capped = TRUE, interactive = TRUE)
    
    htmlwidgets::saveWidget(as_widget(p), paste0(result_folder_path, output_folder_path, "/gostplot.html"))

    # Generate GEM file for use with Cytoscape
    gem <- multi_gostres1$result[,c("query", "term_id", "term_name", "p_value", "intersection")]
    colnames(gem) <- c("query", "GO.ID", "Description", "p.Val", "Genes")
    gem$FDR <- gem$p.Val

    ## From: https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files
    ## The Generic Enrichment Results file needs:
    ### gene-set ID (must match the gene-set ID in the GMT file)
    ### gene-set name or description
    ### p-value
    ### FDR correction value
    ### Phenotype: +1 or -1, to identify enrichment in up- and down-regulation, or, more in general, in either of the two phenotypes being compared in the two-class analysis
    #### +1 maps to red
    #### -1 maps to blue
    ### gene list separated by commas

    ## Assign phenotype based on "UP" or "DOWN" regulated values
    gem$Phenotype <- ifelse(gem$query == "UP", "+1", "-1")

    ## Write single file with both queries
    write.table(data.frame(gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]),
                file = paste0(result_folder_path, output_folder_path, "/up_down_gem.txt"),
                sep = "\t", quote = F, row.names = F)

    ## Write separate files for queries
    # gem %>% group_by(query) %>%
    #   group_walk(~
    #     write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]), 
    #                 file = paste0("gProfiler_", unique(.y$query), "_gem.txt"),
    #                 sep = "\t", quote = F, row.names = F))


    # Generate table with the result statistics PDF image (publish_gosttable()) and CSV file:
    ## Get UP and DOWN regulated results:
    UP_multi_gostres1 = multi_gostres1$result[ which(multi_gostres1$result$query =='UP'),]
    DOWN_multi_gostres1 = multi_gostres1$result[ which(multi_gostres1$result$query =='DOWN'),]

    ## Order by p_value (increasing):
    UP_multi_gostres1 = UP_multi_gostres1[order(UP_multi_gostres1$p_value),]
    DOWN_multi_gostres1 = DOWN_multi_gostres1[order(DOWN_multi_gostres1$p_value),]

    ## Generate CSV files (remove last column "parents"):
    write.csv(UP_multi_gostres1[1:ncol(UP_multi_gostres1)-1], paste0(result_folder_path, output_folder_path, "/GO_analysis_UP_table.csv"), row.names=F)
    write.csv(DOWN_multi_gostres1[1:ncol(DOWN_multi_gostres1)-1], paste0(result_folder_path, output_folder_path, "/GO_analysis_DOWN_table.csv"), row.names=F)

    ## Generate PDF files:
    if (nrow(UP_multi_gostres1) > 100) {
        UP_multi_gostres1_forPDF = UP_multi_gostres1[1:100,]
    } else {
        UP_multi_gostres1_forPDF = UP_multi_gostres1
    }

    if (nrow(DOWN_multi_gostres1) > 100) {
        DOWN_multi_gostres1_forPDF = DOWN_multi_gostres1[1:100,]
    } else {
        DOWN_multi_gostres1_forPDF = DOWN_multi_gostres1
    }

    if (nrow(UP_multi_gostres1_forPDF) > 0) {
        publish_gosttable(UP_multi_gostres1_forPDF,
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_sizes"),
                        filename = paste0(result_folder_path, output_folder_path, "/GO_analysis_UP_table.pdf"))
    }

    if (nrow(DOWN_multi_gostres1_forPDF) > 0) {
        publish_gosttable(DOWN_multi_gostres1_forPDF,
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_sizes"),
                        filename = paste0(result_folder_path, output_folder_path, "/GO_analysis_DOWN_table.pdf"))
    }
}
