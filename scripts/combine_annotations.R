require(AnnotationForge)
require(tidyr)
require(dplyr)
require(foreach)
require(doMC)
require(gdata)
require(vroom)
require(splitstackshape)
registerDoMC(12)  # Allow to use several processors in parallel. 

save.image(file = "make_annotation_table.rda")

###############################################################
###### Combine the annotations from Blast2GO and BioMart ######
###############################################################

# Blast2GO

##Reads the Blast2GO file with the annotation of all genes
blast2GO_complete <- vroom(snakemake@input[[1]])

## Changes the Enzyme codes (E.C) for the same pattern used by BioMart
blast2GO_complete$Enzyme.Code <- gsub("EC:", "", blast2GO_complete$Enzyme.Code)

##Filters the annotations columns which are not informative.
blast2GO <- blast2GO_complete[, -c(3:15, 18)]

# BioMart

## Reads the BioMart file with the annotation of all genes
biomart_complete <- vroom(snakemake@input[[2]])

## Filters the annotations columns which are not informative.
biomart_data <- biomart_complete[, -c(9:13, 16, 17, 30:33)]

## Remove Nas
biomart_data <- sapply(biomart_data, as.character)
biomart_data[is.na(biomart_data)] <- ""
biomart_data <- as.data.frame(biomart_data)

# Selects the compatible annotation between the both sources, BioMart and Blast2GO and merge the datas in an unique data.frame based in the genes names.
blast2GO <- blast2GO[, c("Sequence.Name",
                         "Sequence.Description",
                         "Annotation.GO.ID",
                         "Annotation.GO.Term",
                         "Enzyme.Code",
                         "Enzyme.Name")]

colnames(blast2GO) <- c("gene_name",
                        "b2g_gene_description",
                        "b2g_go_id",
                        "b2g_go_term",
                        "b2g_kegg_enzyme_id",
                        "b2g_kegg_enzyme_desc")

biomart_data <- biomart_data[, c("gene_name1",
                                 "gene_description",
                                 "go_id",
                                 "go_desc",
                                 "kegg_enzyme_id",
                                 "kegg_enzyme_desc")]

colnames(biomart_data) <- c("gene_name",
                            "bioM_gene_description",
                            "bioM_go_id",
                            "bioM_go_term",
                            "bioM_kegg_enzyme_id",
                            "bioM_kegg_enzyme_desc")

all_annot <- merge(biomart_data, blast2GO, by = "gene_name")

## Compares the two annotations. The final results are the union between the annotations. Additionaly, the terms which was exclusive for each source are insert in new columns, for manual inspection.

# vector with the columns for comparison.
feature <- c("gene_description",
             "go_id",
             "go_term",
             "kegg_enzyme_id",
             "kegg_enzyme_desc")

annot_final <- data.frame(gene_name = all_annot[, "gene_name"])

count <- as.numeric(0)
for (f in feature){
  
  combine_data_all <- data.frame()
  
  # Extract the columns to comparison using the feature vector. In each loop, one annotation is compared.
  teste_set <- all_annot[, c(1, grep(f, colnames(all_annot)))]
  
  combine_data_all_2 <- foreach(i = 1:nrow(all_annot), .combine = "rbind") %dopar% {
    
    # Determines the terms for each source, Blast2GO and BioMart.
    gene_name <- paste(teste_set[i, 1])
    bioM_set_names <- strsplit(as.character(teste_set[i, grep("bioM", colnames(teste_set))]),
                               ",",
                               fixed = TRUE)[[1]]
    
    b2g_set_names <- strsplit(as.character(teste_set[i, grep("b2g", colnames(teste_set))]),
                              ";",
                              fixed = TRUE)[[1]]
    
    # Compares the terms and build the data.frame.
    comb_merge <- paste(unique(c(bioM_set_names, b2g_set_names)),
                        collapse = ";")
    bioM_uniq <- paste(unique(bioM_set_names[!bioM_set_names %in% b2g_set_names]),
                       collapse = ";")
    b2g_uniq <- paste(unique(b2g_set_names[!b2g_set_names %in% bioM_set_names]),
                      collapse = ";")
    common_terms <- paste(unique(b2g_set_names[b2g_set_names %in% bioM_set_names]),
                          collapse = ";")
    
    combine_data <- data.frame(gene_name,
                               comb_merge,
                               bioM_uniq,
                               b2g_uniq,
                               common_terms)
    
    colnames(combine_data) <- c("gene_name",
                                paste(f, "merged", sep = "_"),
                                paste(f, "unique_in_BioMart", sep = "_"),
                                paste(f, "unique_in_blastGO", sep = "_"),
                                paste(f, "common_to_both", sep = "_"))
    
    combine_data_all <- rbind(combine_data_all, combine_data)
    
  }
  
  annot_final <- merge(annot_final, combine_data_all_2, by = "gene_name")
}

annot_final <- data.frame(lapply(annot_final, function(x) {gsub("NA ;", "", x)} ))
annot_final <- data.frame(lapply(annot_final, function(x) {gsub("; ", ";", x)} ))
annot_final <- data.frame(lapply(annot_final, function(x) {gsub("NA", "", x)} ))

# Insert the annotation unique generated in the blast2GO.
blast2GO_unique <- blast2GO_complete[, c("Sequence.Name",
                                         "InterPro.Accession",
                                         "InterPro.Type",
                                         "InterPro.Name",
                                         "InterPro.Signatures",
                                         "InterPro.GO.ID",
                                         "InterPro.GO.Term",
                                         "InterPro.GO.Category")]

colnames(blast2GO_unique) <- c("gene_name",
                               "InterPro.Accession",
                               "InterPro.Type",
                               "InterPro.Name",
                               "InterPro.Signatures",
                               "InterPro.GO.ID",
                               "InterPro.GO.Term",
                               "InterPro.GO.Category")

annot_final <- merge(annot_final, blast2GO_unique, by = "gene_name")

annot_final <- annot_final[, -grep("unique", colnames(annot_final))]
annot_final <- annot_final[, -grep("common", colnames(annot_final))]

# ## Prepares the table with one term per line
# test <- teste[, c("gene_name", "kegg_enzyme_id_merged")] %>%
#   separate_rows("kegg_enzyme_id_merged", sep = ";")
# 
# test2 <- teste[, c("gene_name", "kegg_enzyme_desc_merged")] %>%
#   separate_rows("kegg_enzyme_desc_merged", sep = ";")
# 
# test_final <- merge(test, test2, by = "gene_name")
# 
# teste_unested <- separate_rows(teste, "gene_description_merged", sep = ";") %>%
#   separate_rows("go_id_merged", sep = ";")  %>%
#   separate_rows("go_term_merged", sep = ";")  %>%
#   separate_rows("kegg_enzyme_id_merged", sep = ";")  %>%
#   separate_rows("kegg_enzyme_desc_merged", sep = ";")  %>%
#   separate_rows("InterPro.Accession", sep = ";")  %>%      
#   separate_rows("InterPro.Type", sep = ";")  %>%           
#   separate_rows("InterPro.Name", sep = ";")  %>%          
#   separate_rows("InterPro.Signatures", sep = ";")  %>%     
#   separate_rows("InterPro.GO.ID", sep = ";")  %>%                    
#   separate_rows("InterPro.GO.Term", sep = ";")  %>%       
#   separate_rows("InterPro.GO.Category", sep = ";")
# 
# teste_unested <- separate_rows(teste, "kegg_enzyme_id_merged", sep = ";")  %>%
#   separate_rows("kegg_enzyme_desc_merged", sep = ";") 
# 
# write.table(annot_final,
#             "anotation_of_all_e.grandis_genes.tst",
#             col.names = F,
#             row.names = F,
#             sep = "\t",
#             quote = F)

# read the file with the distances between the tested sites and the closest gene
distance <- vroom(snakemake@input[[3]], delim = "\t", col_names = F)

# keep only the closest element of each methylation site
distance2 <- distance %>% 
  dplyr::group_by(X4) %>%
  dplyr::slice(which.min(abs(X16)))

distance2 <- distance2[, c(4, 2, 3, 6, 16, 15, 7, 10, 11, 13)]

# Select the columns of interest and split the name of genes.

distance2 <- cSplit(distance2, "X15", "=") %>% 
  dplyr::select(-"X15_1", - "X15_2")

colnames(distance2) <- c("locus",
                         "locus_start",
                         "locus_end",
                         "locus_strand",
                         "distance_to_gene(kb)",
                         "chromosome",
                         "gene_start",
                         "gene_end",
                         "gene_strand",
                         "gene_name")

distance2 <- as.data.frame(distance2)

# For each loci, merge the distance to the closest gene with the annotation.
annot_table <- merge(distance2, annot_final, by = "gene_name")

annot_table$position_to_gene <- ifelse(annot_table$`distance_to_gene(kb)` < 0, "upstream",
                                       ifelse(annot_table$`distance_to_gene(kb)` > 0, "downstram", "overlapping"))

annot_table$`distance_to_gene(kb)` <- abs(annot_table$`distance_to_gene(kb)`)

annot_table <- annot_table %>%
  dplyr::select("chromosome",
                "gene_name",
                "locus",
                "distance_to_gene(kb)",
                "position_to_gene",
                "gene_description_merged",
                "go_id_merged",
                "go_term_merged",
                "kegg_enzyme_id_merged",
                "kegg_enzyme_desc_merged",
                "InterPro.Name",
                "InterPro.GO.ID",
                "InterPro.GO.Term")

# Writes the annotation table with the agregated values 
write.table(annot_table,
            snakemake@output[[2]],
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

save.image("make_annotation_table_2_0.rda")

####### >>>>> Join it with the information of if the MSD-site in UTR or Exon

# # Reads the marks that overlap a gene
# utr_exon <- vroom("marks_in_exons_or_utrs.tst", col_names = F) %>%
#   dplyr::select(X4, X6, X9, X15) %>%
#   dplyr::rename("locus" = X4,
#                 "MDS strand" = X6,
#                 "gene feature" = X9,
#                 "gene_name_iso" = X15)
# 
# # Remove the duplications caused by the multiple isoforms
# utr_exon <- dplyr::select(utr_exon, -gene_name_iso) %>%
#   dplyr::distinct()
# 
# # gene_name <- data.frame(do.call("rbind",
# #                                 strsplit(utr_exon$gene_name_iso,
# #                                          split = "=",
# #                                          fixed = TRUE)))
# # gene_name <- as.character(gene_name$X3)
# # 
# # gene_name <- data.frame(do.call("rbind",
# #                                 strsplit(gene_name,
# #                                          split = ".v2.0",
# #                                          fixed = TRUE)))
# # utr_exon$gene_name <- gene_name$X1
# 
# # Select the set of genes that have a MSD-site within it
# #annot_table_genes <- annot_table[annot_table$position_to_gene == "overlapping", ]
# 
# annot_table_2 <- merge(annot_table,
#                        utr_exon,
#                        by = "locus",
#                        all = T)
# 
# annot_table_2 <- dplyr::select(annot_table_2,
#                                "locus",
#                                "chromosome",
#                                "distance_to_gene(kb)",
#                                "gene_name",
#                                "position_to_gene",
#                                "gene feature",
#                                "gene_description_merged",
#                                "go_id_merged",
#                                "go_term_merged",
#                                "kegg_enzyme_id_merged",
#                                "kegg_enzyme_desc_merged",
#                                "InterPro.Name",
#                                "InterPro.GO.ID",
#                                "InterPro.GO.Term")
# 
# ## Identify the introns (inside genes, but not in exon or UTR)
# for(i in 1:nrow(annot_table_2)) {
#   
#   if(annot_table_2[i, 3] == 0 && is.na(annot_table_2[i, 6])) {
#     
#     annot_table_2[i, 6] <- "intron"
#     
#   }
#   
# }
# 
# 
# 
# # Writes the table with the combined annotation.
# write.table(annot_table,
#             snakemake@output[[2]],
#             col.names = T,
#             row.names = F,
#             sep = "\t",
#             quote = F)
# 
# save.image("distance_annot_final_2.0_final.rda")