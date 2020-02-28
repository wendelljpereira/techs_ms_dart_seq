save.image("transposons_plots.rda")
require(ggplot2)
require(plyr)
require(VennDiagram)
require(gridExtra)
require(gdata)
require(scales)
require(splitstackshape)
require(SuperExactTest)

# LÃª o arquivo com a informaÃ§Ã£o dos transposons que contem marcas.
# transposons fora de genes
transposons <- read.table(snakemake@input[[1]], sep = "\t")
transposon_in_genes <- read.table(snakemake@input[[2]], sep = "\t")

# ExtraÃ? apenas as colunas de interesse
transposons <- transposons[, c(4, 10)]
transposon_in_genes <- transposon_in_genes[, c(4, 10)]
marks_in_transp <- c(unique(as.character(transposons$V4)),
                     unique(as.character(transposon_in_genes$V4)))

super_exact_background <- length(marks_in_transp)

for(c in 3:(length(snakemake@input) - 1 )){
    
    group_name <- strsplit(as.character(snakemake@input[[c]]), "group_")[[1]][2]
    group_name <- strsplit(as.character(group_name), "_consensus")[[1]][1]
    
    
    # LÃª os arquivos contendo as marcas diferencialmente metiladas do grupo 4
    DM_marks <- read.table(snakemake@input[[c]],
                           header = T,
                           na.strings = "NA",
                           colClasses = "character")
    
    # a busca dos transposons relacionados com as marcas diferencialmente metiladas nas amostras do grupo 4. ConstroÃ? um data.frame com as informaÃ§Ãµes dos transposons de cada amostra.
    transposons_all <- data.frame()
    for(i in 1:length(names(DM_marks))){
        if(i == 1){
            
            # aceita apenas as marcas presentes em transposons
            query <- DM_marks[DM_marks[, i] %in% marks_in_transp, i]
            
            # Execute the analysis
            transposons_group <- transposons[transposons$V4 %in% query,][2]
            colnames(transposons_group) <- paste(names(DM_marks)[i])
            transposons_all <- cbind(transposons_group)
            
        }else{
            
            # aceita apenas as marcas presentes em transposons
            query <- DM_marks[DM_marks[,i] %in% marks_in_transp, i]
            
            # Execute the analysis
            transposons_group <- transposons[transposons$V4 %in% query,][2]
            colnames(transposons_group) <- paste(names(DM_marks)[i])
            transposons_all <- cbindX(transposons_all, transposons_group)
        }
    }
    
    # Escreve o arquivo contendo os transposons de cada marca.
    write.table(transposons_all,
                snakemake@output[[c-2]],
                sep = "\t",
                col.names = T,
                row.names = F,
                quote = F)
    
    
    ## By tissues ##
    
    #### Leaves
    transposons_leaves <- transposons_all[, grep("leaf", colnames(transposons_all))]
    
    
    samples_leaves <- list(A = transposons_leaves[, 1],
                           A = transposons_leaves[, 2],
                           A = transposons_leaves[, 3],
                           A = transposons_leaves[, 4],
                           A = transposons_leaves[, 5])
    
    colnames_leaves <- colnames(transposons_leaves)
    
    for (i in seq_along(samples_leaves)){
        
        names(samples_leaves)[i] <- colnames_leaves[i]
    }
    
    # Removes NAs
    samples_leaves <- lapply(samples_leaves, function(x) x[!is.na(x)])
    
    # Executes the test
    samples_leaves_res <- supertest(samples_leaves, n = super_exact_background)
    
    # Make and save the plot
    svg(paste(paste0("images/transposons_plots/",
                     "group"),
              group_name, 
              "leaf",
              "MSDArT_superexact_plot_TE.svg",
              sep = "_"),
        width = 12,
        height = 8)
    
    plot(samples_leaves_res,
         sort.by = "degree",
         degree = c(1:5),
         legend.col = 1,
         track.area.range = 0.3,
         bar.area.range = 0.10)
    
    dev.off()
    
    # Save the summary file
    plot_name <- paste(paste0("genomic_context_files/",
                              "group"),
                       group_name, 
                       "leaf",
                       "TEs_superexact_table.tst",
                       sep = "_")
    
    write.csv(summary(samples_leaves_res)$Table,
              file = plot_name,
              row.names = FALSE)
    
    #### Xylem
    transposons_xylem <- transposons_all[, grep("wood", colnames(transposons_all))]
    
    
    samples_xylem <- list(A = transposons_xylem[, 1],
                          A = transposons_xylem[, 2],
                          A = transposons_xylem[, 3],
                          A = transposons_xylem[, 4],
                          A = transposons_xylem[, 5])
    
    colnames_xylem <- colnames(transposons_xylem)
    
    for (i in seq_along(samples_xylem)){
        
        names(samples_xylem)[i] <- colnames_xylem[i]
    }
    
    # Removes NAs
    samples_xylem <- lapply(samples_xylem, function(x) x[!is.na(x)])
    
    # Executes the test
    samples_xylem_res <- supertest(samples_xylem, n = super_exact_background)
    
    # Make and save the plot
    svg(paste(paste0("images/transposons_plots/",
                     "group"),
              group_name, 
              "wood",
              "MSDArT_superexact_plot_TE.svg",
              sep = "_"),
        width = 16,
        height = 8)
    
    plot(samples_xylem_res,
         sort.by = "degree",
         degree = c(1:5),
         legend.col = 1,
         track.area.range = 0.3,
         bar.area.range = 0.10)
    
    dev.off()
    
    # Save the summary file
    plot_name <- paste(paste0("genomic_context_files/",
                              "group"),
                       group_name, 
                       "wood",
                       "TEs_superexact_table.tst",
                       sep = "_")
    
    write.csv(summary(samples_xylem_res)$Table,
              file = plot_name,
              row.names = FALSE)
}

## Estimativa para o genoma completo
transp_genome <- read.table(snakemake@input[[5]], sep = "\t")[,4]

#Junta os conjuntos de dados
transposons_all <- cbindX(transposons_all, as.data.frame(transp_genome))

# Verify the TEs that are not possible to classify correctly.
TE_not_class <- read.table(snakemake@input[[5]], sep = "\t")
TE_not_class <- cSplit(indt = TE_not_class, splitCols = "V4", sep = ",", drop = F)

TE_not_class <- TE_not_class[is.na(TE_not_class$V4_02) == FALSE, ]
TE_not_class <- as.character(TE_not_class$V4)

##### 

#GrÃ¡fico representando a classificaÃ§Ã£o dos transposons
## ConstroÃ? a tabela de contagens para cada classe de transposons avaliadas
transp_class_all <- data.frame()
for(i in 1:length(names(transposons_all))){
    
    transp_subset <- unique(as.character(transposons_all[, i]))
    
    # check TEs with overlaps
    transp_unclass <- unique(as.character(transp_subset[transp_subset %in% TE_not_class]))
    
    # check if in the regions with more than one TE, the TEs are of different classification
    transp_unclass_final <- character()
    for( t in 1:length(transp_unclass)){
        
        sum_of_class <- (length(grep("DMX", transp_unclass[t])) + 
                             length(grep("DTX", transp_unclass[t])) + 
                             length(grep("DHX", transp_unclass[t])) + 
                             length(grep("DXX-MtTE", transp_unclass[t])) + 
                             length(grep("DXX_Blc", transp_unclass[t])) + 
                             length(grep("RYX", transp_unclass[t])) + 
                             length(grep("RtX", transp_unclass[t])) + 
                             length(grep("RLX", transp_unclass[t])) + 
                             length(grep("RXX-LARD", transp_unclass[t])) + 
                             length(grep("RXX-TRtM", transp_unclass[t])) + 
                             length(grep("RSX", transp_unclass[t])) + 
                             length(grep("RXX_Blc", transp_unclass[t])))
        
        if(sum_of_class > 1){
            transp_unclass_final <- c(transp_unclass_final, transp_unclass[t])
        }
    }
    
    unknow_TE <- length(transp_unclass_final)
    
    transp_subset <- transp_subset[!transp_subset %in% transp_unclass_final]
    
    #MAVERICK | DMX
    maverick <- length(grep("DMX", transp_subset))
    #CACTA | DTX
    cacta <- length(grep("DTX", transp_subset))
    #HELITRON | DHX
    helitron <- length(grep("DHX", transp_subset))
    #MITE | DXX-MITE
    mite <- length(grep("DXX-MITE", transp_subset))
    #DNA_general | DXX_Blc
    dna_general <- length(grep("DXX_Blc", transp_subset))
    #DIRS/VIPER | RYX
    dirs_viper <- length(grep("RYX", transp_subset))
    #LINE | RIX
    line <- length(grep("RIX", transp_subset))
    #LTR | RLX
    ltr <- length(grep("RLX", transp_subset))
    #LARD | RXX-LARD
    lard <- length(grep("RXX-LARD", transp_subset))
    #TRIM | RXX-TRIM
    trim <- length(grep("RXX-TRIM", transp_subset))
    #SINE | RSX
    sine <- length(grep("RSX", transp_subset))
    #general | RXX_Blc
    rna_general <- length(grep("RXX_Blc", transp_subset))
    
    # Junta os dados de contagem para cada amostra em um Ãºnico data.frame
    category <- c("MAVERICK",
                  "CACTA",
                  "HELITRON",
                  "MITE",
                  "DNA_general",
                  "DIRS/VIPER",
                  "LINE",
                  "LTR",
                  "LARD",
                  "TRIM",
                  "SINE",
                  "RNA_general",
                  "Nested TEs")
    
    len <- c(maverick,
             cacta,
             helitron,
             mite,
             dna_general,
             dirs_viper,
             line,ltr,
             lard,
             trim,
             sine,
             rna_general,
             unknow_TE)
    
    transp_class <- data.frame(rep(names(transposons_all)[i], length(category)), category, len)
    transp_class_all <- rbind(transp_class_all ,transp_class)
}

colnames(transp_class_all) <- c("samples", "classif", "quantif")

transp_class_all <- transp_class_all[!transp_class_all$quantif == "0",]

group_plot <- character()
for(i in 1:nrow(transp_class_all)){
    
    class_of_sample <- strsplit(as.character(transp_class_all$samples[i]), "_")[[1]][[2]]
    
    if(class_of_sample == "genome"){
        
        group_plot <- append(group_plot, "E. grandis")
        
    }else if(class_of_sample == "leaf"){
        
        group_plot <- append(group_plot, "Leaves")
        
    }else if(class_of_sample == "wood"){
        
        group_plot <- append(group_plot, "Wood")
        
    }
}

transp_class_all <- cbind(group_plot, transp_class_all)

write.table(transp_class_all,
            snakemake@output[[11]],
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")

## Defines the colors of the bars.
fill <- c("#000000",
          "#E69F00",
          "#56B4E9",
          "#009E73",
          "#F0E442",
          "#0072B2",
          "#D55E00",
          "#CC79A7",
          "#2C7417",
          "#CBE75F",
          "#A877EA",
          "#6BF4ED",
          "#AC6AC9")

# labels of the grids
levels(transp_class_all$group_plot) <- c("E. grandis" = "E. grandis - all TEs",
                                         "Leaves" = "Methylated TEs - Leaves",
                                         "Wood" = "Methylated TEs - Wood")

transp_plot <- ggplot(transp_class_all, aes(x = samples, y = quantif, fill = classif)) +
    geom_bar(stat = "identity", alpha = 0.9, width = 0.4) +
    scale_x_discrete(name="",
                     breaks = c("transp_genome",
                                "A1_leaf",
                                "B2_leaf",
                                "E5_leaf",
                                "Q8_leaf",
                                "D4_leaf",
                                "A1_wood",
                                "B2_wood",
                                "E5_wood",
                                "Q8_wood",
                                "D4_wood"),
                     labels = c("Genome",
                                "A1",
                                "B2",
                                "E5",
                                "Q8",
                                "D4",
                                "A1",
                                "B2",
                                "E5",
                                "Q8",
                                "D4")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    facet_wrap( ~ group_plot, ncol = 3, scales = "free") +
    ylab("Number of transposons") +
    theme_bw() +
    theme(legend.text = element_text(size = 14),
          axis.title.y = element_text(size = 20, vjust = 2),
          axis.title.x=element_text(size = 14, vjust = 0),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 15)) +
    scale_fill_manual(values = fill)

# Get the ggplot grob
transp_plot_grob = ggplotGrob(transp_plot)

# changes the dimension of the the second and third box
transp_plot_grob$widths[7] = 0.8*transp_plot_grob$widths[7]
transp_plot_grob$widths[9] = 2*transp_plot_grob$widths[9]
transp_plot_grob$widths[13] = 2*transp_plot_grob$widths[13]

# Builds and save the bar plot.
svg(snakemake@output[[12]], width = 18, height = 7)
# Draw the plot
grid.newpage()
grid.draw(transp_plot_grob)
dev.off()

transp_plot_perc <- ggplot(transp_class_all, aes(x = samples, y = quantif, fill = classif)) +
    geom_bar(position = "fill", stat = "identity", alpha = 0.9, width = 0.5) +
    scale_x_discrete(name = "",
                     breaks = c("transp_genome",
                                "A1_leaf",
                                "B2_leaf",
                                "E5_leaf",
                                "Q8_leaf",
                                "D4_leaf",
                                "A1_wood",
                                "B2_wood",
                                "E5_wood",
                                "Q8_wood",
                                "D4_wood"),
                     labels = c("Genome",
                                "A1",
                                "B2",
                                "E5",
                                "Q8",
                                "D4",
                                "A1",
                                "B2",
                                "E5",
                                "Q8",
                                "D4")) +
    scale_y_continuous(labels = percent_format()) +
    facet_wrap( ~ group_plot, ncol = 3, scales = "free") +
    ylab("Relative frequency") +
    theme_bw() +
    theme(legend.text=element_text(size = 14),
          axis.title.y=element_text(size = 20, vjust = 2),
          axis.title.x=element_text(size= 14, vjust = 0),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 15)) +
    scale_fill_manual(values = fill)

# Get the ggplot grob
transp_plot_perc_grob = ggplotGrob(transp_plot_perc)

# changes the dimension of the the second and third box
# transp_plot_perc_grob$widths[7] = 0.8*transp_plot_perc_grob$widths[7]
# transp_plot_perc_grob$widths[11] = 2*t    ransp_plot_perc_grob$widths[11]
# transp_plot_perc_grob$widths[19] = 2*transp_plot_perc_grob$widths[19]

# Builds and save the bar plot.
svg(snakemake@output[[13]], width = 18, height = 7)
# Draw the plot
grid.newpage()
grid.draw(transp_plot_perc_grob)
dev.off()

