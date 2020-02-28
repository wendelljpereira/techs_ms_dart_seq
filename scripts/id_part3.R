"Usage: id_part3.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) <input1> <input2>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
--out3   name3    specify the name for the third output file
id_part3.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)

######## Parte 3 #########
require(splitstackshape)
marcas_in_genes <- read.table(opts$`<input1>`, colClasses = "character")
marcas_in_exons_feature <- read.table(opts$`<input2>`, colClasses = "character")

marcas_in_exons_feature_names <- unique(marcas_in_exons_feature$V4)

methylated_features_all <- data.frame()
methylated_features_problems <- data.frame()

for (i in marcas_in_exons_feature_names){
    
    target <- i
    target_data <- marcas_in_exons_feature[marcas_in_exons_feature$V4 == target, ][15]
    target_exon <- cSplit(indt = target_data,
                          splitCols = c("V15"),
                          sep = ";",
                          type.convert = F)
    
    target_exon <- cSplit(indt = target_exon,
                          splitCols = c("V15_1"),
                          sep = ".v2.0.",
                          type.convert = F)
    
    target_exon <- cSplit(indt = target_exon,
                          splitCols = c("V15_1_1"),
                          sep = "=",
                          type.convert = F)
    
    target_gene <- unique(substr(target_exon$V15_1_1_2, 0, 12))
    target_exon <- unique(substr(target_exon$V15_1_2, 0, 12))
    
    if (length(target_gene) == 1){
        
        methylated_features <- cbind(marcas_in_genes[marcas_in_genes$V4 == i, ],
                                     as.data.frame(target_gene),
                                     paste(t(as.data.frame(target_exon)), sep = "", collapse = ";"))
        
        colnames(methylated_features) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
        methylated_features_all <- rbind(methylated_features_all, methylated_features)
        
    }else{
        
        methylated_features <- cbind(marcas_in_genes[marcas_in_genes$V4 == i, ],
                                     as.data.frame(target_gene),
                                     paste(t(as.data.frame(target_exon)), sep = "", collapse = ";"))
        
        colnames(methylated_features) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
        methylated_features_problems <- rbind(methylated_features_problems,
                                              methylated_features)
    }
}

# Escreve o arquivo bed com as marcas em exons de genes únicos
write.table(methylated_features_all,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Escreve o arquivo bed com as marcas em exons de genes que sobrepõem outros genes
write.table(methylated_features_problems,
            paste(opts$O2),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Marcas que estão em genes, mas fora de exons
marcas_fora_de_exons <- marcas_in_genes[!marcas_in_genes$V4 %in% marcas_in_exons_feature_names, ]

# Escreve o bed das marcas em genes, mas fora de exons
write.table(marcas_fora_de_exons,
            paste(opts$O3),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")