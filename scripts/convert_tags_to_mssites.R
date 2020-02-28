# comparação entre Deseq2 e edgeR
require(gdata)
require(VennDiagram)
require(gridExtra)
save.image("snakemake.Rda")

# Define a funcao que removera as redundancias entre os sitios suportados por mais de um fragmento
correct_sites <- function(sites_bed_file = sites_bed_file, data = data){
    
    all_sites <- read.table(sites_bed_file, sep = "\t")
    all_sites <- unique(as.character(all_sites[, 4]))
    
    # para lidar com sitios suportados por dois fragmentos é necessário comparar com o arquivo das posicoes. Substituir o nome do sitio e depois remover a redundancia.
    dup_sites <- all_sites[grep(pattern = "|", all_sites, fixed = T)]
    
    dup_sites_names <- strsplit(x = dup_sites, split = "|", fixed = T)
    dup_sites_names <- unique(c(do.call("rbind", dup_sites_names)))
    
    new_data_corrected <- data.frame()
    for (i in 1:ncol(data)){
        new_data_elements <- data[i]
        new_data_renamed_all <- data.frame()
        for (j in 1:nrow(new_data_elements)){
            
            if (as.character(new_data_elements[j, 1]) %in% dup_sites_names){
                
                x <- grep(pattern = paste0("\\<", new_data_elements[j, 1], "\\|"), x = dup_sites)
                y <- grep(pattern = paste0("\\|", new_data_elements[j, 1], "\\>"), x = dup_sites)
                
                if (length(x) > 0){
                    z <- x
                }else if (length(y) > 0){
                    z <- y
                }else{
                    print("Erro!")
                }
                
                new_data_renamed <- as.data.frame(dup_sites[z])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all, new_data_renamed) # Apenas o x ou o y terão valor, nunca os dois.
            }else{
                new_data_renamed <- as.data.frame(new_data_elements[j, 1])
                colnames(new_data_renamed) <- "new_data_renamed"
                new_data_renamed_all <- rbind(new_data_renamed_all,
                                              new_data_renamed)
            }
            
        } 
        if (i == 1){
            colnames(new_data_renamed_all) <- colnames(data)[i]
            new_data_corrected <- as.data.frame(unique(new_data_renamed_all))
            
        }else{
            colnames(new_data_renamed_all) <- colnames(data)[i]
            new_data_corrected <- cbindX(new_data_corrected,
                                         as.data.frame(unique(new_data_renamed_all)))
        }
    }
    
    return(new_data_corrected)
}

deseq_dm_marks_g2 <- read.table(snakemake@input[[1]],
                                header = T,
                                colClasses = "character")

deseq_dm_marks_g3 <- read.table(snakemake@input[[2]],
                                header = T,
                                colClasses = "character")

# testa se os sítios metilados são um dos que possuem suporte de mais de um fragmento e corrige a redundancia
sites_bed_file <- snakemake@input[[3]]
deseq_dm_marks_g2_corrected <- correct_sites(sites_bed_file = sites_bed_file,
                                             data = deseq_dm_marks_g2)

deseq_dm_marks_g3_corrected <- correct_sites(sites_bed_file = sites_bed_file,
                                             data = deseq_dm_marks_g3)

write.table(deseq_dm_marks_g2_corrected,
            snakemake@output[[1]],
            row.names = F,
            col.names = T,
            sep = "\t",
            quote = F)

write.table(deseq_dm_marks_g3_corrected,
            snakemake@output[[2]],
            row.names = F,
            col.names = T,
            sep = "\t",
            quote = F)
