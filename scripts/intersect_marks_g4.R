# load the docopt library
require(gdata)

save.image(file = "intersect_snakemake.rda")

# Carrega a função que verifica quais marcas possuem contagem em MspI para ambas as amostras que serão comparadas.
intersect_marks <- function(data_intersect = data_intersect,
                            sample_name = "all", # List the clones that should be compared
                            tissue = "all",
                            group = "all",
                            prefix = "file",
                            enzyme = "ms",
                            int_name = "all",
                            non_overlap=non_overlap,
                            non_overlap_is_file=non_overlap_is_file
){
    
    # Loads the parameters
    sample_name <- as.character(c(sample_name))
    tissue <- as.character(c(tissue))
    group <- as.character(c(group))
    prefix <- as.character(prefix)
    enzyme <- as.character(c(enzyme))
    int_name <- as.character(int_name)
    
    # Realiza a leitura dos dados
    clone <- data_intersect
    
    # Subset the samples to keep only the ones listed
    if (sample_name == "all"){
        
        print("all clones will be use!")
        
    } else {
        for (i in 1:length(sample_name)){
            
            y <- clone[, grep(pattern = sample_name[i], colnames(clone))]
            
            if (i == 1){
                
                x <- as.data.frame(y)
                
            } else{
                
                x <- cbindX(x, y)
                
            }
        }
        # renomeando o conjunto de dados
        clone <- x
    }
    
    # Select the samples based on the tissue.
    if (tissue == "all"){
        
        print("all tissues will be use!")
        
    }else{
        for (i in 1:length(tissue)){
            
            y <- clone[, grep(pattern = tissue[i], colnames(clone))]
            
            if (i == 1){
                
                x <- as.data.frame(y)
                
            }else{
                
                x <- cbindX(x, y)
                
            }
        }
        # renomeando o conjunto de dados
        clone <- x
    }
    
    # Select the samples based on the group
    if (group == "all"){
        
        print("all groups will be use!")
        
    }else{
        
        for (g in 1:length(group)){
            
            y <- clone[, substr(names(clone), 1, 2) == paste(group[g], "_", sep = "")]
            
            if (g == 1){
                
                x <- as.data.frame(y)
                
            }else{
                
                x <- cbindX(x, y)
                
            }
        }
        # renomeando o conjunto de dados
        clone <- x
    }
    
    # Select the samples based on the enzyme used to generate the libraries
    if (enzyme == "names"){
        
        print("all enzymes data will be use!")
        
    }else{
        
        for (i in 1:length(enzyme)){
            
            y <- clone[, grep(pattern = enzyme[i], colnames(clone))]
            
            if (i == 1){
                
                x <- as.data.frame(y)
                
            }else{
                
                x <- cbindX(x, y)
                
            }
        }
        # renomeando o conjunto de dados
        clone <- x
    }
    
    # Removes the missing data
    x <- as.list(clone)
    x <- lapply(x, function(x) x[!is.na(x)])
    
    # Creates the intersection among the selected samples
    Int <- Reduce(intersect, x)
    I <- as.data.frame(Int)
    
    # name the column with the data
    colnames(I) <- paste(int_name)
    
    #testa se o arquivo com as intersecções existe e adiciona os novos dados. 
    if (file.exists(paste(prefix, "intersect_marks.txt", sep = "_")) == "FALSE"){
        
        print(paste("O arquivo",
                    file = paste(prefix, "intersect_marks.txt", sep = "_"),
                    "não existe! Um novo arquivo será gerado."))
        
        write.table(I,
                    file = paste(prefix, "intersect_marks.txt", sep = "_"),
                    sep = "\t",
                    quote = F,
                    row.names = F,
                    col.names = T)
        
    }else if (file.exists(paste(prefix, "intersect_marks.txt", sep = "_")) == "TRUE"){
        print(paste("O arquivo", file = paste(prefix, "intersect_marks.txt", sep = "_"), "foi detectado! A saída será adicionada."))
        ## Lê o conjunto de marcas da interseccao e grava o novo conjunto.)
        ## Lê o arquivo do grupo.
        Data_int_g <- read.table(paste(prefix, "intersect_marks.txt", sep = "_"), header = T, sep = "\t")
        ## Adiciona o novo conjunto de dados e reescreve o arquivo no diretório
        Data_int_g <- cbindX(Data_int_g, as.data.frame(I))
        write.table(Data_int_g, file = paste(prefix, "intersect_marks.txt", sep = "_"), sep = "\t", quote = F, row.names = F, col.names = T)
    }
}

# Lê todos os sitios dos multiplos arquivos e salva em um unico data.frame
all_sites <- data.frame()
for (i in 2:length(snakemake@input)){
    if (i == 2){
        
        sites <- read.table(snakemake@input[[i]],
                            header = T,
                            sep = "\t",
                            quote = "\"",
                            check.names = F,
                            na.strings = "NA",
                            colClasses = "character")
        all_sites <- sites
        
    }else{
        
        sites <- read.table(snakemake@input[[i]],
                            header = T,
                            sep = "\t",
                            quote = "\"",
                            check.names = F,
                            na.strings = "NA",
                            colClasses = "character")
        
        all_sites <- cbindX(all_sites, sites)
    }
}

# Calculates the intersection between the groups, for each sample and tissue.
for (n in 1:length(snakemake@params[["names"]])){
        for (t in 1:length(snakemake@params[["tissue"]])){
            intersect_marks(data_intersect = all_sites,
                            sample_name = snakemake@params[["names"]][n],
                            tissue = snakemake@params[["tissue"]][t],
                            prefix = snakemake@params[["prefix_samples"]],
                            enzyme = snakemake@params[["enzyme"]][1],
                            int_name = paste(snakemake@params[["names"]][n],
                                             snakemake@params[["tissue"]][t],
                                             sep = "_"))
        }
}

# Calculates the intersection between the clones, for each group and tissue.
for (g in 1:length(snakemake@params[["groups"]])){
    for (t in 1:length(snakemake@params[["tissue"]])){
        intersect_marks(data_intersect = all_sites,
                        tissue = snakemake@params[["tissue"]][t],
                        prefix = snakemake@params[["prefix_groups"]],
                        enzyme = snakemake@params[["enzyme"]][1],
                        group = snakemake@params[["groups"]][g],
                        int_name = paste("group",
                                         snakemake@params[["groups"]][g],
                                         snakemake@params[["tissue"]][t],
                                         sep = "_"))
    }
}

# Calculates the intersection between all samples of each tissue
for (t in 1:length(snakemake@params[["tissue"]])){
    
    intersect_marks(data_intersect = all_sites,
                    tissue = snakemake@params[["tissue"]][t],
                    prefix = paste("sites_on_the_intersection/all_samples", snakemake@params[["tissue"]][t], sep = "_"),
                    enzyme = snakemake@params[["enzyme"]][1],
                    int_name = paste("all",
                                     snakemake@params[["tissue"]][t],
                                     sep = "_"))
    
}
