
# Script para análise no DEseq2 dos dados de contagem DArT.
## O script possuí uma função generica que deve atender a demanda de análise para todas as amostras dos 5 grupos propostos no delineamento experimental.

# Carrega os pacotes de interesse e define o diretório de análise
require(ggplot2)
require(DESeq2)
require(plyr)
require(gdata)
require(tidyr)

save.image(file = "deseq2_snakemake.rda")
#load(file = "deseq2_snakemake.rda")

# Carrega a função para expressão diferencial
deseq2_with_counts <- function(file = file,
                               prefix = prefix,
                               clone_name = clone_name,
                               tissue = tissue,
                               fdr = 0.001,
                               log_fc = 1,
                               sep_into = sep_into,
                               subset_model = "normal",
                               no_bio_rep = "FALSE",
                               dispersion = NULL,
                               filter = "posit",
                               min_msp = min_msp,
                               intersect_file = intersect_file){
    
    # Leitura do conjunto de dados.
    data <- read.table(as.character(file),
                       header = T,
                       sep = ",",
                       check.names = F)
    
    #Seleciona o nome das marcas
    marcas <- as.character(data[, 1])
    prefixo <- as.character(prefix)
    clone_name <- as.character(clone_name)
    tecido <- as.character(tissue)
    fdr <- as.numeric(fdr)
    log_fc <- as.numeric(log_fc)
    sep_into <- as.numeric(sep_into)
    subset_model <- as.character(subset_model)
    filter <- as.character(filter)
    min_msp <- as.numeric(min_msp)
    
    # Seleciona o clone, tecido e grupo amostral que será utilizado
    clone <- data[, grep(pattern = clone_name, colnames(data))]
    clone <- clone[, grep(pattern = tecido, colnames(clone))]
    clone <- clone[, grep(pattern = paste(g, clone_name, sep = "_"), colnames(clone))]
    
    #Renomeia as linhas com os nomes das marcas
    clone <- as.data.frame(clone, row.names = marcas)
    
    print(paste("These are the selected samples:", colnames(clone)))
    
    # Seleciona o modelo de análise a ser realizado
    if (subset_model == "intersect"){
        
        intersect <- read.table(intersect_file,
                                header = T,
                                quote = "\"",
                                sep = "\t",
                                check.names = F,
                                colClasses = "character")
        
        intersect <- as.character(intersect[, 1])
        
        clone <- clone[rownames(clone) %in% intersect, ]
        print("subset_model = intersect")
        
    } else {
        
        print("subset_model = normal")
        
    }
    
    # Escreve um arquivo contendo a lista de sítios testados para metilação
    tested_sites <- rownames(clone)
    
    write.table(tested_sites,
                "tested_sites_to_methylation.txt",
                col.names = F,
                row.names = F,
                quote = F,
                sep = "\t")
    
    # Aplica o filtro de contagem minima aceitavel para MspI, aceitando apenas marcas cujo a média de contagens nas amostras com mspI é maior do que o determinado pelo usuario
    if (filter == "posit"){
        if (no_bio_rep == "FALSE"){
            
            filtred_clone <- data.frame()
            
            for (i in 1:nrow(clone)){
                
                count_average <- clone[i, ]
                count_average <- count_average[, grep(pattern = "ms", colnames(clone))]
                
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                    
                }
            }
        }else if (no_bio_rep == "TRUE"){
            
            filtred_clone <- data.frame()
            
            for (i in 1:nrow(clone)){
                
                count_average <- clone[i, ]
                count_average <- as.data.frame(count_average[, grep(pattern = "ms", colnames(clone))])
                
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                    
                }
            }
        }
    }else if (filter == "negat"){
        if (no_bio_rep == "FALSE"){
            
            filtred_clone <- data.frame()
            
            for (i in 1:nrow(clone)){
                
                count_average <- clone[i, ]
                count_average <- count_average[, grep(pattern = "hp", colnames(clone))]
                
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                }
            }
        }else if (no_bio_rep == "TRUE"){
            filtred_clone <- data.frame()
            for (i in 1:nrow(clone)){
                count_average <- clone[i, ]
                count_average <- as.data.frame(count_average[, grep(pattern = "hp", colnames(clone))])
                if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp){
                    
                    filtred_clone <- rbind(filtred_clone, clone[i, ])
                    
                }
            }
        }
    }
    
    print(paste("For sample",
                paste("'", clone_name, "_", tecido, "'", sep = ""),
                nrow(filtred_clone),
                "marks passed for all filters!"))
    
    msp_10_sites <- as.data.frame(rownames(filtred_clone))
    
    write.table(msp_10_sites,
                paste(clone_name, "_", tecido, "_", "msp_bigger_than_threshold.txt", sep = ""),
                row.names = F,
                col.names = F,
                quote = F)
    
    # # Transfere os clones filtrados para a variavel clone
    # clone <- lapply(filtred_clone, as.integer)
    # clone <- as.data.frame(clone, check.names = FALSE)
    # 
    # #Renomeia as linhas com os nomes das marcas
    # clone <- as.data.frame(clone, row.names = rownames(filtred_clone))
    clone <- filtred_clone
    head(clone)
    
    
    # Define a função que implementará as estatísticas descritivas
    descritivas <- function(y) {
        data.frame(n = length(y),
                   mean = mean(y),
                   IC95inf = t.test(y)$conf.int[1],
                   IC95sup = t.test(y)$conf.int[2],
                   var = var(y),
                   std.dev. = sd(y),
                   std.err. = sd(y) / sqrt(length(y)),
                   CV = 100 * sd(y) / mean(y),
                   min = min(y),
                   max = max(y),
                   Q1 = quantile(y, 0.25),
                   median = median(y),
                   Q3 = quantile(y, 0.75),
                   sum(y))
    }
    
    # Determina os grupos aos quais pertencem cada amostra, baseando-se no último fator após o separador "_".
    ##Lê os nomes das amostras, salva em um data frame e nomeia a coluna para o próximo comando.
    names <- as.data.frame(names(clone))
    colnames(names) <- "nome"
    
    # Quebra o texto em cada separador "_" e salva em colunas distintas no dataframe
    sep <- paste0("x", 1:(sep_into + 1))
    
    names_group <- separate(names,
                            nome,
                            into = c(paste0("null", 1:(sep_into)), "Enzime"),
                            sep = "_",
                            remove = T,
                            extra = "merge")
    
    # Salva a coluna com os nomes das enzimas como um vetor de caracteres
    names_group <- as.character(names_group$Enzime)
    
    # Remove o ".1" das réplicas e salva os nomes das enzimas como um vetor que poderá ser passado ao deseq2
    x <- as.character()
    for (i in 1:length(names_group)){
        x <- rbind(x, substr(names_group[i], 0, 2))
    }
    groups <- as.character(x)
    groups <- gsub(".", "", groups, fixed = T)
    
    ##########
    # construção do modelo para o DEseq2
    coldata <- data.frame(condition = groups,
                          type = rep("single-read", length(groups)),
                          row.names = paste(colnames(clone)))
    
    # Carrega os dados para o objeto adequado ao deseq2
    data_deseq <- DESeqDataSetFromMatrix(countData = clone,
                                         colData = coldata,
                                         design = ~ condition)
    
    # Determina qual o conjunto será utilizado como referência
    data_deseq$condition <- relevel(data_deseq$condition, ref = "hp")
    
    # Executa a análise
    data_deseq <- DESeq(data_deseq)
    
    # Extraí os resultados
    raw_results <- as.data.frame(results(data_deseq))
    
    # Filtra pelo FDR e log do fold change passado.
    ## obs. a filtragem será determinada pelo usuário, podendo filtrar apenas os valores de foldchange positivos (maiores com a enzima MspI), negativos (maiores na enzima HpaII) ou o conjunto completo, contendo todas as marcas.
    if (filter == "posit"){
        
        results_sig <- subset(raw_results, raw_results$padj <= fdr)
        results_sig <- subset(results_sig, results_sig$log2FoldChange >= log_fc)
        
    }else if (filter == "negat"){
        
        results_sig <- subset(raw_results, raw_results$padj <= fdr)
        results_sig <- subset(results_sig, results_sig$log2FoldChange <= 0)
        results_sig <- subset(results_sig, abs(results_sig$log2FoldChange) >= log_fc)
        
    }else if (filter == "full"){
        
        results_sig <- subset(raw_results, raw_results$padj <= fdr)
        results_sig <- subset(results_sig, abs(results_sig$log2FoldChange) >= log_fc)
        
    }else{
        
        print("parâmetro filter incorreto!")
        
    }
    
    # Testa se foram encontradas marcas diferencialmente expressas para os parâmetros implementados.
    if (nrow(results_sig) == 0){
        
        print(paste("Não existem marcas diferencialmente expressas para o",
                    prefix,
                    clone_name,
                    "tecido:",
                    tissue,
                    "com os parâmetros implementados!",
                    sep = " "))
        
    }else{
        
        ## Armazena as marcas diferencialmente expressas em um vetor.
        marcks_sig <- rownames(results_sig)
        Data_DE <- as.data.frame(marcks_sig)
        colnames(Data_DE) <- paste(clone_name, tecido, sep = "_")
        
        #testa se o arquivo de contagens existe e adiciona os novos dados.
        if (file.exists(paste(prefixo, "DE_marks.txt", sep = "_")) == "FALSE"){
            
            print(paste("O arquivo",
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        "não existe! Um novo arquivo será gerado."))
            
            write.table(Data_DE,
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = F,
                        col.names = T)
            
        }else if (file.exists(paste(prefixo, "DE_marks.txt", sep = "_")) == "TRUE"){
            
            print(paste("O arquivo",
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        "foi detectado! A saída será adicionada."))
            
            ## Lê o conjunto de marcas diferencialmente expressas existente e grava o novo conjunto.)
            ## Lê o arquivo do grupo.
            Data_DE_g <- read.table(paste(prefixo, "DE_marks.txt", sep = "_"),
                                    header = T,
                                    sep = "\t")
            
            ## Adiciona o novo conjunto de dados e reescreve o arquivo no diretório
            Data_DE_g <- cbindX(Data_DE_g, Data_DE)
            
            write.table(Data_DE_g,
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = F,
                        col.names = T)
        }
        
        # Estatítiscas descritivas das contagens obtidas das marcas diferencialmente expressas ###
        # Recuperando contagens dos clones que foram considerados diferencialmente expressos
        factor <- rownames(results_sig)
        counts <- data.frame()
        for (i in 1:length(factor)){
            counts <- rbind(counts, subset(clone, row.names(clone) == factor[i]))
        }
        
        #Aplica a função descritivas e armazena os resultados em um data frame
        data_stats <- apply(counts, 2, descritivas)
        data_stats <- do.call("rbind", data_stats)
        
        # Testa se o arquivo com as estatisticas do grupo já existe e adiciona os resultados.
        if (file.exists(paste(prefixo, "DE_stats.txt", sep = "_")) == "FALSE"){
            
            print(paste("O arquivo",
                        file = paste(prefixo, "DE_stats.txt", sep = "_"),
                        "não existe! Um novo arquivo será gerado."))
            
            write.table(data_stats,
                        file = paste(prefixo, "DE_stats.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)
            
        }else if (file.exists(paste(prefixo, "DE_stats.txt", sep = "_")) == "TRUE"){
            
            print(paste("O arquivo",
                        file = paste(prefixo, "DE_stats.txt", sep = "_"),
                        "foi detectado! A saída será adicionada."))
            
            # Carrega o arquivo do grupo, adiciona as estatistícas e reescreve o arquivo
            data_stats_g <- read.table(paste(prefixo, "DE_stats.txt", sep = "_"),
                                       header = T,
                                       sep = "\t")
            
            data_stats_g <- rbind(data_stats_g, data_stats)
            
            write.table(data_stats_g,
                        file = paste(prefixo, "DE_stats.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = T,
                        col.names = T)
        }
    }
}

samples_without_rep <- snakemake@params[["samples_without_rep"]]

# executa a função realizando a expressão diferencial de todas as amostras no grupo 4
for (g in as.numeric(snakemake@params[["groups"]])){
    for (c in 1:length(snakemake@params[["genotypes"]])){
        for (t in 1:length(snakemake@params[["tissues"]])){
            deseq2_with_counts(file = snakemake@input[1],
                               prefix = paste(snakemake@params[["prefix"]], "_", g, "_",
                                                       snakemake@params[["subset_model"]], sep = ""),
                               clone_name = snakemake@params[["genotypes"]][c],
                               tissue = snakemake@params[["tissues"]][t],
                               fdr = snakemake@params[["fdr"]],
                               log_fc = snakemake@params[["log_fold_change"]],
                               sep_into = snakemake@params[["sep_into"]],
                               subset_model = snakemake@params[["subset_model"]],
                               no_bio_rep = snakemake@params[["no_bio_rep"]],
                               dispersion = snakemake@params[["dispersion"]],
                               filter = snakemake@params[["filtration_mode"]],
                               min_msp = as.numeric(snakemake@params[["min_msp"]]),
                               intersect_file = paste(snakemake@input[2])
            )
        }
    }
}
