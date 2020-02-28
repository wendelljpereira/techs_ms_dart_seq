save.image("tec_rep_snakemake.rda")

require(tidyverse)

# calculate the average of the counts of samples with techinical replicates.
counts_df <- read.table(snakemake@input[[1]],
                        sep = "\t",
                        header = T,
                        check.names = F)

# Names of samples with techinical replicates
samples <- snakemake@params[["samples_with_tec_reps"]]

cor_table  <- data.frame()
for (i in samples){
    
    # subset containing the counts of the sample
    counts_sub <- counts_df[, grep(i, colnames(counts_df))]
    
    ## Estimates the correlation between technical replicates
    cor_counts_sub <- cor.test(counts_sub[, 1], counts_sub[,  2], method = "pearson")
    
    cor_table_sub <- data.frame(sample_name = i,
                                Pearson_corr = cor_counts_sub[[4]])
    
    cor_table <- rbind(cor_table, cor_table_sub)
    
    # Estimate the average of the replicates
    counts_sub$aver <- round((as.numeric(counts_sub[,1]) + 
                                  as.numeric(counts_sub[,2])) / 2,)
    
    # Exclude the original replicates of the dataset
    counts_df <- counts_df[, -grep(i, colnames(counts_df))]
    
    # Add the average to the dataset
    counts_df <- cbind(counts_df, counts_sub[, ncol(counts_sub)])
    
    colnames(counts_df)[ncol(counts_df)] <- paste(i)
}

write.table(cor_table,
            snakemake@output[[1]],
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)

# Change the names for the standard format (4 pieces)
col_names <- colnames(counts_df)

colnames_final <- character()
for(c in 2:length(col_names)){
    
    colnames_sub <- str_split(col_names[c], "_")
    
    colnames_sub_corrected <- paste(colnames_sub[[1]][1],
                                    colnames_sub[[1]][3],
                                    colnames_sub[[1]][4],
                                    colnames_sub[[1]][5],
                                    sep = "_")
    
    colnames_final <- c(colnames_final,
                        colnames_sub_corrected)
}
colnames(counts_df)[2:length(counts_df)] <- colnames_final

write.table(counts_df,
            snakemake@output[[2]],
            col.names = T,
            row.names = F,
            sep = "\t",
            quote = F)
