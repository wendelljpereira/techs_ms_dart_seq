require(gdata)

save.image("snakemake.rda")

dados <- read.table(snakemake@input[[1]],
                    header = T,
                    sep = ",",
                    check.names = F)

# Select all sites with counts bigger than one in each column.

for (i in 2:ncol(dados)){
    if (i == 2){
        data_without_zero_all_samples <- dados[dados[, i] > 0, ]
        all_marks <- as.data.frame(data_without_zero_all_samples[, 1])
        colnames(all_marks) <- paste(names(dados)[i])
    }else{
        data_without_zero_sub <- dados[dados[, i] > 0, ]
        all_marks_sub <- as.data.frame(data_without_zero_sub[, 1])
        colnames(all_marks_sub) <- paste(names(dados)[i])
        all_marks <- cbindX(all_marks, all_marks_sub)
    }
}

# Separates the columns (samples) accordingly with the experimental groups (first name before the separator "_").
## Writes a file for each group.
y <- data.frame()
for (c in 1:length(names(all_marks))){
    x <- as.data.frame(paste(strsplit(names(all_marks)[c], "_")[[1]][1]))
    colnames(x) <- "group"
    y <- rbind(y, x)
}
y$group <- as.factor(y$group)

for (l in 1:length(levels(y$group))){
    group <- all_marks[y$group == sort(levels(y$group))[l]]
    
    write.table(group,
                snakemake@output[[l]],
                row.names = F,
                quote = F,
                sep = "\t")
}
