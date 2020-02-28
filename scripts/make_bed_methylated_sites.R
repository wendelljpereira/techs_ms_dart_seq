save.image("snakemake.rda")

# Reads the set of all MSD-Sites
sites <- read.table(snakemake@input[[1]], sep = "\t")

# Reads the MSD-Methylated sites
marks_unique_final <- data.frame()
for(i in 2:length(snakemake@input)){
    
    marks_sub <- read.table(snakemake@input[[i]], header = T, sep = "\t")
    
    marks_unique <- as.data.frame(unique(as.character(as.matrix(marks_sub))))
    marks_unique <- na.omit(marks_unique)
    
    marks_pos <- sites[sites$V4 %in% marks_unique[, 1], ]
    
    marks_unique_final <- rbind(marks_unique_final, marks_unique)
}

marks_unique <- as.data.frame(unique(as.character(as.matrix(marks_unique_final))))
marks_unique <- na.omit(marks_unique)

# ExtrÃ¡i as posiÃ§Ãµes das marcas metiladas
marks_pos <- sites[sites$V4 %in% marks_unique[, 1], ]

write.table(marks_pos,
            snakemake@output[[1]],
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

# Select the positions in plus strand as representative of the methylation position
marks_pos_plus <- marks_pos[marks_pos$V6 == "+", ]

write.table(marks_pos_plus,
            snakemake@output[[2]],
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")