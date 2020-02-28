"Usage: id_part4.R (--out1 <O1>) (--out2 <O2>) (--out3 <O3>) (--out4 <O4>) <input1> <input2> <input3> <input4>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
--out3   name3    specify the name for the third output file
--out4   name4    specify the name for the third output file
id_part4.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)

####### Parte 4 ########
marcas_exons <- read.table(opts$`<input1>`, colClasses = "character")

marcas_exons_com_transposons <- read.table(opts$`<input2>`,
                                           colClasses = "character")

marcas_exons_com_transposons_names <- unique(marcas_exons_com_transposons$V4)
marcas_exons_com_transposons <- marcas_exons[marcas_exons$V4 %in% marcas_exons_com_transposons_names, ]

marcas_exons_sem_transposons <- marcas_exons[!marcas_exons$V4 %in% marcas_exons_com_transposons_names, ]

write.table(marcas_exons_sem_transposons,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(marcas_exons_com_transposons,
            paste(opts$O2),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

marcas_fora_de_exons <- read.table(opts$`<input3>`, colClasses = "character")

marcas_fora_de_exons_com_transposons <- read.table(opts$`<input4>`, colClasses = "character")
marcas_fora_de_exons_com_transposons_names <- unique(marcas_fora_de_exons_com_transposons$V4)

marcas_nao_condantes_sem_transposons <- marcas_fora_de_exons[!marcas_fora_de_exons$V4 %in% marcas_fora_de_exons_com_transposons_names, ]

marcas_nao_condantes_com_transposons <- marcas_fora_de_exons[marcas_fora_de_exons$V4 %in% marcas_fora_de_exons_com_transposons_names, ]

write.table(marcas_nao_condantes_com_transposons,
            paste(opts$O3),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(marcas_nao_condantes_sem_transposons,
            paste(opts$O4),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")