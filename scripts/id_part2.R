"Usage: id_part2.R (--out1 <O1>) (--out2 <O2>) <input1> <input2>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
id_part2.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)
save.image()

######## Parte 2 #########
# Lê os arquivos bed com as features
marcas <- read.table(opts$`<input1>`, colClasses = "character")

# Lê as marcas que não estão em genes e possuem intersecção com transposons
marcas_transposon_int <- read.table(opts$`<input2>`,
                                    sep = "\t",
                                    colClasses = "character")

marcas_transposons_int_names <- unique(marcas_transposon_int$V4)

#Cria um arquivo com as marcas que estão em transposons.
marcas_nos_transposons <- marcas[marcas$V4 %in% marcas_transposons_int_names, ]

write.table(marcas_nos_transposons,
            paste(opts$O1),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)

# Remove as marcas que estão em transposons do conjunto restante. 
marcas_fora_de_genes_e_transposons <- marcas[!marcas$V4 %in% marcas_transposons_int_names, ]

write.table(marcas_fora_de_genes_e_transposons,
            paste(opts$O2),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)