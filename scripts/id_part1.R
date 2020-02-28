"Usage: id_part1.R (--out1 <O1>) (--out2 <O2>) <input1> <input2>
-h --help    show this
--out1   name1    specify the name for the first output file
--out2   name2    specify the name for the second output file
id_part1.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)
save.image()

# Verificação do contexto de cada marca no genoma de E. grandis
# Lê os arquivos bed com as features
marcas <- read.table(opts$`<input1>`, colClasses = "character")

# Lê as marcas com intersecção com genes
marcas_genes_int <- read.table(opts$`<input2>`,
                               sep = "\t",
                               colClasses = "character")

marcas_genes_int_names <- unique(marcas_genes_int$V4)

# Seleciona as marcas que estão em genes e escreve um novo bed
marcas_em_genes <- marcas[marcas$V4 %in% marcas_genes_int_names, ]

write.table(marcas_em_genes,
            paste(opts$O1),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)

# Remove as marcas que estão em genes do conjunto total e escreve um novo bed
marcas_fora_dos_genes <- marcas[!marcas$V4 %in% marcas_genes_int_names, ]

write.table(marcas_fora_dos_genes,
            paste(opts$O2),
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F)