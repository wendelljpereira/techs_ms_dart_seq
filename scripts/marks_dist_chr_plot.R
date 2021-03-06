"Usage: marks_dist_chr_plot.R (--out1 <O1>) <input1>
-h --help    show this
--out1  name1   bed file with the fragments generated by double digestion
input1  input1  bed file with the position of the restriction sites

marks_dist_chr_plot.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
# retrieve the command-line arguments
opts <- docopt(doc)

######################
#Carregar bibliotecas#
######################
require(GenomicRanges)
require(ggbio)
require(Gviz)
require(rtracklayer)

# Remove os arquivos que devem ser criados
if (file.exists("file_to_plot")) file.remove("file_to_plot")
if (file.exists("euk.txt")) file.remove("euk.txt")

#################################
#Criar arquivo no formato Grange#
#################################

data <- read.table(opts$`<input1>`, header = F)
colnames(data) <- c("chr", "start", "end", "id", "score", "strand")

bed <- with(data, GRanges(chr, IRanges(start + 1, end), strand, id = id))
bed_merged <- reduce(bed)

##################################################
#Contar transposons que estao dentro de uma range#
##################################################

## O tamanho dos intervalos e n?mero de vezes deve ser sempre setado para englobar
## todo o tamanho do genoma de esp?cie.

vezes <- seq(1, 360)
ini <- -249999
fim <- 0
cromo <- sprintf("%02d", 1:11)
for (i in vezes) {
    ini <- ini + 250000
    fim <- fim + 250000
    result <- restrict(bed_merged, start = ini, end = fim)
    for (j in cromo) {
        atual <- paste("Chr", j, sep = "")
        scor <- length(grep(atual, result))
        line <- paste(atual, ini, fim, "*", scor, sep = "\t")
        write(line, file = "file_to_plot", append = TRUE)
    }
}

#################################
#Utilizar file no Gviz para plot#
#################################
data2 <- read.table("file_to_plot", header = F)
colnames(data2) <- c("chr", "start", "end", "strand", "score")
bed2 <- with(data2, GRanges(chr, IRanges(start, end), strand, score))


###################################
#Setar tamanho do eixo x analisado#
###################################
gtrack <- GenomeAxisTrack()
displayPars(gtrack)$size <- 1.2
displayPars(gtrack)$cex <- 6
displayPars(gtrack)$col.line <- "black"
displayPars(gtrack)$col <- "black"

# Cria o arquivo para dimencionar os cromossomos
pos <- rep("0", 11)
chr <- c("Chr01",
         "Chr02",
         "Chr03",
         "Chr04",
         "Chr05",
         "Chr06",
         "Chr07",
         "Chr08",
         "Chr09",
         "Chr10",
         "Chr11")

len <- c(40297282,
         64237462,
         80088348,
         41978404,
         74731017,
         53893726,
         52447651,
         74330457,
         39019482,
         39359118,
         45510589)

chr_file <- data.frame(chr, pos, len)

write.table(chr_file,
            "euk.txt",
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

########
#Plotar#
########

chr <- import.bed("euk.txt")
ncols <- 3
nrows <- 4
grid.newpage()
png(paste(opts$O1), width = 5000, height = 3000)

## Pushviewport usado para plotar os gr?ficos um ao lado do outro.
pushViewport(viewport(layout = grid.layout(nrows, ncols)))
i <- 0
ii <- 0
k <- 1
for (j in cromo) {
    i <- i + 1
    ii <- ii + 1
    atual <- paste("Chr", j, sep = "")
    nome <- paste("Bands_distribution_in_", atual, sep = "")
    chr_idx <- levels(seqnames(chr)) == atual
    f <- 1
    t <- end(chr[chr_idx])
    
    dTrack <- DataTrack(bed2,
                        name = atual,
                        chromosome = atual,
                        ylim = c(0, 100),
                        fill.histogram = c("#0015ff", "#556B2F"),
                        cex.legend = 10)
    
    displayPars(dTrack)$cex.axis <- 4
    displayPars(dTrack)$fontsize.legend <- 4
    displayPars(dTrack)$cex <- 4
    displayPars(dTrack)$fontcolor.legend <- "black"
    pushViewport(viewport(layout.pos.col = ii, layout.pos.row = k))
    
    plotTracks(list(dTrack, gtrack),
               littleTicks = FALSE,
               add = TRUE,
               from = f,
               to = t,
               chromosome = atual,
               type = "histogram",
               col.axis = "black",
               background.panel = "transparent",
               showTitle = TRUE,
               fontcolor = "black",
               background.title = "transparent",
               main = atual,
               cex.main = 8,
               fontcolor.legend = "black",
               legend = FALSE)
    popViewport(1)
    if (i %% 3 == 0){
        k <- k + 1
        ii <- 0
    }
}
dev.off()