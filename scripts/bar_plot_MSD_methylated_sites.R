save.image("snakemake_bar_plot.rda")

require(ggplot2)
require(ggthemes)

## A1
A1_MG_and_PR_Leaf <- read.table(snakemake@input[[1]])
A1_MG_Leaf <- read.table(snakemake@input[[6]])
A1_PR_Leaf <- read.table(snakemake@input[[11]])
A1_MG_and_PR_Xylem <- read.table(snakemake@input[[16]])
A1_MG_Xylem <- read.table(snakemake@input[[21]])
A1_PR_Xylem <- read.table(snakemake@input[[26]])

## B2
B2_MG_and_PR_Leaf <- read.table(snakemake@input[[2]])
B2_MG_Leaf <- read.table(snakemake@input[[7]])
B2_PR_Leaf <- read.table(snakemake@input[[12]])
B2_MG_and_PR_Xylem <- read.table(snakemake@input[[17]])
B2_MG_Xylem <- read.table(snakemake@input[[22]])
B2_PR_Xylem <- read.table(snakemake@input[[27]])

## E5 
E5_MG_and_PR_Leaf <- read.table(snakemake@input[[3]])
E5_MG_Leaf <- read.table(snakemake@input[[8]])
E5_PR_Leaf <- read.table(snakemake@input[[13]])
E5_MG_and_PR_Xylem <- read.table(snakemake@input[[18]])
E5_MG_Xylem <- read.table(snakemake@input[[23]])
E5_PR_Xylem <- read.table(snakemake@input[[28]])

## E5
Q8_MG_and_PR_Leaf <- read.table(snakemake@input[[4]])
Q8_MG_Leaf <- read.table(snakemake@input[[9]])
Q8_PR_Leaf <- read.table(snakemake@input[[14]])
Q8_MG_and_PR_Xylem <- read.table(snakemake@input[[19]])
Q8_MG_Xylem <- read.table(snakemake@input[[24]])
Q8_PR_Xylem <- read.table(snakemake@input[[29]])

## Q8
D4_MG_and_PR_Leaf <- read.table(snakemake@input[[5]])
D4_MG_Leaf <- read.table(snakemake@input[[10]])
D4_PR_Leaf <- read.table(snakemake@input[[15]])
D4_MG_and_PR_Xylem <- read.table(snakemake@input[[20]])
D4_MG_Xylem <- read.table(snakemake@input[[25]])
D4_PR_Xylem <- read.table(snakemake@input[[30]])

#Defines the levels of the analysis
clones <- c("A1", "B2", "D4", "E5", "Q8")
local <- c("MG", "MG_and_PR", "PR")
tissue <- c("Leaf", "Xylem")

# Builds the table with cumulative values for bar plot.
graph_table <- data.frame()
for(c in 1:length(clones)){
    for(t in 1:length(tissue)){
        for(l in 1:length(local)){
            if(l == 1){
                l1 <- get(paste(clones[c], local[l], tissue[t], sep="_"))
                y <- as.data.frame(t(c(clones[c], tissue[t], local[l], length(l1$V1), length(l1$V1), length(l1$V1)/2)))
                names(y) <- c("nomes", "tissue", "local", "valor", "label", "posicao")
                graph_table <- rbind(graph_table, y)
            }else if(l == 2){
                l2 <- get(paste(clones[c], local[l], tissue[t], sep="_"))
                y <- as.data.frame(t(c(clones[c], tissue[t], local[l], (length(l1$V1)+ length(l2$V1)), length(l2$V1), (length(l2$V1)/2+length(l1$V1)))))
                names(y) <- c("nomes", "tissue", "local", "valor", "label", "posicao")
                graph_table <- rbind(graph_table, y)            
            }else if(l == 3){
                l3 <- get(paste(clones[c], local[l], tissue[t], sep="_"))
                y <- as.data.frame(t(c(clones[c], tissue[t], local[l], length(l1$V1)+ length(l2$V1) + length(l3$V1),length(l3$V1), (length(l3$V1)/2+(length(l1$V1)+length(l2$V1))))))
                names(y) <- c("nomes", "tissue", "local", "valor", "label", "posicao")
                graph_table <- rbind(graph_table, y)
            }
        }
    }
}

# Change factors to numeric or integer accordingly with the variable. 
graph_table$local <- factor(graph_table$local, levels = c("PR", "MG_and_PR", "MG"))
graph_table$posicao <- as.integer(as.character(graph_table$posicao))
graph_table$label <- as.numeric(as.character(graph_table$label))
graph_table$valor <- as.numeric(as.character(graph_table$valor))

# Builds and save the bar plot.
svg(filename=snakemake@output[[1]], width=12, height=8)

## Defines the colors of the bars.
fill <- c("#F0E442", "#56B4E9", "#009E73")

ggplot(graph_table, aes(x = nomes, y = label, fill = local)) + 
    geom_bar(stat = "identity", alpha = 0.8, width = 0.7) + 
    geom_text(size = 6, aes(y = posicao, label = label), colour = "black") +
    facet_wrap(~tissue, ncol = 2) +
    scale_x_discrete(name = "") +
    ylab("Number of methylated sites") +
    scale_y_continuous(breaks = seq(0, 10000, 500)) +
    theme_bw() +
    theme(legend.text = element_text(size = 14),
          axis.title.y = element_text(size = 16, vjust = 2),
          axis.title.x = element_text(size = 14, vjust = 0),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 15)) +
    scale_fill_manual(values=fill,
                      labels = c("Wet site (PR)", "Both sites", "Dry site (MG)"))

dev.off()

