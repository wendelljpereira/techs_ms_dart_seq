### Comparação dos dados entre os dois tecidos ###
save.image("snakemake_venn_diag_msdart_sites.rda")

require(plyr)
require(VennDiagram)
require(gridExtra)

# Define a função que criará o gráfico
venn_2_samples <- function(sample1 = sample1, sample2 = sample2, name1 = name1, name2 = name2, clone_name = clone_name, save_ids = "FALSE"){
    
    # Determina os parâmetros para a criação do Diagrama de Venn
    inter <- intersect(sample1[complete.cases(sample1)], sample2[complete.cases(sample2)]) # interseção
    unic_s1 <- setdiff(sample1[complete.cases(sample1)], sample2[complete.cases(sample2)]) # exclusivos no conjunto do primeiro argumento. 
    unic_s2 <- setdiff(sample2[complete.cases(sample2)], sample1[complete.cases(sample1)]) # exclusivos no conjunto do segundo argumento.
    
    # Escreve os arquivos com os Ids pertencentes a cada grupo do diagrama de venn.
    if(save_ids == "TRUE"){
        write.table(inter, paste0("intersections_to_bar_plot/", paste(clone_name, name1, "vs", name2, "intersection_marks.txt", sep = "_")), row.names = F, col.names = F, quote = F)
        write.table(unic_s1, paste0("intersections_to_bar_plot/", paste(clone_name, name1, "vs", name2, name1, "unique_marks.txt", sep = "_")), row.names = F, col.names = F, quote = F)
        write.table(unic_s2, paste0("intersections_to_bar_plot/", paste(clone_name, name1, "vs", name2, name2, "unique_marks.txt", sep = "_")), row.names = F, col.names = F, quote = F)
    }else if(save_ids == "FALSE"){
        
    }else{
        print("Invalid save_ids option!")
    }
    
    # Constrói o Diagrama de Venn a partir dos conjuntos gerados acima.
    if(length(unic_s1) > length(unic_s2)){
        grid.newpage();
        venn.plot <- draw.pairwise.venn(
            area1 = length(unic_s1) + length(inter),
            area2 = length(unic_s2) + length(inter),
            cross.area = length(inter),
            alpha = 0.65,
            category = c(deparse(name1), deparse(name2)),
            fill = c("#8470FF", "#FF7F24"),
            lty = "blank",
            cex = 3,
            cat.cex = 3,
            cat.pos = c(0, 0),
            cat.dist = 0.055,
            cat.just = list(c(1, 0), c(0,0 )),
            ext.pos = 0,
            ext.dist = -0.05,
            ext.length = 0.85,
            ext.line.lwd = 2,
            ext.line.lty = "dashed",
            scaled = T,
            print.mode = c("raw", "percent"),
            rotation.degree = 0)
    }else{
        grid.newpage();
        venn.plot <- draw.pairwise.venn(
            area1 = length(unic_s1) + length(inter),
            area2 = length(unic_s2) + length(inter),
            cross.area = length(inter),
            alpha = 0.65,
            category = c(deparse(name1), deparse(name2)),
            fill = c("#8470FF", "#FF7F24"),
            lty = "blank",
            cex = 3,
            cat.cex = 3,
            cat.pos = c(0, 0),
            cat.dist = 0.055,
            cat.just = list(c(1, 0), c(0,0 )),
            ext.pos = 0,
            ext.dist = -0.05,
            ext.length = 0.85,
            ext.line.lwd = 2,
            ext.line.lty = "dashed",
            scaled = T,
            print.mode = c("raw", "percent"),
            rotation.degree = 180)
    }
}

## Análise com as marcas diferencialmente metiladas ##

# Carregando os arquivos com as marcas "diferencialmente expressas" para aplica-los na função
G2_int <- read.table(snakemake@input[[1]], header = T, na.strings = "NA", colClasses = "character")
G3_int <- read.table(snakemake@input[[2]], header = T, na.strings = "NA", colClasses = "character")

## A1
C <- venn_2_samples(G2_int$A1_leaf, G3_int$A1_leaf,"Leaf-MG","Leaf-PR", "A1", save_ids = "TRUE")
D <- venn_2_samples(G2_int$A1_wood, G3_int$A1_wood,"Wood-MG","Wood-PR", "A1", save_ids = "TRUE")

## Salvando como svg para edição posterior
svg(filename = snakemake@output[[1]], width = 20, height = 10, pointsize = 12)
grid.arrange(grobTree(C), grobTree(D), ncol=2, top=textGrob("A1", gp=gpar(fontsize=40,font=8)))
dev.off()

## D4
C <- venn_2_samples(G2_int$D4_leaf, G3_int$D4_leaf,"Leaf-MG","Leaf-PR","D4", save_ids = "TRUE")
D <- venn_2_samples(G2_int$D4_wood, G3_int$D4_wood,"Wood-MG","Wood-PR","D4", save_ids = "TRUE")

## Salvando como svg para edição posterior
svg(filename = snakemake@output[[3]], width = 20, height = 10, pointsize = 12)
grid.arrange(grobTree(C), grobTree(D), ncol=2, top=textGrob("D4", gp=gpar(fontsize=40,font=8)))
dev.off()

## Q8
C <- venn_2_samples(G2_int$Q8_leaf, G3_int$Q8_leaf,"Leaf-MG","Leaf-PR","Q8", save_ids = "TRUE")
D <- venn_2_samples(G2_int$Q8_wood, G3_int$Q8_wood,"Wood-MG","Wood-PR","Q8", save_ids = "TRUE")

## Salvando como svg para edição posterior
svg(filename = snakemake@output[[5]], width = 20, height = 10, pointsize = 12)
grid.arrange(grobTree(C), grobTree(D), ncol=2, top=textGrob("Q8", gp=gpar(fontsize=40,font=8)))
dev.off()

## E5
C <- venn_2_samples(G2_int$E5_leaf, G3_int$E5_leaf,"Leaf-MG","Leaf-PR","E5", save_ids = "TRUE")
D <- venn_2_samples(G2_int$E5_wood, G3_int$E5_wood,"Wood-MG","Wood-PR","E5", save_ids = "TRUE")

## Salvando como svg para edição posterior
svg(filename = snakemake@output[[4]], width = 20, height = 10, pointsize = 12)
grid.arrange(grobTree(C), grobTree(D), ncol=2, top=textGrob("E5", gp=gpar(fontsize=40,font=8)))
dev.off()

## B2
C <- venn_2_samples(G2_int$B2_leaf, G3_int$B2_leaf,"Leaf-MG","Leaf-PR", "B2", save_ids = "TRUE")
D <- venn_2_samples(G2_int$B2_wood, G3_int$B2_wood,"Wood-MG","Wood-PR", "B2", save_ids = "TRUE")

## Salvando como svg para edição posterior
svg(filename = snakemake@output[[2]], width = 20, height = 10, pointsize = 12)
grid.arrange(grobTree(C), grobTree(D), ncol=2, top=textGrob("B2", gp=gpar(fontsize=40,font=8)))
dev.off()
