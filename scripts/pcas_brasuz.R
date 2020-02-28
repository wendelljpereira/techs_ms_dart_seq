save.image("PCAs.Rda")
require(edgeR)
require(stringr)
require(ggfortify)
require(gridExtra)
require(VennDiagram)
require(ggthemes)
require(tidyverse)

###########################################################################################
############################### Using only methylated sites ###############################
###########################################################################################

# All samples, but only DM marks
counts_data <- read.table(snakemake@input[[1]],
                          header = T,
                          sep = ",",
                          check.names = F,
                          colClasses = c("character", rep("numeric", 80)))

# File with information about DM marks in BRASUZ1 samples.
DM_sites <- read_tsv(snakemake@input[[2]],
                     col_names =  F) %>%
    select(X4)

DM_sites_names <- character()
for(i in 1:nrow(DM_sites)){
    
    str_sub <- DM_sites[i, 1]
    
    new_char <- str_split(str_sub, "\\|") %>% 
        unlist()
    
    DM_sites_names <- c(DM_sites_names, new_char)
}

counts_data <- counts_data[counts_data$Geneid %in% DM_sites_names, ]
counts_data_filt <- counts_data[, -1]

##########################################################
######## PCAs comparing sites, clones and tissues ########
##########################################################

## For each tissue:

tissues <- c("leaf", "wood")
for( i in 1:length(tissues)){
    
    ## Plots including both enzymes##
    
    # Filter the data of each tissue
    counts_data_sub <- counts_data[, grep(tissues[i], x = colnames(counts_data)) ]
    counts_data_t <- as.data.frame(t(counts_data_sub))
    
    # Create the index to color the plot
    clones_names <- as.data.frame(row.names(counts_data_t))
    
    colnames(clones_names) <- "clone"
    clones_names <- as.data.frame(str_split_fixed(clones_names$clone, "_", 4))
    colnames(clones_names) <- c("group", "clone", "tissue", "enzime")
    
    clones_names$clone_group <- paste(clones_names$clone,
                                      clones_names$group,
                                      sep = "_")
    
    clones_names$enzime <- substr(clones_names$enzime, 1, 2)
    
    # Remove the marks with counts equal to zero for all samples
    subset_index <- as.data.frame(lapply(counts_data_t, function(x) sum(x) > 0))
    counts <- counts_data_t[, subset_index == TRUE]
    
    # Normalize the data and execute the PCA
    res.pca <- prcomp(counts, scale = TRUE, center = TRUE, retx = T)
    
    # Store the data with the variation
    x <- summary(res.pca)
    
    col_blind_pal <- c("#D81B60",
                       "#1E88E5",
                       "#FFC107",
                       "#004D40",
                       "#401F18",
                       "#EA82CE",
                       "#7EC4CA",
                       "#B51DB7",
                       "#31E2E5",
                       "#9F9E8B")
    
    # PCA colored by clone/grupo (With frame)
    A <- autoplot(res.pca,
                  data = clones_names,
                  colour = "clone_group",
                  shape = "enzime",
                  size = 3,
                  frame = T,
                  frame.type = "norm") +
        theme(legend.text = element_text(size = 14),
              axis.title.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.title = element_blank(),
              strip.text = element_text(size = 15)) +
        scale_colour_manual(values = col_blind_pal) +
        scale_fill_manual(values = col_blind_pal) +
        theme_bw()
    
    B <- autoplot(res.pca,
                  data = clones_names,
                  colour = "clone_group",
                  shape = "enzime",
                  size = 3,
                  frame = F,
                  frame.type = "norm") +
        theme(legend.text = element_text(size = 14),
              axis.title.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.title = element_blank(),
              strip.text = element_text(size = 15)) +
        scale_colour_manual(values = col_blind_pal) +
        scale_fill_manual(values = col_blind_pal) +
        theme_bw()
    
    svg(filename = paste0("images/PCAs/", paste0("PCA_by_", tissues[i], "_with_frame.svg")), width = 12, height = 12, pointsize = 12)
    print(A)
    dev.off()
    
    svg(filename = paste0("images/PCAs/", paste0("PCA_by_", tissues[i], "_no_frame.svg")), width = 12, height = 12, pointsize = 12)
    print(B)
    dev.off()
    
    ### Using only mspI
    
    counts_data_sub <- counts_data_sub[, grep("ms", x = colnames(counts_data_sub)) ]
    
    counts_data_t <- as.data.frame(t(counts_data_sub))
    
    # Create the index to color the plot
    clones_names <- as.data.frame(row.names(counts_data_t))
    
    colnames(clones_names) <- "clone"
    clones_names <- as.data.frame(str_split_fixed(clones_names$clone, "_", 4))
    colnames(clones_names) <- c("group", "clone", "tissue", "enzime")
    
    clones_names$clone_group <- paste(clones_names$clone,
                                      clones_names$group,
                                      sep = "_")
    
    clones_names$enzime <- substr(clones_names$enzime, 1, 2)
    
    # Remove the marks with counts equal to zero for all samples
    subset_index <- as.data.frame(lapply(counts_data_t, function(x) sum(x) > 0))
    counts <- counts_data_t[, subset_index == TRUE]
    
    # Normalize the data and execute the PCA
    res.pca <- prcomp(counts, scale = TRUE, center = TRUE, retx = T)
    
    # Store the data with the variation
    x <- summary(res.pca)
    
    col_blind_pal <- c("#D81B60",
                       "#1E88E5",
                       "#FFC107",
                       "#004D40",
                       "#401F18",
                       "#EA82CE",
                       "#7EC4CA",
                       "#B51DB7",
                       "#31E2E5",
                       "#9F9E8B")
    
    C <- autoplot(res.pca,
                  data = clones_names,
                  colour = "clone_group",
                  shape = "enzime",
                  size = 3,
                  frame = F,
                  frame.type = "norm") +
        theme(legend.text = element_text(size = 14),
              axis.title.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.title = element_blank(),
              strip.text = element_text(size = 15)) +
        scale_colour_manual(values = col_blind_pal) +
        scale_fill_manual(values = col_blind_pal) +
        theme_bw()
    
    svg(filename = paste0("images/PCAs/", paste0("PAC_by_", tissues[i], "_only_mspI_no_frame.svg")), width = 12, height = 12, pointsize = 12)
    print(C)
    dev.off()
}

## For each clone compare how the tissues separate in each environment
clone_names <- c("A1", "B2", "D4", "E5", "Q8")
for( i in 1:length(clone_names)){
    
    counts_data_sub_both_t <- counts_data[, grep(clone_names[i], x = colnames(counts_data)) ]
    
    for( t in 1:length(tissues) ){
        
        ## Plots including both enzymes##
        
        # Filter the data of each tissue
        counts_data_sub <- counts_data_sub_both_t[, grep(tissues[t], x = colnames(counts_data_sub_both_t)) ]
        counts_data_t <- as.data.frame(t(counts_data_sub))
        
        # Create the index to color the plot
        clones_names <- as.data.frame(row.names(counts_data_t))
        
        colnames(clones_names) <- "clone"
        clones_names <- as.data.frame(str_split_fixed(clones_names$clone, "_", 4))
        colnames(clones_names) <- c("group", "clone", "tissue", "enzime")
        
        clones_names$enzime <- substr(clones_names$enzime, 1, 2)
        
        # Remove the marks with counts equal to zero for all samples
        subset_index <- as.data.frame(lapply(counts_data_t, function(x) sum(x) > 0))
        counts <- counts_data_t[, subset_index == TRUE]
        
        # Normalize the data and execute the PCA
        res.pca <- prcomp(counts, scale = TRUE, center = TRUE, retx = T)
        
        # Store the data with the variation
        x <- summary(res.pca)
        
        # PCA colored by clone/grupo (With frame)
        A <- autoplot(res.pca,
                      data = clones_names,
                      colour = "group",
                      shape = "enzime",
                      size = 3,
                      frame = T,
                      frame.type = "norm") +
            theme(legend.text = element_text(size = 14),
                  axis.title.y = element_text(size = 15),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14),
                  legend.title = element_blank(),
                  strip.text = element_text(size = 15)) +
            scale_color_colorblind() +
            scale_fill_colorblind() +
            theme_bw()
        
        B <- autoplot(res.pca,
                      data = clones_names,
                      colour = "group",
                      shape = "enzime",
                      size = 3,
                      frame = F,
                      frame.type = "norm") +
            theme(legend.text = element_text(size = 14),
                  axis.title.y = element_text(size = 15),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14),
                  legend.title = element_blank(),
                  strip.text = element_text(size = 15)) +
            scale_color_colorblind() +
            scale_fill_colorblind() +
            theme_bw()
        
        svg(filename = paste0("images/PCAs/", paste("PCA_by", clone_names[i], tissues[t], "with_frame.svg", sep = "_")), width = 12, height = 12, pointsize = 12)
        print(A)
        dev.off()
        
        svg(filename = paste0("images/PCAs/", paste("PCA_by", clone_names[i], tissues[t], "no_frame.svg", sep = "_")), width = 12, height = 12, pointsize = 12)
        print(B)
        dev.off()
        
        ### Using only mspI
        
        counts_data_sub <- counts_data_sub[, grep("ms", x = colnames(counts_data_sub)) ]
        
        counts_data_t <- as.data.frame(t(counts_data_sub))
        
        # Create the index to color the plot
        clones_names <- as.data.frame(row.names(counts_data_t))
        
        colnames(clones_names) <- "clone"
        clones_names <- as.data.frame(str_split_fixed(clones_names$clone, "_", 4))
        colnames(clones_names) <- c("group", "clone", "tissue", "enzime")
        
        clones_names$clone_group <- paste(clones_names$clone,
                                          clones_names$group,
                                          sep = "_")
        
        clones_names$enzime <- substr(clones_names$enzime, 1, 2)
        
        # Remove the marks with counts equal to zero for all samples
        subset_index <- as.data.frame(lapply(counts_data_t, function(x) sum(x) > 0))
        counts <- counts_data_t[, subset_index == TRUE]
        
        # Normalize the data and execute the PCA
        res.pca <- prcomp(counts, scale = F, center = TRUE, retx = T)
        
        # Store the data with the variation
        x <- summary(res.pca)
        
        # PCA colored by clone/grupo (With frame)
        C <- autoplot(res.pca,
                      data = clones_names,
                      colour = "clone_group",
                      shape = "enzime",
                      size = 3,
                      frame = F,
                      frame.type = "norm") +
            theme(legend.text = element_text(size = 14),
                  axis.title.y = element_text(size = 15),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14),
                  legend.title = element_blank(),
                  strip.text = element_text(size = 15)) +
            scale_colour_manual(values = col_blind_pal) +
            scale_fill_manual(values = col_blind_pal) +
            theme_bw()
        
        svg(filename = paste0("images/PCAs/", paste("PCA_by", clone_names[i], tissues[t], "only_mspI_no_frame.svg", sep = "_")), width = 12, height = 12, pointsize = 12)
        print(C)
        dev.off()
    }
}
