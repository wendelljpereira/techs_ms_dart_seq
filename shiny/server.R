if (!is.element("org.Egrandis.eg.db", installed.packages()[, 1])){
    
    install.packages("data/org.Egrandis.eg.db",
                     repos = NULL,
                     type = "source")
} else {
    
    require("org.Egrandis.eg.db")
    
}

if (!is.element("lfmm", installed.packages()[, 1])){
    
    if (!is.element("devtools", installed.packages()[, 1])){
        
        install.packages("devtools", dependencies = T)
    }
    
    devtools::install_github("bcm-uga/lfmm")
}

use_package_cran <- function(p){
    
    if (!is.element(p, installed.packages()[, 1])){
        
        install.packages(p, dep = TRUE)
        require(p, character.only = TRUE)
        
    } else {
        
        require(p, character.only = TRUE)
        
    }
}

use_package_bioc <- function(p) {
    
    if (!is.element(p, installed.packages()[, 1])) {
        
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        
        BiocManager::install(p)
        
        require(p, character.only = TRUE)
        
    } else {
        
        require(p, character.only = TRUE)
        
    }
}

# list of necessary packages to run the application
packages_cran <- c(
    "factoextra",
    "FactoMineR",
    "gdata",
    "ggthemes",
    "lfmm",
    "qqman",
    "circlize",
    "RColorBrewer",
    "scales",
    "shiny",
    "tidyverse",
    "vroom",
    "reactable",
    "VennDiagram",
    "gplots",
    "rclipboard",
    "reshape2",
    "purrr",
    "cluster",
    "fpc")

packages_bioc <- c("ComplexHeatmap",
                   "clusterProfiler")

# installs and/or load the packages from CRAN
for (i in packages_cran){
    use_package_cran(i)
}

for (i in packages_bioc){
    use_package_bioc(i)
}

## Loading phenotype data -- Tree level##
clone_tech_names <- c("A1", "B2", "D4", "E5", "Q8")

## Loading phenotype data -- Plot level##
prod_g2 <- read_csv("data/TECHS_Data_Plot_Level_19_06_2019.csv") %>%
    dplyr::filter(date_for_measurement == "10/1/15") %>%
    dplyr::filter(Measurement_date == "9/24/15") %>%
    dplyr::filter(Throughfall_exclusion == "não") %>%
    dplyr::group_by(Clone) %>%
    dplyr::summarise(Basal_area_mean = mean(Basal_area),
                     Volume_mean = mean(Volume),
                     Stem_dry_mass_mean = mean(Stem_dry_mass))

prod_g3 <- read_csv("data/TECHS_Data_Plot_Level_19_06_2019.csv",
                    col_types = cols(Site = "f",
                                     Plot = "f")) %>%
    dplyr::select(-c("PV50", "Missing_tree", "Dead_tree")) %>%
    dplyr::filter(date_for_measurement == "10/1/15") %>%
    dplyr::filter(Measurement_date == "10/15/15") %>%
    dplyr::filter(Throughfall_exclusion == "não") %>%
    dplyr::group_by(Clone) %>%
    dplyr::summarise(Basal_area_mean = mean(Basal_area),
                     Volume_mean = mean(Volume),
                     Stem_dry_mass_mean = mean(Stem_dry_mass))

# Reads the data files with Differential Methylated marks
g2_int <- read.table("data/deseq_group_2_consensus_methylated_sites_intersect.tst",
                     header = T,
                     sep = "\t",
                     colClasses = "character")

g3_int <- read.table("data/deseq_group_3_consensus_methylated_sites_intersect.tst",
                     header = T,
                     sep = "\t",
                     colClasses = "character")

flmm_clones <- c("A1",
                 "B2",
                 "D4",
                 "E5",
                 "Q8")

shinyServer(function(input, output) {
    
    ###### Tab 3 - LFMM analysis ######
    
    # Creates a panel that reacts to the clones choices to allows the user to select the comparisons only among the clones already choosed 
    output$comps_to_annot_ui = renderUI({
        
        div(class = "option-group",
            checkboxGroupInput("comps_to_annot",
                               "Select the clones participating the comparison that you want to annotate", 
                               choices = input$annot_clones,
                               selected = input$annot_clones
            ))
        
        
    })
    
    # Help functions 
    help_manhat_plot <- function(input_df = input_df, trait = trait, plot_title = plot_title) {
        
        pvalues_df_chr_trait <- input_df %>%
            dplyr::select("MSD-methylated",
                          paste0("p-value.", trait),
                          "Chr",
                          "start") %>%
            dplyr::rename(pvalue = paste0("p-value.", trait)) %>%
            dplyr::filter(!is.nan(pvalue))
        
        trait_sig_df <- pvalues_df_chr_trait %>%
            dplyr::filter(pvalue < input$lfmm_pvalue)
        
        trait_sig_names <- trait_sig_df$`MSD-methylated`
        
        manhattan(pvalues_df_chr_trait,
                  chr="Chr",
                  bp="start",
                  snp="MSD-methylated",
                  p="pvalue",
                  main = plot_title,
                  col = c("blue4", "orange3"),
                  suggestiveline = F,
                  highlight = trait_sig_names,
                  genomewideline = -log10(input$lfmm_pvalue),
                  ylim = c(0, round(max(-log10(pvalues_df_chr_trait$pvalue)) + 1)))
    }
    
    help_bar_plot <- function(input_df = input_df,
                              epigenotype_df = epigenotype_df,
                              trait = trait,
                              plot_title = plot_title){
        
        pvalues_df_chr_trait <- input_df %>%
            dplyr::select("MSD-methylated",
                          paste0("p-value.", trait),
                          "Chr",
                          "start") %>%
            dplyr::rename(pvalue = paste0("p-value.", trait)) %>%
            dplyr::filter(!is.nan(pvalue))
        
        sig_pvalues <- pvalues_df_chr_trait[pvalues_df_chr_trait$pvalue < input$lfmm_pvalue, ]
        sig_lfmm_MSD <- sig_pvalues$`MSD-methylated`
        
        epigenotype_sig <- epigenotype_df %>%
            dplyr::mutate(clones = paste(epigenotype_df$clones,
                                         epigenotype_df$site,
                                         sep = "_")) %>%
            dplyr::select(- site) %>%
            t()
        
        colnames(epigenotype_sig) <- epigenotype_sig[1, ]
        epigenotype_sig <- as.data.frame(epigenotype_sig)
        epigenotype_sig <- epigenotype_sig[rownames(epigenotype_sig) %in% sig_lfmm_MSD, ]
        
        clone_names <- c("A1",
                         "B2",
                         "D4",
                         "E5",
                         "Q8")
        
        comp_btw_sites <- data.frame()
        for( c in clone_names){
            
            clones_set <- epigenotype_sig[grep(c, colnames(epigenotype_sig))]
            clones_set[, 1] <- as.numeric(as.character(clones_set[, 1]))
            clones_set[, 2] <- as.numeric(as.character(clones_set[, 2]))
            
            prop <- ifelse(clones_set[, 1] == clones_set[, 2] & clones_set[, 1] == 0, "Unmethylated in both localities", 
                           ifelse(clones_set[, 1] == clones_set[, 2], "Methylated in both localities",
                                  ifelse(clones_set[, 1] == 1, "Methylated only in Wet (PR)",
                                         "Methylated only in Dry (MG)")))
            
            if (nrow(comp_btw_sites) == 0){
                
                comp_btw_sites <- as.data.frame(prop)
                colnames(comp_btw_sites) <- c
                
            } else {
                
                comp_btw_sites_sub <- as.data.frame(prop)
                colnames(comp_btw_sites_sub) <- c
                comp_btw_sites <- cbind(comp_btw_sites, comp_btw_sites_sub)
                
            }
        }
        
        ggplot_input <- as.data.frame(summary(comp_btw_sites))
        
        # Removes multiple spaces
        ggplot_input$Freq <- str_replace(gsub("\\s+", " ", str_trim(ggplot_input$Freq)), "B", "b")
        ggplot_input$Freq <- gsub(" :", ":", ggplot_input$Freq)
        
        ggplot_input2 <- strsplit(as.character(ggplot_input$Freq), split = ":")
        ggplot_input2 <- as.data.frame(do.call(rbind, ggplot_input2))
        
        ggplot_input3 <- data.frame("Clone" = ggplot_input$Var2,
                                    "Class" = ggplot_input2$V1,
                                    "Counts" = as.numeric(as.character(ggplot_input2$V2)))
        
        ggplot(ggplot_input3, aes(x = Clone, y = Counts, fill = Class)) +
            geom_bar(position = "fill", stat = "identity", alpha = 0.9, width = 0.5) +
            scale_y_continuous(labels = percent_format())+
            coord_flip() + 
            ggtitle(plot_title) +
            ylab("Relative frequency") +
            xlab("") +
            theme_bw() +
            theme(legend.position="bottom")+
            scale_fill_colorblind()
        
    }
    
    help_tables <- function(input_df = input_df,
                            trait = trait,
                            epigenotype_df = epigenotype_df){
        
        # Selects the sites that has a pvalue for the trait
        pvalues_df_chr_trait <- input_df %>%
            dplyr::select("MSD-methylated",
                          paste0("p-value.", trait),
                          "Chr",
                          "start") %>%
            dplyr::rename(pvalue = paste0("p-value.", trait)) %>%
            dplyr::filter(!is.nan(pvalue))
        
        # Selects the significative p-values
        sig_pvalues <- pvalues_df_chr_trait[pvalues_df_chr_trait$pvalue < input$lfmm_pvalue, ]
        sig_lfmm_MSD <- sig_pvalues$`MSD-methylated`
        
        # Selects the epigenotype of this sites
        epigenotype_sig <- epigenotype_df %>%
            dplyr::mutate(clones = paste(epigenotype_df$clones,
                                         epigenotype_df$site,
                                         sep = "_")) %>%
            dplyr::select(- site) %>%
            t()
        
        colnames(epigenotype_sig) <- epigenotype_sig[1, ]
        epigenotype_sig <- as.data.frame(epigenotype_sig)
        epigenotype_sig <- epigenotype_sig[rownames(epigenotype_sig) %in% sig_lfmm_MSD, ]
        
        clone_names <- c("A1",
                         "B2",
                         "D4",
                         "E5",
                         "Q8")
        
        meth_status_by_locality <- data.frame()
        for( c in clone_names) {
            
            clones_set <- epigenotype_sig[grep(c, colnames(epigenotype_sig))]
            clones_set[, 1] <- as.numeric(as.character(clones_set[, 1]))
            clones_set[, 2] <- as.numeric(as.character(clones_set[, 2]))
            
            clones_set$clone <- paste(c)
            
            clones_set$class <- ifelse(clones_set[, 1] == clones_set[, 2] & clones_set[, 1] == 0,
                                       "Unmethylated in both localities", 
                                       ifelse(clones_set[, 1] == clones_set[, 2],
                                              "Methylated in both localities",
                                              ifelse(clones_set[, 1] == 1,
                                                     "Methylated only in Wet (PR)",
                                                     "Methylated only in Dry (MG)")))
            
            clones_set[, 1] <- rownames(clones_set)
            clones_set <- clones_set[,-2]
            
            colnames(clones_set) <- c("MSD-site", "clone", "class")
            meth_status_by_locality <- rbind(meth_status_by_locality, clones_set)
        }
        
        colnames(sig_pvalues)[1] <- "MSD-site"
        sig_pvalues2 <- merge(sig_pvalues,
                              meth_status_by_locality,
                              by = "MSD-site")
    }
    
    # Loads and processes the trait measures
    phenotype <- reactive({
        
        phenot_mean <- vroom("data/TECHS_Data_Tree_Level_19_06_2019.csv") %>%
            dplyr::filter(date_for_measurement == "10/1/15") %>%
            dplyr::filter(Clone %in% clone_tech_names) %>%
            dplyr::filter(Throughfall_exclusion == "Não") %>%
            dplyr::mutate(location = ifelse(Site == 30, "Dry", "Wet")) %>%
            dplyr::mutate(clone_site = paste(Clone, Site, sep = "_")) %>%
            dplyr::filter(!is.na(Volume)) %>%
            dplyr::filter(!is.na(Stem_dry_mass)) %>%
            dplyr::filter(!is.na(Total_height)) %>%
            dplyr::filter(!is.na("DBH")) %>%
            dplyr::filter(!is.na(Height_bottom_live_crown)) %>%
            dplyr::group_by(Clone, Site) %>%
            dplyr::summarise_if(funs(mean), .predicate = is.numeric) %>%
            dplyr::select(c("Clone", "Site", input$traits )) %>%
            dplyr::mutate(Site = as.character(Site))
        
        phenotype <- phenot_mean[phenot_mean$Clone %in% flmm_clones, ]
        
        phenotype <- phenotype[order(phenotype$Clone, phenotype$Site ), ]
        phenotype <- phenotype[, -c(1:2)]
        
        #    Scales and center the phenotypes
        if( input$lfmm_phenot_cent == TRUE) {
            
            phenotype <- scale(phenotype, center = T)
            
        }
        
    })
    
    epigenotype <- reactive({
        
        ## Reads the epigenotype of clones in each location
        if (input$flmm_tissue == "leaves") {
            
            epigenotype <- vroom("data/leaves_epigenotype_techs.tst")
            
        } else if (input$flmm_tissue == "xylem") {
            
            epigenotype <- vroom("data/xylem_epigenotype_techs.tst")
            
        }
        
        # Selects the clones
        epigenotype <- epigenotype[epigenotype$clones %in% flmm_clones, ]
        
        # Check if should use all MSD-Sites or only the ones with methylation in at least one sample. Also, order both datasets
        if (input$flmm_marks_selection == "all") {
            
            epigenotype <- epigenotype[order(epigenotype$clones, epigenotype$site ), ]
            
        } else if (input$flmm_marks_selection == "methy_true"){
            
            epigenotype <- epigenotype[order(epigenotype$clones, epigenotype$site ), ]
            epigenotype_id <- epigenotype[, c(1,2)]
            epigenotype <- epigenotype[, -c(1:2)]
            epigenotype <- as.data.frame(epigenotype[,  colSums(epigenotype) > 0])
            
            epigenotype <- cbind(epigenotype_id, epigenotype)
        }
        
        epigenotype
    })
    
    output$scree_plot <- renderPlot({
        
        pc <- prcomp(epigenotype()[, -c(1,2)])
        fviz_eig(pc, addlabels = TRUE)
        
    })
    
    pvalues_df_chr <- reactive({
        
        ## Fit an LFMM, i.e, compute B, U, V estimates
        mod.lfmm <- lfmm_ridge(Y = epigenotype()[, -c(1,2)],
                               X = phenotype(),
                               K = input$flmm_k)
        
        ## performs association testing using the fitted model:
        pv <- lfmm_test(Y = epigenotype()[, -c(1,2)],
                        X = phenotype(),
                        lfmm = mod.lfmm,
                        calibrate = "gif")
        
        pvalues_df <- data.frame("MSD-methylated" = rownames(pv$calibrated.pvalue) ,
                                 "p-value" = pv$calibrated.pvalue,
                                 check.names = FALSE)
        
        position <- read_tsv("data/msdartseq_methylation_sites_of_sequenced_fragments_merged.bed",
                             col_names = F) %>%
            dplyr::filter(X6 == "+") %>%
            dplyr::rename("Chr" = X1,
                          "start" = X2,
                          "end" = X3,
                          "MSD-methylated" = X4) %>%
            dplyr::select(-c(X5, X6))
        
        pvalues_df <- merge(pvalues_df, position, by = "MSD-methylated")
        pvalues_df
        # Removes marks in scaffolds
        pvalues_df_chr <- pvalues_df[grep(pattern = "Chr", x = as.character(pvalues_df$Chr)), ]
        pvalues_df_scaffold <- pvalues_df[-grep(pattern = "Chr",
                                                x = as.character(pvalues_df$Chr)), ]
        
        pvalues_df_chr$Chr <- gsub(pattern = "Chr0",replacement = "", x =  pvalues_df_chr$Chr)
        pvalues_df_chr$Chr <- gsub(pattern = "Chr",replacement = "", x =  pvalues_df_chr$Chr)
        pvalues_df_chr$Chr <- as.numeric(pvalues_df_chr$Chr)
        
        pvalues_df_chr <- pvalues_df_chr %>% 
            dplyr::filter_all(all_vars(!is.infinite(.)))
        
        pvalues_df_chr
    })
    
    ## Generates the manhattan plots
    
    observe({
        
        if ("Volume" %in% input$traits) {
            
            output$manhattan_vol <- renderPlot({
                
                help_manhat_plot(input_df = pvalues_df_chr(),
                                 trait = "Volume",
                                 plot_title = "MSD-Methylated sites associated with volume")
                
            })
            output$bar_vol <- renderPlot({
                
                help_bar_plot(input_df = pvalues_df_chr(),
                              epigenotype_df = epigenotype(),
                              trait = "Volume",
                              plot_title = "Methylation status of the significative marks - Volume")
                
            })
            output$table_vol <- renderReactable({
                
                table_vol <- help_tables(input_df = pvalues_df_chr(),
                                         trait = "Volume",
                                         epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_vol) > 0 ) {
                    
                    table_vol <- table_vol[table_vol$clone %in% input$flmm_res_filter_by_clone_vol, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_vol) > 0 ) {
                    
                    table_vol <- table_vol[table_vol$class %in% input$flmm_res_filter_by_class_vol, ]
                    
                }
                
                reactable(table_vol,
                          defaultColDef = colDef(align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 5,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
                
            })
            output$clip_volume <- renderUI({
                
                table_vol <- help_tables(input_df = pvalues_df_chr(),
                                         trait = "Volume",
                                         epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_vol) > 0 ) {
                    
                    table_vol <- table_vol[table_vol$clone %in% input$flmm_res_filter_by_clone_vol, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_vol) > 0 ) {
                    
                    table_vol <- table_vol[table_vol$class %in% input$flmm_res_filter_by_class_vol, ]
                    
                }
                
                table_vol <- table_vol[, colnames(table_vol) == "MSD-site"]
                
                table_vol <- paste(unique(table_vol), collapse = ";")
                
                rclipButton("clip_volume",
                            "Copy the IDs of significative associated MSD-Site to clipboard",
                            table_vol,
                            icon("clipboard"))
            })
            
        }
        
        if ("Stem_dry_mass" %in% input$traits) {
            
            output$manhattan_mass <- renderPlot({
                
                help_manhat_plot(input_df = pvalues_df_chr(),
                                 trait = "Stem_dry_mass",
                                 plot_title = "MSD-Methylated sites associated with volume")
                
            })
            output$bar_mass <- renderPlot({
                
                help_bar_plot(input_df = pvalues_df_chr(),
                              epigenotype_df = epigenotype(),
                              trait = "Stem_dry_mass",
                              plot_title = "Methylation status of the significative marks - Volume")
                
            })
            output$table_mass <- renderReactable({
                
                table_mass <- help_tables(input_df = pvalues_df_chr(),
                                          trait = "Stem_dry_mass",
                                          epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_mass) > 0 ) {
                    
                    table_mass <- table_mass[table_mass$clone %in% input$flmm_res_filter_by_clone_mass, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_mass) > 0 ) {
                    
                    table_mass <- table_mass[table_mass$class %in% input$flmm_res_filter_by_class_mass, ]
                    
                }
                
                reactable(table_mass,
                          defaultColDef = colDef(align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 5,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
                
            })
            output$clip_mass <- renderUI({
                
                table_mass <- help_tables(input_df = pvalues_df_chr(),
                                          trait = "Stem_dry_mass",
                                          epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_mass) > 0 ) {
                    
                    table_mass <- table_mass[table_mass$clone %in% input$flmm_res_filter_by_clone_mass, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_mass) > 0 ) {
                    
                    table_mass <- table_mass[table_mass$class %in% input$flmm_res_filter_by_class_mass, ]
                    
                }
                
                table_mass <- table_mass[, colnames(table_mass) == "MSD-site"]
                
                table_mass <- paste(unique(table_mass), collapse = ";")
                
                rclipButton("clip_steam",
                            "Copy the IDs of significative associated MSD-Site to clipboard",
                            table_mass,
                            icon("clipboard"))
            })
            
        }
        
        if ("Total_height" %in% input$traits) {
            
            output$manhattan_height <- renderPlot({
                
                help_manhat_plot(input_df = pvalues_df_chr(),
                                 trait = "Total_height",
                                 plot_title = "MSD-Methylated sites associated with total height")
                
            })
            output$bar_height <- renderPlot({
                
                help_bar_plot(input_df = pvalues_df_chr(),
                              epigenotype_df = epigenotype(),
                              trait = "Total_height",
                              plot_title = "Methylation status of the significative marks - Total height")
                
            })
            output$table_height <- renderReactable({
                
                table_height <- help_tables(input_df = pvalues_df_chr(),
                                            trait = "Total_height",
                                            epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_height) > 0 ) {
                    
                    table_height <- table_height[table_height$clone %in% input$flmm_res_filter_by_clone_height, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_height) > 0 ) {
                    
                    table_height <- table_height[table_height$class %in% input$flmm_res_filter_by_class_height, ]
                    
                }
                
                reactable(table_height,
                          defaultColDef = colDef(align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 5,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
                
            })
            output$clip_height <- renderUI({
                
                table_to_clip <- help_tables(input_df = pvalues_df_chr(),
                                             trait = "Total_height",
                                             epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_height) > 0 ) {
                    
                    table_to_clip <- table_to_clip[table_to_clip$clone %in% input$flmm_res_filter_by_clone_height, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_height) > 0 ) {
                    
                    table_to_clip <- table_to_clip[table_to_clip$class %in% input$flmm_res_filter_by_class_height, ]
                    
                }
                
                table_to_clip <- table_to_clip[, colnames(table_to_clip) == "MSD-site"]
                
                table_to_clip <- paste(unique(table_to_clip), collapse = ";")
                
                rclipButton("clip_height",
                            "Copy the IDs of significative associated MSD-Site to clipboard",
                            table_to_clip,
                            icon("clipboard"))
            })
        }
        
        if ("DBH" %in% input$traits){
            
            output$manhattan_dbh <- renderPlot({
                
                help_manhat_plot(input_df = pvalues_df_chr(),
                                 trait = "DBH",
                                 plot_title = "MSD-Methylated sites associated with diameter at breast height")
                
            })
            output$bar_dbh <- renderPlot({
                
                help_bar_plot(input_df = pvalues_df_chr(),
                              epigenotype_df = epigenotype(),
                              trait = "DBH",
                              plot_title = "Methylation status of the significative marks - Diameter at Breast Height")
                
            })
            output$table_dbh <- renderReactable({
                
                table_dbh <- help_tables(input_df = pvalues_df_chr(),
                                         trait = "DBH",
                                         epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_dbh) > 0 ) {
                    
                    table_dbh <- table_dbh[table_dbh$clone %in% input$flmm_res_filter_by_clone_dbh, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_dbh) > 0 ) {
                    
                    table_dbh <- table_dbh[table_dbh$class %in% input$flmm_res_filter_by_class_dbh, ]
                    
                }
                
                reactable(table_dbh,
                          defaultColDef = colDef(align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 5,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
            })
            output$clip_dbh <- renderUI({
                
                table_to_clip <- help_tables(input_df = pvalues_df_chr(),
                                             trait = "DBH",
                                             epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_dbh) > 0 ) {
                    
                    table_to_clip <- table_to_clip[table_to_clip$clone %in% input$flmm_res_filter_by_clone_dbh, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_bottom) > 0 ) {
                    
                    table_to_clip <- table_to_clip[table_to_clip$class %in% input$flmm_res_filter_by_class_dbh, ]
                    
                }
                
                table_to_clip <- table_to_clip[, colnames(table_to_clip) == "MSD-site"]
                
                table_to_clip <- paste(unique(table_to_clip), collapse = ";")
                
                rclipButton("clip_dbh",
                            "Copy the IDs of significative associated MSD-Site to clipboard",
                            table_to_clip,
                            icon("clipboard"))
            })
        }
        
        if ("Height_bottom_live_crown" %in% input$traits){
            
            output$manhattan_bottom <- renderPlot({
                
                help_manhat_plot(input_df = pvalues_df_chr(),
                                 trait = "Height_bottom_live_crown",
                                 plot_title = "MSD-Methylated sites associated with Height_bottom_live_crown")
                
            })
            output$bar_bottom <- renderPlot({
                
                help_bar_plot(input_df = pvalues_df_chr(),
                              epigenotype_df = epigenotype(),
                              trait = "Height_bottom_live_crown",
                              plot_title = "Methylation status of the significative marks - Height_bottom_live_crown")
                
            })
            output$table_bottom <- renderReactable({
                
                table_bottom <- help_tables(input_df = pvalues_df_chr(),
                                            trait = "Height_bottom_live_crown",
                                            epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_bottom) > 0 ) {
                    
                    table_bottom <- table_bottom[table_bottom$clone %in% input$flmm_res_filter_by_clone_bottom, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_bottom) > 0 ) {
                    
                    table_bottom <- table_bottom[table_bottom$class %in% input$flmm_res_filter_by_class_bottom, ]
                    
                }
                
                reactable(table_bottom,
                          defaultColDef = colDef(align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 5,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
                
            })
            output$clip_tblc <- renderUI({
                
                table_to_clip <- help_tables(input_df = pvalues_df_chr(),
                                             trait = "Height_bottom_live_crown",
                                             epigenotype_df = epigenotype())
                
                if (length(input$flmm_res_filter_by_clone_bottom) > 0 ) {
                    
                    table_to_clip <- table_to_clip[table_to_clip$clone %in% input$flmm_res_filter_by_clone_bottom, ]
                    
                }
                
                if (length(input$flmm_res_filter_by_class_bottom) > 0 ) {
                    
                    table_to_clip <- table_to_clip[table_to_clip$class %in% input$flmm_res_filter_by_class_bottom, ]
                    
                }
                
                table_to_clip <- table_to_clip[, colnames(table_to_clip) == "MSD-site"]
                
                table_to_clip <- paste(unique(table_to_clip), collapse = ";")
                
                rclipButton("clipbtn",
                            "Copy the IDs of significative associated MSD-Site to clipboard",
                            table_to_clip,
                            icon("clipboard"))
            })
        }
        
    })
    
    ###############################################
    ###                 Tab 2                   ###
    ###############################################
    
    ## Help functions ##
    list_help <- function(clones_samples = annot_clones,
                          data_set = data_set,
                          tissues = tissues) {
        
        if ( length(clones_samples) == 5 ) {
            
            list_clones <- list(c1 = data_set[paste(clones_samples[1], tissues, sep = "_")],
                                c2 = data_set[paste(clones_samples[2], tissues, sep = "_")],
                                c3 = data_set[paste(clones_samples[3], tissues, sep = "_")],
                                c4 = data_set[paste(clones_samples[4], tissues, sep = "_")],
                                c5 = data_set[paste(clones_samples[5], tissues, sep = "_")])
            
        } else if ( length(clones_samples) == 4 ) {
            
            list_clones <- list(c1 = data_set[paste(clones_samples[1], tissues, sep = "_")],
                                c2 = data_set[paste(clones_samples[2], tissues, sep = "_")],
                                c3 = data_set[paste(clones_samples[3], tissues, sep = "_")],
                                c4 = data_set[paste(clones_samples[4], tissues, sep = "_")])
            
        } else if ( length(clones_samples) == 3 ) {
            
            list_clones <- list(c1 = data_set[paste(clones_samples[1], tissues, sep = "_")],
                                c2 = data_set[paste(clones_samples[2], tissues, sep = "_")],
                                c3 = data_set[paste(clones_samples[3], tissues, sep = "_")])
            
        } else if ( length(clones_samples) == 2 ) {
            
            list_clones <- list(c1 = data_set[paste(clones_samples[1], tissues, sep = "_")],
                                c2 = data_set[paste(clones_samples[2], tissues, sep = "_")])
            
        } else if ( length(clones_samples) == 1 ) {
            
            list_clones <- list(c1 = data_set[paste(clones_samples[1], tissues, sep = "_")])
            
        }
        
        list_clones
    }
    
    ## Help function to generate the UpSet plots.
    upset_plot <- function(data_set = data_set,
                           col_fun = col_fun,
                           col_val = col_val,
                           ymax = ymax,
                           sets_order = sets_order){
        
        ht <- UpSet(data_set,
                    set_order = sets_order,
                    pt_size = unit(5, "mm"),
                    lwd = 3,
                    top_annotation = upset_top_annotation(data_set,
                                                          bar_width = 0.5,
                                                          ylim = c(0, ymax),
                                                          annotation_name_rot = 90,
                                                          annotation_name_side = "left",
                                                          axis_param = list(side = "left",
                                                                            at = seq(0, ymax, 500),
                                                                            gp = gpar(fontsize = 12)),
                                                          height = unit(12, "cm")),
                    right_annotation = upset_right_annotation(data_set,
                                                              bar_width = 0.5,
                                                              ylim = c(0, 4000),
                                                              gp = gpar(fill = "grey70"),
                                                              annotation_name_side = "top",
                                                              width = unit(5, "cm"),
                                                              axis_param = list(side = "top",
                                                                                at = seq(0, 4000, 500),
                                                                                gp = gpar(fontsize = 12))),
                    #bottom_annotation = HeatmapAnnotation(
                    # volume = seq(1,26)),
                    left_annotation = rowAnnotation(foo = col_val,
                                                    col = list(foo = col_fun),
                                                    show_annotation_name = FALSE,
                                                    annotation_legend_param = list(
                                                        foo = list(title = expression(paste("Yield (m"^"3", ")")),
                                                                   legend_height = unit(4.5, "cm"),
                                                                   gp = gpar(fontsize = 12)))))
        
        cs <- comb_size(data_set)
        ht <- draw(ht)
        co <- column_order(ht)
        nc <- ncol(data_set)
        decorate_annotation("Intersection\nsize", {
            grid.text(cs[co],
                      x = 1:nc,
                      y = unit(cs[co], "native") + unit(1, "mm"),
                      gp = gpar(fontsize = 9),
                      just = "bottom",
                      default.units = "native")})
        
    }
    
    # loads the annotation function
    anotation_search <- function(query = query,
                                 annotation_file = annotation_file,
                                 list_by = list_by,
                                 output_name = output_name,
                                 query_is_file = "FALSE",
                                 query_is = query_is) {
        
        # Realiza a leitura do conjunto de marcas. Checa se será carregado um arquivo com as ids ou um vetor
        if (query_is_file == "TRUE") {
            
            lista <- vroom(query)
            lista <- as.character(lista$V1)
            
        } else if (query_is_file == "FALSE") {
            
            lista <- as.character(query[!is.na(query)])
            
        } else {
            
            print("Invalid option for 'query_is_file' argument!")
            
        }
        
        # Realiza a leitura do csv com as anotações.
        anot <- vroom(annotation_file, quote = "\"")
        
        # Verifica-se o output será fornecido com uma marca por linha ou uma marca por gene
        if (list_by == "mark") {
            
            # Separa marcas que possuem anotação das que não possuem.
            anotation <- anot[anot$locus %in% lista, ]
            without_anot <- lista[lista %in% anot$locus == "FALSE"]
            
            anotation
            
        } else if (list_by == "gene") {
            
            # Separa marcas que possuem anotação das que não possuem.
            anotation <- anot[anot$locus %in% lista, ]
            without_anot <- lista[lista %in% anot$locus == "FALSE"]
            
            # Junta o nome das marcas que estão relacionadas ao mesmo gene, bem como a strand na qual a marca está inserida e a distância de cada marca ao gene.
            gene_anot_complete <- data.frame()
            for (i in unique(anotation$gene_name)) {
                
                gene <- anotation[anotation$gene_name == i, ]
                locus_gene <- paste(gene$locus, collapse = ";")
                locus_distance <- paste(gene$`distance_to_gene(kb)`, collapse = ";")
                position_to_gene <- paste(gene$position_to_gene, collapse = ";")
                gene_feature <- paste(gene$`gene feature`, collapse = ";")
                
                gene_anot <- gene[1, ]
                gene_anot$locus <- locus_gene
                gene_anot$`distance_to_gene(kb)` <- locus_distance
                gene_anot$position_to_gene <- position_to_gene
                gene_anot$`gene feature` <- gene_feature
                
                gene_anot_complete <- rbind(gene_anot_complete, gene_anot)
            }
            
            gene_anot_complete <- gene_anot_complete %>%
                dplyr::rename("MDS-sites" = "locus",
                              "Chr" = "chromosome",
                              "Dist to gene (Kb)" = "distance_to_gene(kb)",
                              "Gene name" = "gene_name",
                              "Position to gene" = "position_to_gene",
                              "Gene feature" = "gene feature",
                              "Gene description" = "gene_description_merged",
                              "GO ID" = "go_id_merged",
                              "GO description" = "go_term_merged",
                              "KEEG ID" = "kegg_enzyme_id_merged",
                              "KEEG description" = "kegg_enzyme_desc_merged",
                              "Interpro name" = "InterPro.Name",
                              "Interpro GO ID" = "InterPro.GO.ID",
                              "Interpro GO description" = "InterPro.GO.Term")
            
            gene_anot_complete
            
            # Se existirem marcas não anotadas emite um warning e escreve um arquivo com as marcas não anotadas
            
        } else{
            
            print("Invalid 'list_by' argument!")
            print("'list_by' accepts only 'gene' or 'mark'")
            
        }
    }
    
    list_upset_g2 <- reactive({
        
        list_upset <- list_help(input$annot_clones, g2_int, input$annot_tissue)
        
        list_upset <- lapply(list_upset, function(x) x[!is.na(x)])
        
        names(list_upset) <- input$annot_clones
        
        list_upset
        
    })
    
    list_upset_g3 <- reactive({
        
        list_upset <- list_help(input$annot_clones, g3_int, input$annot_tissue)
        
        list_upset <- lapply(list_upset, function(x) x[!is.na(x)])
        
        names(list_upset) <- input$annot_clones
        
        list_upset
        
    })
    
    g2_matrix <- reactive({
        
        matrix <- make_comb_mat(list_upset_g2(),
                                mode = input$upset_mod,
                                remove_empty_comb_set = FALSE)
        
        matrix
    })
    
    g3_matrix <- reactive({
        
        matrix <- make_comb_mat(list_upset_g3(),
                                mode = input$upset_mod,
                                remove_empty_comb_set = FALSE)
        
        matrix
    })
    
    output$upset_g2_leaves <- renderPlot({
        
        prod_g2 <- prod_g2[prod_g2$Clone %in% input$annot_clones, ]
        
        max_value_g2 = max(prod_g2$Volume_mean)
        mean_value_g2 = max_value_g2 / 2
        
        col_fun_g2 = colorRamp2(c(0, mean_value_g2, max_value_g2),
                                c("blue", "yellow", "darkgreen"))
        prod_values_g2 <- prod_g2$Volume_mean
        
        upset_plot(data_set = g2_matrix(),
                   col_fun = col_fun_g2,
                   col_val = prod_values_g2,
                   ymax = 3000,
                   input$annot_clones)
        
    })
    
    output$upset_g3_leaves <- renderPlot({
        
        prod_g3 <- prod_g3[prod_g3$Clone %in% input$annot_clones, ]
        
        max_value_g3 = max(prod_g3$Volume_mean)
        mean_value_g3 = max_value_g3 / 2
        
        col_fun_g3 = colorRamp2(c(0, mean_value_g3, max_value_g3),
                                c("blue", "yellow", "darkgreen"))
        prod_values_g3 <- prod_g3$Volume_mean
        
        upset_plot(data_set = g3_matrix(),
                   col_fun = col_fun_g3,
                   col_val = prod_values_g3,
                   ymax = 3000,
                   input$annot_clones)
        
    })
    
    annotation_g2 <- reactive({
        
        select_group <- ifelse(input$annot_clones %in% input$comps_to_annot, 1, 0)
        select_group <- paste(select_group, collapse = "")
        
        filtered_MSD <- suppressWarnings(extract_comb(g2_matrix(), select_group))
        
        annotated_elements <- anotation_search(query = filtered_MSD,
                                               annotation_file = "data/combined_annotation_blast2GO_biomart.txt",
                                               list_by = "gene",
                                               query_is_file = FALSE,
                                               query_is = "mark")
    })
    
    output$annot_g2 <- renderReactable({
        reactable(annotation_g2(),
                  defaultColDef = colDef(align = "center"),
                  filterable = TRUE,
                  bordered = TRUE,
                  highlight = TRUE,
                  searchable = TRUE,
                  showPageSizeOptions = TRUE,
                  defaultPageSize = 5,
                  pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                  showPagination = TRUE)
    })
    
    output$clip_annot_MG <- renderUI({
        
        data_to_clip <- annotation_g2()
        data_to_clip <- data_to_clip[, colnames(data_to_clip) == input$clip_annot_MG_filter]
        
        data_to_clip <- paste(unique(data_to_clip), collapse = ";")
        
        rclipButton("clipbtn",
                    "Copy the selected ID to clipboard",
                    data_to_clip,
                    icon("clipboard"))
    })
    
    annotation_g3 <- reactive({
        
        select_group <- ifelse(input$annot_clones %in% input$comps_to_annot, 1, 0)
        select_group <- paste(select_group, collapse = "")
        
        filtered_MSD <- suppressWarnings(extract_comb(g3_matrix(), select_group))
        
        annotated_elements <- anotation_search(query = filtered_MSD,
                                               annotation_file = "data/combined_annotation_blast2GO_biomart.txt",
                                               list_by = "gene",
                                               query_is_file = FALSE,
                                               query_is = "mark")
        
        
    })
    
    output$annot_g3 <- renderReactable({
        reactable(annotation_g3(),
                  defaultColDef = colDef(align = "center"),
                  filterable = TRUE,
                  bordered = TRUE,
                  highlight = TRUE,
                  searchable = TRUE,
                  showPageSizeOptions = TRUE,
                  defaultPageSize = 5,
                  pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                  showPagination = TRUE)
    })
    
    output$clip_annot_PR <- renderUI({
        
        data_to_clip <- annotation_g3()
        data_to_clip <- data_to_clip[, colnames(data_to_clip) == input$clip_annot_PR_filter]
        
        data_to_clip <- paste(unique(data_to_clip), collapse = ";")
        
        rclipButton("clipbtn",
                    "Copy the selected ID to clipboard",
                    data_to_clip,
                    icon("clipboard"))
    })
    
    # Starts the enrichment only when the respective botton is used
    observeEvent(input$run_enrich_mg, {
        
        enriched_terms_g2 <- reactive({
            # Reads the list of genes from gff3 file
            query <- annotation_g2()
            query <- query$`Gene name`
            
            # Verifies if all genes or a subset of then should be used as universe in enrichment analysis.
            if (input$enrich_mod == "msd_only") {
                
                universe <- read.table("data/names_of_genes_with_sites_mspI_sequenced.txt", sep = "\t")
                universe <- universe$V1
                
            } else if (input$enrich_mod == "all") {
                
                universe <- read.table("data/all_genes_names_of_e.grandis.tst", sep = "\t")
                universe <- universe$V1
                
            }
            
            # Enrichment parameters
            fdr_cutoff <- input$enr_fdr_cutoff
            
            # Cutoff to remove the redundancy
            sim_cutoff <- input$enr_sim_cutoff  # cutoff to removes similarity
            
            ont <- c("BP", "MF", "CC")
            all_rich_ont_terms <- data.frame()
            
            for (o in ont){
                
                comparison <- clusterProfiler::enrichGO(gene = query,
                                                        universe = universe,
                                                        OrgDb = "org.Egrandis.eg.db",
                                                        keyType = "GID", # Nome da coluna com os genes no OrgDb
                                                        ont = paste(o),
                                                        pAdjustMethod = "fdr",
                                                        qvalueCutoff  = fdr_cutoff,
                                                        readable = F)
                
                # Remove redundância
                comparison_simp <- clusterProfiler::simplify(x = comparison,
                                                             cutoff = sim_cutoff,
                                                             by = "p.adjust",
                                                             select_fun = min)
                
                rich_ont_terms <- as.data.frame(comparison_simp)
                rich_ont_terms$Ontology <- rep(o, nrow(rich_ont_terms))
                
                rich_ont_terms <- dplyr::select(rich_ont_terms,
                                                c("ID",
                                                  "Ontology",
                                                  "Description",
                                                  "GeneRatio",
                                                  "qvalue",
                                                  "Count"
                                                ))
                
                all_rich_ont_terms <- rbind(all_rich_ont_terms,
                                            rich_ont_terms)
                
            }
            
            all_rich_ont_terms
            
            
        })
        
        output$enrich_g2 <- renderReactable({
            reactable(enriched_terms_g2(),
                      defaultColDef = colDef(align = "center"),
                      filterable = TRUE,
                      bordered = TRUE,
                      highlight = TRUE,
                      searchable = TRUE,
                      showPageSizeOptions = TRUE,
                      defaultPageSize = 5,
                      pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                      showPagination = TRUE)
        })
        
        output$clip_enrich_MG <- renderUI({
            
            data_to_clip <- enriched_terms_g2()
            data_to_clip <- data_to_clip[, colnames(data_to_clip) == "ID"]
            
            data_to_clip <- paste(unique(data_to_clip), collapse = ";")
            
            rclipButton("clipbtn",
                        "Copy the IDs of enriched GO terms to clipboard",
                        data_to_clip,
                        icon("clipboard"))
        })
        
        output$enrich_g2_download <- downloadHandler(
            filename = function() {
                paste('enriched_terms_dry_mg_', Sys.Date(), '.csv', sep='')
            },
            content = function(con) {
                write.csv(enriched_terms_g2(), con)
            }
        )
    })
    
    observeEvent(input$run_enrich_pr, {
        enriched_terms_g3 <- reactive({
            
            # Reads the list of genes from gff3 file
            query <- annotation_g3()
            query <- query$`Gene name`
            
            # Verifies if all genes or a subset of then should be used as universe in enrichment analysis.
            if (input$enrich_mod == "msd_only") {
                
                universe <- read.table("data/names_of_genes_with_sites_mspI_sequenced.txt", sep = "\t")
                universe <- universe$V1
                
            } else if (input$enrich_mod == "all") {
                
                universe <- read.table("data/all_genes_names_of_e.grandis.tst", sep = "\t")
                universe <- universe$V1
                
            }
            
            # Enrichment parameters
            fdr_cutoff <- input$enr_fdr_cutoff
            
            # Cutoff to remove the redundancy
            sim_cutoff <- input$enr_sim_cutoff  # cutoff to removes similarity
            
            ont <- c("BP", "MF", "CC")
            all_rich_ont_terms <- data.frame()
            
            for (o in ont){
                
                comparison <- clusterProfiler::enrichGO(gene = query,
                                                        universe = universe,
                                                        OrgDb = "org.Egrandis.eg.db",
                                                        keyType = "GID", # Nome da coluna com os genes no OrgDb
                                                        ont = paste(o),
                                                        pAdjustMethod = "fdr",
                                                        qvalueCutoff  = fdr_cutoff,
                                                        readable = F)
                
                # Remove redundância
                comparison_simp <- clusterProfiler::simplify(x = comparison,
                                                             cutoff = sim_cutoff,
                                                             by = "p.adjust",
                                                             select_fun = min)
                
                rich_ont_terms <- as.data.frame(comparison_simp)
                rich_ont_terms$Ontology <- rep(o, nrow(rich_ont_terms))
                
                rich_ont_terms <- dplyr::select(rich_ont_terms,
                                                c("ID",
                                                  "Ontology",
                                                  "Description",
                                                  "GeneRatio",
                                                  "qvalue",
                                                  "Count"
                                                  #,"geneID"
                                                ))
                
                all_rich_ont_terms <- rbind(all_rich_ont_terms,
                                            rich_ont_terms)
                
            }
            
            all_rich_ont_terms
            
        })
        
        output$enrich_g3 <- renderReactable({
            reactable(enriched_terms_g3(),
                      defaultColDef = colDef(align = "center"),
                      filterable = TRUE,
                      bordered = TRUE,
                      highlight = TRUE,
                      searchable = TRUE,
                      showPageSizeOptions = TRUE,
                      defaultPageSize = 5,
                      pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                      showPagination = TRUE)
        })
        
        output$clip_enrich_PR <- renderUI({
            
            data_to_clip <- enriched_terms_g3()
            data_to_clip <- data_to_clip[, colnames(data_to_clip) == "ID"]
            
            data_to_clip <- paste(unique(data_to_clip), collapse = ";")
            
            rclipButton("clipbtn",
                        "Copy the IDs of enriched GO terms to clipboard",
                        data_to_clip,
                        icon("clipboard"))
        })
        
        output$enrich_g3_download <- downloadHandler(
            filename = function() {
                paste('enriched_terms_wet_pr_', Sys.Date(), '.csv', sep='')
            },
            content = function(con) {
                write.csv(enriched_terms_g3(), con)
            }
        )
    })
    
    output$annot_g2_download <- downloadHandler(
        filename = function() {
            paste('genes_annotation_dry_mg_', Sys.Date(), '.csv', sep='')
        },
        content = function(con) {
            write.csv(annotation_g2(), con)
        }
    )
    
    output$annot_g3_download <- downloadHandler(
        filename = function() {
            paste('genes_annotation_wet_pr_', Sys.Date(), '.csv', sep='')
        },
        content = function(con) {
            write.csv(annotation_g3(), con)
        }
    )
    
    ##------- Tab 2 - Venn and comparison between environmments -------##
    
    output$sites_venn <- renderPlot({
        
        ## G2
        list_upset <- list_help(input$venn_clones, g2_int, input$venn_tissue)
        list_upset <- lapply(list_upset, function(x) x[!is.na(x)])
        names(list_upset) <- input$venn_clones
        MG <- as.data.frame(list_upset[1])
        
        ## G3
        list_upset <- list_help(input$venn_clones, g3_int, input$venn_tissue)
        list_upset <- lapply(list_upset, function(x) x[!is.na(x)])
        names(list_upset) <- input$venn_clones
        PR <- as.data.frame(list_upset[1])
        
        venn_info <- venn.diagram(list(Dry = MG[,1],
                                       Wet = PR[,1]),
                                  fill = c('green', 'yellow'),
                                  alpha = 0.3,
                                  filename = NULL)
        grid.draw(venn_info)
        
    })
    
    venn_df <- reactive({
        
        ## G2
        list_upset <- list_help(input$venn_clones, g2_int, input$venn_tissue)
        list_upset <- lapply(list_upset, function(x) x[!is.na(x)])
        names(list_upset) <- input$venn_clones
        MG <- as.data.frame(list_upset[1])
        
        ## G3
        list_upset <- list_help(input$venn_clones, g3_int, input$venn_tissue)
        list_upset <- lapply(list_upset, function(x) x[!is.na(x)])
        names(list_upset) <- input$venn_clones
        PR <- as.data.frame(list_upset[1])
        
        venn_data <- venn(list(MG = MG, PR = PR), show.plot = F)
        venn_data <- attr(venn_data,"intersections")
        
        if (input$venn_enrich_set == "Dry") {
            
            venn_df <- as.data.frame(venn_data$MG)
            
        } else if (input$venn_enrich_set == "Wet") {
            
            venn_df <- as.data.frame(venn_data$PR)
            
        } else if (input$venn_enrich_set == "Intersection") {
            
            venn_df <- as.data.frame(venn_data$`MG:PR`)
        }
        
        venn_df_annotated_elements <- anotation_search(query = venn_df,
                                                       annotation_file = "data/combined_annotation_blast2GO_biomart.txt",
                                                       list_by = "gene",
                                                       query_is_file = FALSE,
                                                       query_is = "mark")
        venn_df_annotated_elements
        
    })
    
    ## Tab 3 - Comparison between environments ####
    output$venn_table <- renderReactable({
        
        reactable(venn_df(),
                  defaultColDef = colDef(align = "center"),
                  filterable = TRUE,
                  bordered = TRUE,
                  highlight = TRUE,
                  searchable = TRUE,
                  showPageSizeOptions = TRUE,
                  defaultPageSize = 5,
                  pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                  showPagination = TRUE)
        
    })
    
    output$clip_venn <- renderUI({
        
        data_to_clip <- venn_df()
        
        if (input$clip_venn_filter == "MDS-sites") {
            
            data_to_clip <- data_to_clip$`MDS-sites`
            
        } else if (input$clip_venn_filter == "Gene name") {
            
            data_to_clip <- data_to_clip$`Gene name`
            
        } else if (input$clip_venn_filter == "GO ID") {
            
            data_to_clip <- data_to_clip$`GO ID`
            
        }
        
        data_to_clip <- paste(data_to_clip, collapse = ";")
        
        rclipButton("clip_venn",
                    "Copy the selected IDs to clipboard",
                    data_to_clip,
                    icon("clipboard"))
        
    })
    
    observeEvent(input$run_enrich_venn, {
        
        enriched_terms_venn <- reactive({
            
            query <- venn_df()
            query <- query$`Gene name`
            
            # Verifies if all genes or a subset of then should be used as universe in enrichment analysis.
            if (input$venn_enrich_mod == "msd_only") {
                
                universe <- read.table("data/names_of_genes_with_sites_mspI_sequenced.txt", sep = "\t")
                universe <- universe$V1
                
            } else if (input$venn_enrich_mod == "all") {
                
                universe <- read.table("data/all_genes_names_of_e.grandis.tst", sep = "\t")
                universe <- universe$V1
                
            }
            
            # Enrichment parameters
            fdr_cutoff <- input$venn_enr_fdr_cutoff
            
            # Cutoff to remove the redundancy
            sim_cutoff <- input$venn_enr_sim_cutoff  # cutoff to removes similarity
            
            ont <- c("BP", "MF", "CC")
            all_rich_ont_terms <- data.frame()
            
            for (o in ont){
                
                comparison <- clusterProfiler::enrichGO(gene = query,
                                                        universe = universe,
                                                        OrgDb = "org.Egrandis.eg.db",
                                                        keyType = "GID", # Nome da coluna com os genes no OrgDb
                                                        ont = paste(o),
                                                        pAdjustMethod = "fdr",
                                                        qvalueCutoff  = fdr_cutoff,
                                                        readable = F)
                
                # Remove redundância
                comparison_simp <- clusterProfiler::simplify(x = comparison,
                                                             cutoff = sim_cutoff,
                                                             by = "p.adjust",
                                                             select_fun = min)
                
                rich_ont_terms <- as.data.frame(comparison_simp)
                rich_ont_terms$Ontology <- rep(o, nrow(rich_ont_terms))
                
                rich_ont_terms <- dplyr::select(rich_ont_terms,
                                                c("ID",
                                                  "Ontology",
                                                  "Description",
                                                  "GeneRatio",
                                                  "qvalue",
                                                  "Count"
                                                ))
                
                all_rich_ont_terms <- rbind(all_rich_ont_terms,
                                            rich_ont_terms)
                
            }
            
            all_rich_ont_terms
            
        })
        
        output$clip_venn_enrich <- renderUI({
            
            data_to_clip <- enriched_terms_venn()
            data_to_clip <- data_to_clip[, colnames(data_to_clip) == "ID"]
            
            data_to_clip <- paste(unique(data_to_clip), collapse = ";")
            
            rclipButton("clipbtn",
                        "Copy the enriched GO terms to clipboard",
                        data_to_clip,
                        icon("clipboard"))
        })
        
        output$enrich_venn <- renderReactable({
            reactable(enriched_terms_venn(),
                      defaultColDef = colDef(align = "center"),
                      filterable = TRUE,
                      bordered = TRUE,
                      highlight = TRUE,
                      searchable = TRUE,
                      showPageSizeOptions = TRUE,
                      defaultPageSize = 5,
                      pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                      showPagination = TRUE)
        })
        
        output$enrich_venn_download <- downloadHandler(
            filename = function() {
                paste('enriched_terms_in_',
                      input$venn_tissue,
                      "_",
                      input$venn_clones,
                      "_", input$venn_enrich_set,
                      '.csv', sep='')
            },
            content = function(con) {
                enriched_terms_venn <- enriched_terms_venn()
                write.csv(enriched_terms_venn, con)
            }
        )
        
    })
    
    # All annotation
    
    all_genes_annot <- reactive({
        all_genes_annot <- vroom("data/all_genes_annotation_including_MSD_sites.tst")
        
        all_genes_annot <- all_genes_annot[!all_genes_annot$`Gene name` == ".", ] %>%
            dplyr::select(c("Gene name",
                            "Gene chr",
                            "Gene strand",
                            "MSD-site",
                            "MSD-site start",
                            "Distance to the closest gene (kb)",
                            "MSD position in relation to the gene",
                            "Gene description",
                            "GO ID",
                            "GO description",
                            "Enzime code",
                            "Enzime description",
                            "Interpro",
                            "source"))
        
        if (input$annot_filter_by_id_options == "Gene name") {
            
            id_to_filter <- strsplit(as.character(input$annot_filter_by_id), 
                                     split = input$annot_filter_by_id_options_delimiter)
            
            id_to_filter <- id_to_filter[[1]]
            
            all_genes_annot <- all_genes_annot[all_genes_annot$`Gene name` %in% id_to_filter, ]
            
            all_genes_annot
            
        } else if (input$annot_filter_by_id_options == "MDS-Site") {
            
            id_to_filter <- strsplit(as.character(input$annot_filter_by_id), 
                                     split = input$annot_filter_by_id_options_delimiter)
            
            id_to_filter <- id_to_filter[[1]]
            
            all_genes_annot <- anotation_search(query = id_to_filter,
                                                annotation_file = "data/combined_annotation_blast2GO_biomart.txt",
                                                list_by = "gene",
                                                query_is_file = FALSE,
                                                query_is = "mark")
            
            all_genes_annot 
            
        } else if (input$annot_filter_by_id_options == "GO ID") {
            
            id_to_filter <- strsplit(as.character(input$annot_filter_by_id), 
                                     split = input$annot_filter_by_id_options_delimiter)
            
            id_to_filter <- id_to_filter[[1]]
            
            all_genes_annot <- all_genes_annot[all_genes_annot$`GO ID` %in% id_to_filter, ]
            
            all_genes_annot 
            
        }
    })
    
    output$all_genes_annotation <- renderReactable({
        
        reactable(all_genes_annot(),
                  defaultColDef = colDef(align = "center"),
                  filterable = TRUE,
                  bordered = TRUE,
                  highlight = TRUE,
                  searchable = TRUE,
                  showPageSizeOptions = TRUE,
                  defaultPageSize = 5,
                  pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                  showPagination = TRUE,
                  groupBy = "Gene name")
        
    })
    
    observeEvent(input$run_enrich_all_genes, {
        
        enriched_terms_venn_all_genes <- reactive({
            
            query <- all_genes_annot()
            query <- query$`Gene name`
            
            # Verifies if all genes or a subset of then should be used as universe in enrichment analysis.
            if (input$venn_enrich_all_genes_mod == "msd_only") {
                
                universe <- read.table("data/names_of_genes_with_sites_mspI_sequenced.txt", sep = "\t")
                universe <- universe$V1
                
            } else if (input$venn_enrich_all_genes_mod == "all") {
                
                universe <- read.table("data/all_genes_names_of_e.grandis.tst", sep = "\t")
                universe <- universe$V1
                
            }
            
            # Enrichment parameters
            fdr_cutoff <- input$venn_enr_fdr_all_genes_cutoff
            
            # Cutoff to remove the redundancy
            sim_cutoff <- input$venn_enr_sim_all_genes_cutoff  # cutoff to removes similarity
            
            ont <- c("BP", "MF", "CC")
            all_rich_ont_terms <- data.frame()
            
            for (o in ont){
                
                comparison <- clusterProfiler::enrichGO(gene = query,
                                                        universe = universe,
                                                        OrgDb = "org.Egrandis.eg.db",
                                                        keyType = "GID", # Nome da coluna com os genes no OrgDb
                                                        ont = paste(o),
                                                        pAdjustMethod = "fdr",
                                                        qvalueCutoff  = fdr_cutoff,
                                                        readable = F)
                
                # Remove redundância
                comparison_simp <- clusterProfiler::simplify(x = comparison,
                                                             cutoff = sim_cutoff,
                                                             by = "p.adjust",
                                                             select_fun = min)
                
                rich_ont_terms <- as.data.frame(comparison_simp)
                rich_ont_terms$Ontology <- rep(o, nrow(rich_ont_terms))
                
                rich_ont_terms <- dplyr::select(rich_ont_terms,
                                                c("ID",
                                                  "Ontology",
                                                  "Description",
                                                  "GeneRatio",
                                                  "qvalue",
                                                  "Count",
                                                  "geneID"
                                                ))
                
                all_rich_ont_terms <- rbind(all_rich_ont_terms,
                                            rich_ont_terms)
                
            }
            
            all_rich_ont_terms
            
        })
        
        output$clip_venn_enrich_all_genes <- renderUI({
            
            data_to_clip <- enriched_terms_venn_all_genes()
            data_to_clip <- data_to_clip[, colnames(data_to_clip) == "ID"]
            
            data_to_clip <- paste(unique(data_to_clip), collapse = ";")
            
            rclipButton("clip_venn_enrich_all_genes",
                        "Copy the enriched GO terms to clipboard",
                        data_to_clip,
                        icon("clipboard"))
        })
        
        output$enrich_venn_all_genes <- renderReactable({
            reactable(enriched_terms_venn_all_genes(),
                      defaultColDef = colDef(align = "center"),
                      filterable = TRUE,
                      bordered = TRUE,
                      highlight = TRUE,
                      searchable = TRUE,
                      showPageSizeOptions = TRUE,
                      defaultPageSize = 5,
                      pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                      showPagination = TRUE)
        })
        
    })
    
    observe({
        
        if (input$heat_map_tissue == "leaf") {
            
            all_marks_leaves_var <- reactive({
                
                all_marks_leaves <- vroom("data/leaves_epigenotype_techs_transposed.tst",
                                          col_names = T, 
                                          col_types = "fffffffffff")
                
                id_to_filter <- strsplit(input$mds_filter_by_id, 
                                         split = input$mds_filter_by_id_options_delimiter)
                
                id_to_filter <- id_to_filter[[1]]
                
                all_marks_leaves <- all_marks_leaves[all_marks_leaves$MSD_sites %in% id_to_filter, ]
                
                all_marks_leaves <- all_marks_leaves[ , order(names(all_marks_leaves))]
                all_marks_leaves <- all_marks_leaves[, c(9, 1:8, 10:ncol(all_marks_leaves))]
                
            })
            output$all_marks_leaves <- renderReactable({
                
                reactable(all_marks_leaves_var(),
                          defaultColDef = colDef(
                              cell = function(value) {
                                  if (value == 0) paste0("+", value) else value
                              },
                              style = function(value) {
                                  color <- if (value == "Methylated") {
                                      "#63be7b"
                                  } else if (value == "Unmethylated") {
                                      "#ffeb84"
                                  }
                                  list(fontWeight = 600, background = color)
                              },
                              align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 10,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
                
            })
            
            gower.dist <- reactive({
                all_marks_leaves <- all_marks_leaves_var()
                daisy(all_marks_leaves[ ,2:ncol(all_marks_leaves)], metric = c("gower"))
                
            })
            aggl.clust.c <- reactive({
                
                hclust(gower.dist(), method = "complete")
                
            })
            
            output$selec_k_leaves_p1 <- renderPlot({
                
                cstats.table <- function(dist, tree, k) {
                    
                    clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between", "wb.ratio","dunn2","avg.silwidth")
                    
                    clust.size <- c("cluster.size")
                    stats.names <- c()
                    row.clust <- c()
                    
                    if(k >= length(tree$order)){
                        k <- (length(tree$order)-1)
                    }
                    
                    output.stats <- matrix(ncol = k, nrow = length(clust.assess))
                    cluster.sizes <- matrix(ncol = k, nrow = k)
                    
                    for(i in c(1:k)){
                        
                        row.clust[i] <- paste("Cluster-", i, " size")
                    }
                    
                    for(i in c(2:k)){
                        stats.names[i] <- paste("Test", i-1)
                        
                        for(j in seq_along(clust.assess)){
                            output.stats[j, i] <- unlist(cluster.stats(d = dist,
                                                                       clustering = cutree(tree, k = i))[clust.assess])[j]
                            
                        }
                        
                        for(d in 1:k) {
                            
                            cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
                            dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
                            cluster.sizes[d, i]
                            
                        }
                    }
                    
                    output.stats.df <- data.frame(output.stats)
                    
                    cluster.sizes <- data.frame(cluster.sizes)
                    cluster.sizes[is.na(cluster.sizes)] <- 0
                    
                    rows.all <- c(clust.assess, row.clust)
                    output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
                    colnames(output) <- stats.names[2:k]
                    rownames(output) <- rows.all
                    
                    is.num <- sapply(output, is.numeric)
                    output[is.num] <- lapply(output[is.num], round, 2)
                    
                    output
                }
                
                all_marks_df <- all_marks_leaves_var()
                
                # Agglomerative clustering,provides a more ambiguous picture
                ggplot(data = data.frame(t(cstats.table(gower.dist(), aggl.clust.c(), 20))), 
                       aes(x=cluster.number, y=within.cluster.ss)) + 
                    geom_point()+
                    geom_line()+
                    ggtitle("Agglomerative clustering") +
                    labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
                    theme(plot.title = element_text(hjust = 0.5))
                
            })
            output$selec_k_leaves_p2 <- renderPlot({
                
                ## help function
                cstats.table <- function(dist, tree, k) {
                    
                    clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between", "wb.ratio","dunn2","avg.silwidth")
                    
                    clust.size <- c("cluster.size")
                    stats.names <- c()
                    row.clust <- c()
                    
                    if(k >= length(tree$order)){
                        k <- (length(tree$order)-1)
                    }
                    
                    output.stats <- matrix(ncol = k, nrow = length(clust.assess))
                    cluster.sizes <- matrix(ncol = k, nrow = k)
                    
                    for(i in c(1:k)){
                        
                        row.clust[i] <- paste("Cluster-", i, " size")
                    }
                    
                    for(i in c(2:k)){
                        stats.names[i] <- paste("Test", i-1)
                        
                        for(j in seq_along(clust.assess)){
                            output.stats[j, i] <- unlist(cluster.stats(d = dist,
                                                                       clustering = cutree(tree, k = i))[clust.assess])[j]
                            
                        }
                        
                        for(d in 1:k) {
                            
                            cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
                            dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
                            cluster.sizes[d, i]
                            
                        }
                    }
                    
                    output.stats.df <- data.frame(output.stats)
                    
                    cluster.sizes <- data.frame(cluster.sizes)
                    cluster.sizes[is.na(cluster.sizes)] <- 0
                    
                    rows.all <- c(clust.assess, row.clust)
                    
                    output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
                    colnames(output) <- stats.names[2:k]
                    rownames(output) <- rows.all
                    
                    is.num <- sapply(output, is.numeric)
                    output[is.num] <- lapply(output[is.num], round, 2)
                    
                    output
                }
                
                all_marks_df <- all_marks_leaves_var()
                
                # Agglomerative clustering,provides a more ambiguous picture
                # Agglomerative clustering
                ggplot(data = data.frame(t(cstats.table(gower.dist(), aggl.clust.c(), 20))), 
                       aes(x=cluster.number, y=avg.silwidth)) + 
                    geom_point()+
                    geom_line()+
                    ggtitle("Agglomerative clustering") +
                    labs(x = "Num.of clusters", y = "Average silhouette width") +
                    theme(plot.title = element_text(hjust = 0.5))
                
            })
            
            observeEvent(input$run_clustered_heatmap_leaves, {
                
                output$clustered_heatmap_leaves <- renderPlot({
                    
                    all_marks_leaves <- all_marks_leaves_var()
                    
                    ## help function
                    cstats.table <- function(dist, tree, k) {
                        
                        clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between", "wb.ratio","dunn2","avg.silwidth")
                        
                        clust.size <- c("cluster.size")
                        stats.names <- c()
                        row.clust <- c()
                        
                        if(k >= length(tree$order)){
                            k <- (length(tree$order)-1)
                        }
                        
                        output.stats <- matrix(ncol = k, nrow = length(clust.assess))
                        cluster.sizes <- matrix(ncol = k, nrow = k)
                        
                        for(i in c(1:k)){
                            
                            row.clust[i] <- paste("Cluster-", i, " size")
                        }
                        
                        for(i in c(2:k)){
                            stats.names[i] <- paste("Test", i-1)
                            
                            for(j in seq_along(clust.assess)){
                                output.stats[j, i] <- unlist(cluster.stats(d = dist,
                                                                           clustering = cutree(tree, k = i))[clust.assess])[j]
                                
                            }
                            
                            for(d in 1:k) {
                                
                                cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
                                dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
                                cluster.sizes[d, i]
                                
                            }
                        }
                        
                        output.stats.df <- data.frame(output.stats)
                        
                        cluster.sizes <- data.frame(cluster.sizes)
                        cluster.sizes[is.na(cluster.sizes)] <- 0
                        
                        rows.all <- c(clust.assess, row.clust)

                        output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
                        colnames(output) <- stats.names[2:k]
                        rownames(output) <- rows.all
                        
                        is.num <- sapply(output, is.numeric)
                        output[is.num] <- lapply(output[is.num], round, 2)
                        
                        output
                    }
                    
                    ### --- HeatMap --- ###
                    
                    # Set the cluster number
                    clust.num <- cutree(aggl.clust.c(), k = input$n_of_clusters)
                    
                    relation_mark_cluster <- data.frame('MSD_sites' = all_marks_leaves$MSD_sites,
                                                        "cluster_number" = clust.num)
                    
                    dat3 <- melt(all_marks_leaves, id.var = 'MSD_sites')
                    dat3 <- merge(dat3, relation_mark_cluster, by = 'MSD_sites')
                    
                    if(input$filter_cluster_heatmap == "yes") {
                        
                        id_to_filter <- strsplit(as.character(input$clusters_to_keep), 
                                                 split = ",")
                        
                        clusters_to_keep <- id_to_filter[[1]]
                        
                        dat3_filt <- dat3[dat3$cluster_number %in% clusters_to_keep, ]
                        
                    } else if (input$filter_cluster_heatmap == "no") {
                        
                        dat3_filt <- dat3
                        
                    }
                    
                    ggplot(dat3_filt, aes(variable, MSD_sites)) +
                        geom_tile(aes(fill = value), colour = "white") +
                        scale_fill_manual(values=c("darkorange", "blue")) +
                        facet_grid(cluster_number ~ .) +
                        theme_bw()
                    
                    
                })
                
            })
            
        } else if (input$heat_map_tissue == "wood") {
            
            all_marks_leaves_x_var <- reactive({
                
                all_marks_leaves_x <- vroom("data/xylem_epigenotype_techs_transposed.tst",
                                            col_names = T, 
                                            col_types = "fffffffffff") %>%
                    dplyr::rename( "MSD_sites" = "MSD-site")
                
                id_to_filter <- strsplit(input$mds_filter_by_id, 
                                         split = input$mds_filter_by_id_options_delimiter)
                
                id_to_filter <- id_to_filter[[1]]
                
                all_marks_leaves_x <- all_marks_leaves_x[all_marks_leaves_x$MSD_sites %in% id_to_filter, ]
                
                all_marks_leaves_x <- all_marks_leaves_x[ , order(names(all_marks_leaves_x))]
                all_marks_leaves_x <- all_marks_leaves_x[, c(9, 1:8, 10:ncol(all_marks_leaves_x))]
                
            })
            output$all_marks_leaves_x <- renderReactable({
                
                reactable(all_marks_leaves_x_var(),
                          defaultColDef = colDef(
                              cell = function(value) {
                                  if (value == 0) paste0("+", value) else value
                              },
                              style = function(value) {
                                  color <- if (value == "Methylated") {
                                      "#63be7b"
                                  } else if (value == "Unmethylated") {
                                      "#ffeb84"
                                  }
                                  list(fontWeight = 600, background = color)
                              },
                              align = "center"),
                          filterable = TRUE,
                          bordered = TRUE,
                          highlight = TRUE,
                          searchable = TRUE,
                          showPageSizeOptions = TRUE,
                          defaultPageSize = 10,
                          pageSizeOptions = c(5, 10, 25, 50, 100, 200),
                          showPagination = TRUE)
                
            })
            
            gower.dist <- reactive({
                all_marks_leaves_x <- all_marks_leaves_x_var()
                daisy(all_marks_leaves_x[ ,2:ncol(all_marks_leaves_x)], metric = c("gower"))
                
            })
            aggl.clust.c <- reactive({
                
                hclust(gower.dist(), method = "complete")
                
            })
            
            output$selec_k_leaves_p1_x <- renderPlot({
                
                cstats.table <- function(dist, tree, k) {
                    
                    clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between", "wb.ratio","dunn2","avg.silwidth")
                    
                    clust.size <- c("cluster.size")
                    stats.names <- c()
                    row.clust <- c()
                    
                    if(k >= length(tree$order)){
                        k <- (length(tree$order)-1)
                    }
                    
                    output.stats <- matrix(ncol = k, nrow = length(clust.assess))
                    cluster.sizes <- matrix(ncol = k, nrow = k)
                    
                    for(i in c(1:k)){
                        
                        row.clust[i] <- paste("Cluster-", i, " size")
                    }
                    
                    for(i in c(2:k)){
                        stats.names[i] <- paste("Test", i-1)
                        
                        for(j in seq_along(clust.assess)){
                            output.stats[j, i] <- unlist(cluster.stats(d = dist,
                                                                       clustering = cutree(tree, k = i))[clust.assess])[j]
                            
                        }
                        
                        for(d in 1:k) {
                            
                            cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
                            dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
                            cluster.sizes[d, i]
                            
                        }
                    }
                    
                    output.stats.df <- data.frame(output.stats)
                    
                    cluster.sizes <- data.frame(cluster.sizes)
                    cluster.sizes[is.na(cluster.sizes)] <- 0
                    
                    rows.all <- c(clust.assess, row.clust)

                    output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
                    colnames(output) <- stats.names[2:k]
                    rownames(output) <- rows.all
                    
                    is.num <- sapply(output, is.numeric)
                    output[is.num] <- lapply(output[is.num], round, 2)
                    
                    output
                }
                
                all_marks_df <- all_marks_leaves_x_var()
                
                # Agglomerative clustering,provides a more ambiguous picture
                ggplot(data = data.frame(t(cstats.table(gower.dist(), aggl.clust.c(), 20))), 
                       aes(x=cluster.number, y=within.cluster.ss)) + 
                    geom_point()+
                    geom_line()+
                    ggtitle("Agglomerative clustering") +
                    labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
                    theme(plot.title = element_text(hjust = 0.5))
                
            })
            output$selec_k_leaves_p2_x <- renderPlot({
                
                ## help function
                cstats.table <- function(dist, tree, k) {
                    
                    clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between", "wb.ratio","dunn2","avg.silwidth")
                    
                    clust.size <- c("cluster.size")
                    stats.names <- c()
                    row.clust <- c()
                    
                    if(k >= length(tree$order)){
                        k <- (length(tree$order)-1)
                    }
                    
                    output.stats <- matrix(ncol = k, nrow = length(clust.assess))
                    cluster.sizes <- matrix(ncol = k, nrow = k)
                    
                    for(i in c(1:k)){
                        
                        row.clust[i] <- paste("Cluster-", i, " size")
                    }
                    
                    for(i in c(2:k)){
                        stats.names[i] <- paste("Test", i-1)
                        
                        for(j in seq_along(clust.assess)){
                            output.stats[j, i] <- unlist(cluster.stats(d = dist,
                                                                       clustering = cutree(tree, k = i))[clust.assess])[j]
                            
                        }
                        
                        for(d in 1:k) {
                            
                            cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
                            dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
                            cluster.sizes[d, i]
                            
                        }
                    }
                    
                    output.stats.df <- data.frame(output.stats)
                    
                    cluster.sizes <- data.frame(cluster.sizes)
                    cluster.sizes[is.na(cluster.sizes)] <- 0
                    
                    rows.all <- c(clust.assess, row.clust)

                    output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
                    colnames(output) <- stats.names[2:k]
                    rownames(output) <- rows.all
                    
                    is.num <- sapply(output, is.numeric)
                    output[is.num] <- lapply(output[is.num], round, 2)
                    
                    output
                }
                
                all_marks_df <- all_marks_leaves_x_var()
                
                # Agglomerative clustering,provides a more ambiguous picture
                # Agglomerative clustering
                ggplot(data = data.frame(t(cstats.table(gower.dist(), aggl.clust.c(), 20))), 
                       aes(x=cluster.number, y=avg.silwidth)) + 
                    geom_point()+
                    geom_line()+
                    ggtitle("Agglomerative clustering") +
                    labs(x = "Num.of clusters", y = "Average silhouette width") +
                    theme(plot.title = element_text(hjust = 0.5))
                
            })
            
            observeEvent(input$run_clustered_heatmap_leaves_x, {
                
                output$clustered_heatmap_leaves_x <- renderPlot({
                    
                    all_marks_leaves_x <- all_marks_leaves_x_var()
                    
                    ## help function
                    cstats.table <- function(dist, tree, k) {
                        
                        clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between", "wb.ratio","dunn2","avg.silwidth")
                        
                        clust.size <- c("cluster.size")
                        stats.names <- c()
                        row.clust <- c()
                        
                        if(k >= length(tree$order)){
                            k <- (length(tree$order)-1)
                        }
                        
                        output.stats <- matrix(ncol = k, nrow = length(clust.assess))
                        cluster.sizes <- matrix(ncol = k, nrow = k)
                        
                        for(i in c(1:k)){
                            
                            row.clust[i] <- paste("Cluster-", i, " size")
                        }
                        
                        for(i in c(2:k)){
                            stats.names[i] <- paste("Test", i-1)
                            
                            for(j in seq_along(clust.assess)){
                                output.stats[j, i] <- unlist(cluster.stats(d = dist,
                                                                           clustering = cutree(tree, k = i))[clust.assess])[j]
                                
                            }
                            
                            for(d in 1:k) {
                                
                                cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
                                dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
                                cluster.sizes[d, i]
                                
                            }
                        }
                        
                        output.stats.df <- data.frame(output.stats)
                        
                        cluster.sizes <- data.frame(cluster.sizes)
                        cluster.sizes[is.na(cluster.sizes)] <- 0
                        
                        rows.all <- c(clust.assess, row.clust)

                        output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
                        colnames(output) <- stats.names[2:k]
                        rownames(output) <- rows.all
                        
                        is.num <- sapply(output, is.numeric)
                        output[is.num] <- lapply(output[is.num], round, 2)
                        
                        output
                    }
                    
                    ### --- HeatMap --- ###
                    
                    # Set the cluster number
                    clust.num <- cutree(aggl.clust.c(), k = input$n_of_clusters_x)
                    
                    relation_mark_cluster <- data.frame('MSD_sites' = all_marks_leaves_x$MSD_sites,
                                                        "cluster_number" = clust.num)
                    
                    dat3 <- melt(all_marks_leaves_x, id.var = 'MSD_sites')
                    dat3 <- merge(dat3, relation_mark_cluster, by = 'MSD_sites')
                    
                    if(input$filter_cluster_heatmap_x == "yes") {
                        
                        id_to_filter <- strsplit(as.character(input$clusters_to_keep_x), 
                                                 split = ",")
                        
                        clusters_to_keep_x <- id_to_filter[[1]]
                        
                        dat3_filt <- dat3[dat3$cluster_number %in% clusters_to_keep_x, ]
                        
                    } else if (input$filter_cluster_heatmap_x == "no") {
                        
                        dat3_filt <- dat3
                        
                    }
                    
                    ggplot(dat3_filt, aes(variable, MSD_sites)) +
                        geom_tile(aes(fill = value), colour = "white") +
                        scale_fill_manual(values=c("darkorange", "blue")) +
                        facet_grid(cluster_number ~ .) +
                        theme_bw()
                    
                    
                })
                
            })
            
        }
        
    })
})

