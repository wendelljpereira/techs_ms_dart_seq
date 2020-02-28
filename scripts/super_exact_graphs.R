# carrega o pacote para o teste
save.image("SuperExactTest.rda")

require("SuperExactTest")
require("gdata")

# Reads the data files with Differential Methylated marks
g2_int <- read.table(snakemake@input[[1]],
                     header = T,
                     sep = "\t",
                     colClasses = "character")

g3_int <- read.table(snakemake@input[[2]],
                     header = T,
                     sep = "\t",
                     colClasses = "character")

# Define the population size.
total <- snakemake@params[["population_size_marks"]]
tissue <- snakemake@params[["tissues"]]
place <- snakemake@params[["places"]]
sufix_exclusive_in_each_place <- snakemake@params[["sufix_intersect_in_each_place"]]
sufix_same_clone <- snakemake@params[["sufix_same_clone"]]
sufix_same_place <- snakemake@params[["sufix_same_place"]]
groups_graphs_3 <- snakemake@params[["groups_graphs_3"]]

# Builds the groups for comparision between the samples of the same clone.
clones <- snakemake@params[["genotypes"]]
for(c in clones){
    
    clone <- list(Leaf_MG = g2_int[paste(c, "leaf", sep = "_")],
                  Leaf_PR = g3_int[paste(c, "leaf", sep = "_")],
                  Wood_MG = g2_int[paste(c, "wood", sep = "_")],
                  Wood_PR = g3_int[paste(c, "wood", sep = "_")])
    
    # Removes NAs
    clone <- lapply(clone, function(x) x[!is.na(x)])
    
    # Set names to the graph
    names(clone) <- c(paste(c, "Leaf MG"),
                      paste(c, "Leaf PR"),
                      paste(c, "Wood MG"),
                      paste(c, "Wood PR"))
    
    # Executes the test
    clone_res <- supertest(clone, n = total)
    
    # Make and save the plot
    svg(filename = paste("images/SuperExactTest/", c, sufix_same_clone, ".svg", sep = ""), width = 12, height = 8)
    
    plot(clone_res,
         sort.by = "degree",
         degree = c(1:5),
         legend.col = 1,
         track.area.range = 0.3,
         bar.area.range = 0.10)
    
    dev.off()
    
    # Save the summary file
    write.csv(summary(clone_res)$Table,
              file = paste("super_exact_out/summary_", c, sufix_same_clone, ".txt", sep = ""),
              row.names = FALSE)
}

# Builds the groups for comparision between clones in the same place
for(l in 2 : (length(place) + 1)){
    
    # Select the order for the graph by the productivity in each place. 
    if(l == 2){
        
        clones <- snakemake@params[["genotypes_MG_prod_order"]]
        
    }else if(l == 3){
        
        clones <- snakemake@params[["genotypes_PR_prod_order"]]
        
    }
    
    #Builds the groups for comparision
    for(t in tissue){
        
        loc <- get(paste("g", l, "_int", sep = ""))
        
        clone <- list(c1 = loc[paste(clones[1], t, sep = "_")],
                      c2 = loc[paste(clones[2], t, sep = "_")],
                      c3 = loc[paste(clones[3], t, sep = "_")],
                      c4 = loc[paste(clones[4], t, sep = "_")],
                      c5 = loc[paste(clones[5], t, sep = "_")])
        
        # Removes NAs
        clone <- lapply(clone, function(x) x[!is.na(x)])
        
        # Set names to the graph
        names(clone) <- c(paste(clones[1], t, place[l-1], sep = "_"),
                          paste(clones[2], t, place[l-1], sep = "_"),
                          paste(clones[3], t, place[l-1], sep = "_"),
                          paste(clones[4], t, place[l-1], sep = "_"),
                          paste(clones[5], t, place[l-1], sep = "_"))
        
        # Executes the test
        clone_res <- supertest(clone, n = total)
        
        # Make and save the plot
        svg(filename = paste0("images/SuperExactTest/",
                              paste(t, sufix_same_place, place[l - 1], sep = "_"), ".svg"),
            width = 12,
            height = 8)
        
        plot(clone_res,
             sort.by = "degree",
             degree = c(1:5),
             legend.col = 1,
             track.area.range = 0.3,
             bar.area.range = 0.10)
        
        dev.off()
        
        # Save the summary file
        write.csv(summary(clone_res)$Table,
                  file = paste0(paste("super_exact_out/summary", t, sufix_same_place, place[l-1], sep = "_"), ".txt"),
                  row.names = FALSE)
    }
}

# Builds the groups for comparision betwenn marks especificaly DM in each place
for(t in tissue){
    for(g in groups_graphs_3){
        
        if(t == "leaf"){
            
            if(g == "unique_MG"){
                
                inter_sufix <- snakemake@params[["sufix_input_files_leaf"]][1]
                clones <- snakemake@params[["genotypes_MG_prod_order"]]
                
            } else if(g == "intersect_MG_PR"){
                
                inter_sufix <- snakemake@params[["sufix_input_files_leaf"]][2]
                clones <- snakemake@params[["genotypes"]]
                
            } else if(g == "unique_PR"){
                
                inter_sufix <- snakemake@params[["sufix_input_files_leaf"]][3]
                clones <- snakemake@params[["genotypes_PR_prod_order"]]
                
            }
            
        }else if(t == "wood"){
            
            if(g == "unique_MG"){
                
                inter_sufix <- snakemake@params[["sufix_input_files_wood"]][1]
                clones <- snakemake@params[["genotypes_MG_prod_order"]]
                
            }else if(g == "intersect_MG_PR"){
                
                inter_sufix <- snakemake@params[["sufix_input_files_wood"]][2]
                clones <- snakemake@params[["genotypes"]]
                
            }else if(g == "unique_PR"){
                
                inter_sufix <- snakemake@params[["sufix_input_files_wood"]][3]
                clones <- snakemake@params[["genotypes_PR_prod_order"]]
                
            }
        }
        
        # Extract marks of each group and clone
        for(c in 1:length(clones)){
            
            if(c == 1){
                
                marks <- read.table(file = paste("intersections_to_bar_plot/",
                                                 clones[c], inter_sufix, sep = ""),
                                    colClasses = "character")
                
                names(marks) <- paste(clones[c], g, sep="_")
                marks_table <- marks
                
            }else if(c > 1){
                
                marks <- read.table(file = paste("intersections_to_bar_plot/",
                                                 clones[c], inter_sufix, sep = ""),
                                    colClasses = "character")
                names(marks) <- paste(clones[c], g, sep="_")
                marks_table <- cbindX(marks_table, marks)
                
            }
        }  
        
        clone <- list(c1 = marks_table[paste(clones[1], g, sep = "_")],
                      c2 = marks_table[paste(clones[2], g, sep = "_")],
                      c3 = marks_table[paste(clones[3], g, sep = "_")],
                      c4 = marks_table[paste(clones[4], g, sep = "_")],
                      c5 = marks_table[paste(clones[5], g, sep = "_")])
        
        # Removes NAs
        clone <- lapply(clone, function(x) x[!is.na(x)])
        
        # Set names to the graph
        names(clone) <- c(paste(clones[1], g, sep = "_"),
                          paste(clones[2], g, sep = "_"),
                          paste(clones[3], g, sep = "_"),
                          paste(clones[4], g, sep = "_"),
                          paste(clones[5], g, sep = "_"))
        
        # Executes the test
        clone_res <- supertest(clone, n = total)
        
        # Make and save the plot
        svg(filename = paste0("images/SuperExactTest/", 
                              paste(sufix_exclusive_in_each_place, t, g, sep = "_") , ".svg"), width = 12, height= 8)
        plot(clone_res,
             sort.by = "degree",
             degree = c(1:5),
             legend.col = 1,
             track.area.range = 0.3,
             bar.area.range = 0.10)
        
        dev.off()
        
        # Save the summary file
        write.csv(summary(clone_res)$Table,
                  file = paste0(paste("super_exact_out/summary",
                                      sufix_exclusive_in_each_place, t, g, sep = "_"), ".txt"), 
                  row.names=FALSE)
    }
}
