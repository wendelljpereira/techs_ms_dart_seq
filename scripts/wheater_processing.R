require(tidyverse)

wheater_22KLT <- read_csv("TECHS_Weather.csv") %>%
    filter(TECHS == "22KLT") %>%
    mutate(month_year = paste(Month, Year, sep = "_")) %>%
    group_by(month_year) %>%
    summarise_if(funs(mean), .predicate = is.numeric) %>%
    mutate(localite = "22KLT")

wheater_30VMT <- read_csv("TECHS_Weather.csv") %>%
    filter(TECHS == "30VMT") %>%
    mutate(month_year = paste(Month, Year, sep = "_")) %>%
    group_by(month_year) %>%
    summarise_if(funs(mean), .predicate = is.numeric) %>%
    mutate(localite = "30VMT")

time <- unique(wheater_22KLT$month_year)

## 22KLT
new_wheater_22KLT <- data.frame()
for(i in 1:length(time)){
    
    sub_22 <- wheater_22KLT[wheater_22KLT$month_year == time[i], -1]
    colnames(sub_22) <- paste(colnames(sub_22), time[i], sep = "|")
    
    if (i == 1){
        
        new_wheater_22KLT <- sub_22
        
    } else {
        
        new_wheater_22KLT <- cbind(new_wheater_22KLT, sub_22)
        
    }
    
}

## 30 VMT
new_wheater_30VMT <- data.frame()
for(i in 1:length(time)){
    
    sub_22 <- wheater_30VMT[wheater_30VMT$month_year == time[i], -1]
    colnames(sub_22) <- paste(colnames(sub_22), time[i], sep = "|")
    
    if (i == 1){
        
        new_wheater_30VMT <- sub_22
        
    } else {
        
        new_wheater_30VMT <- cbind(new_wheater_30VMT, sub_22)
        
    }
    
}

df_22KLT <- new_wheater_22KLT[rep(seq_len(nrow(new_wheater_22KLT)), each = 5), ]
df_30VMT <- new_wheater_30VMT[rep(seq_len(nrow(new_wheater_30VMT)), each = 5), ]

clones <- c("AEC144", "CNB10", "GG100", "FIB6075", "VER361")

Clones_22KLT <- as.data.frame(paste(clones, "22", sep = "_"))
df_22KLT <- cbind(Clones_22KLT, df_22KLT)
colnames(df_22KLT)[1] <- "clones"

Clones_30VMT <- as.data.frame(paste(clones, "30", sep = "_"))
df_30VMT <- cbind(Clones_30VMT, df_30VMT)
colnames(df_30VMT)[1] <- "clones"
    
wheater_month_ave <- rbind(df_22KLT, df_30VMT)

write.table(wheater_month_ave, 
            "wheater_month_ave.csv",
            col.names = T,
            row.names = F,
            quote = F,
            sep = ",")
