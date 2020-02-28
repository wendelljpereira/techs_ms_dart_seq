################
### Box plot between each environment ###
################
require(tidyverse)
require(ggstatsplot)
require(ggpubr)
require(ggthemes)

## Loading phenotype data ##

phenot_tree_level <- read.table(paste(snakemake@input[[1]]),
                                header = T,
                                sep = ",") %>%
    filter(date_for_measurement == "10/1/15") %>%
    filter(Throughfall_exclusion == "Não") %>%
    mutate(Clone = gsub(pattern = " ", "", Clone)) %>%
    mutate(location = ifelse(Site == 30, "Dry", "Wet")) %>%
    mutate(clone_location = paste(Clone, location, sep = " - "))

## plot - Volume
vol <- ggstatsplot::grouped_ggwithinstats(
    data = dplyr::filter(
        .data = phenot_tree_level,
        Clone %in% c("A1", "B2", "D4", "E5", "Q8")
    ),
    x = clone_location,
    y = Volume,
    k = 2,
    grouping.var = Clone, # grouping variable
    pairwise.comparisons = TRUE, # display significant pairwise comparisons
    pairwise.annotation = "p.value",
    pairwise.display = "all",
    p.adjust.method = "fdr",
    results.subtitle = T,
    sample.size.label = FALSE,
    mean.ci = FALSE,
    conf.level = 0.95, # changing confidence level to 99%
    ylab = "Volume (m3)",
    xlab = "",
    title.prefix = "Clone",
    mean.color = "darkblue",
    messages = T,
    bf.message = F,
    package = "RColorBrewer",
    palette = "Paired",
    mean.label.size = 5,
    mean.size = 3,
    notch = T,
    nrow = 1,
    title.text = expression("Comparison of clone productivity in the wet and dry locations - Volume")) +
    ggplot2::theme(axis.text.x = element_text(color = "grey20",
                                              size = 11),
                   axis.text.y = element_text(color = "grey20",
                                              size = 14),
                   axis.title.y = element_text(size = 18),
                   text = element_text(size = 14))

##height
height <- ggstatsplot::grouped_ggwithinstats(
    data = dplyr::filter(
        .data = phenot_tree_level,
        Clone %in% c("A1", "B2", "D4", "E5", "Q8")
    ),
    x = clone_location,
    y = Total_height,
    k = 2,
    grouping.var = Clone, # grouping variable
    pairwise.comparisons = TRUE, # display significant pairwise comparisons
    pairwise.annotation = "p.value",
    pairwise.display = "all",
    p.adjust.method = "fdr",
    results.subtitle = T,
    sample.size.label = FALSE,
    mean.ci = FALSE,
    conf.level = 0.95, # changing confidence level to 99%
    ylab = "Height (m)",
    xlab = "",
    title.prefix = "Clone",
    mean.color = "darkblue",
    messages = T,
    bf.message = F,
    package = "RColorBrewer",
    palette = "Paired",
    mean.label.size = 5,
    mean.size = 3,
    notch = T,
    nrow = 1,
    title.text = expression("Comparison of clone productivity in the wet and dry locations - Height")) +
    ggplot2::theme(axis.text.x = element_text(color = "grey20",
                                              size = 11),
                   axis.text.y = element_text(color = "grey20",
                                              size = 14),
                   axis.title.y = element_text(size = 18),
                   text = element_text(size = 14))

## DBH
dbh <- ggstatsplot::grouped_ggwithinstats(
    data = dplyr::filter(
        .data = phenot_tree_level,
        Clone %in% c("A1", "B2", "D4", "E5", "Q8")
    ),
    x = clone_location,
    y = DBH,
    k = 2,
    grouping.var = Clone, # grouping variable
    pairwise.comparisons = TRUE, # display significant pairwise comparisons
    pairwise.annotation = "p.value",
    pairwise.display = "all",
    p.adjust.method = "fdr",
    results.subtitle = T,
    sample.size.label = FALSE,
    mean.ci = FALSE,
    conf.level = 0.95, # changing confidence level to 99%
    ylab = "DBH (cm)",
    xlab = "",
    title.prefix = "Clone",
    mean.color = "darkblue",
    messages = T,
    bf.message = F,
    package = "RColorBrewer",
    palette = "Paired",
    mean.label.size = 5,
    mean.size = 3,
    notch = T,
    nrow = 1,
    title.text = expression("Comparison of clone productivity in the wet and dry locations - DBH")) +
    ggplot2::theme(axis.text.x = element_text(color = "grey20",
                                              size = 20),
                   axis.text.y = element_text(color = "grey20",
                                              size = 13),
                   axis.title.y = element_text(size = 18),
                   text = element_text(size = 16))

svg(filename = paste(snakemake@output[[1]]),
    width = 22,
    height = 14,
    pointsize = 12)
ggstatsplot::combine_plots(
    vol,height,dbh,
    nrow = 3,
    labels = c(""),
    title.text = "",
    #    caption.text = "Source: Gapminder Foundation",
    title.size = 14,
    caption.size = 12)
dev.off()

################
### Box plot within each environment ###
################

phenot_tree_level <- read.table(paste(snakemake@input[[1]]),
                                header = T,
                                sep = ",") %>%
    filter(date_for_measurement == "10/1/15") %>%
    filter(Throughfall_exclusion == "Não") %>%
    mutate(Clone = gsub(pattern = " ", "", Clone)) %>%
    mutate(location = ifelse(Site == 30, "Dry", "Wet")) %>%
    mutate(clone_location = paste(Clone, location, sep = " - "))


## plot - Volume

vol <- ggstatsplot::grouped_ggwithinstats(
    data = dplyr::filter(
        .data = phenot_tree_level,
        location %in% c("Dry", "Wet")
    ),
    x = Clone,
    y = Volume,
    k = 2,
    grouping.var = location, # grouping variable
    pairwise.comparisons = TRUE, # display significant pairwise comparisons
    pairwise.annotation = "asterisk",
    pairwise.display = "all",
    p.adjust.method = "fdr",
    results.subtitle = T,
    sample.size.label = FALSE,
    mean.ci = FALSE,
    conf.level = 0.95, # changing confidence level to 99%
    ylab = "Volume (m3)",
    xlab = "",
    title.prefix = "Clone",
    mean.color = "darkblue",
    messages = T,
    bf.message = F,
    package = "RColorBrewer",
    palette = "Paired",
    mean.label.size = 5,
    mean.size = 3,
    notch = T,
    nrow = 1,
    title.text = "Comparison of clones productivity within the dry (left) and wet (right) locations - Volume") +
    ggplot2::theme(axis.text.x = element_text(color = "grey20",
                                              size = 11),
                   axis.text.y = element_text(color = "grey20",
                                              size = 14),
                   axis.title.y = element_text(size = 18),
                   text = element_text(size = 14))

##height
height <- ggstatsplot::grouped_ggwithinstats(
    data = dplyr::filter(
        .data = phenot_tree_level,
        location %in% c("Dry", "Wet")
    ),
    x = Clone,
    y = Total_height,
    k = 2,
    grouping.var = location, # grouping variable
    pairwise.comparisons = TRUE, # display significant pairwise comparisons
    pairwise.annotation = "asterisk",
    pairwise.display = "all",
    p.adjust.method = "fdr",
    results.subtitle = T,
    sample.size.label = FALSE,
    mean.ci = FALSE,
    conf.level = 0.95, # changing confidence level to 99%
    ylab = "Height (m)",
    xlab = "",
    title.prefix = "Clone",
    mean.color = "darkblue",
    messages = T,
    bf.message = F,
    package = "RColorBrewer",
    palette = "Paired",
    mean.label.size = 5,
    mean.size = 3,
    notch = T,
    nrow = 1,
    title.text = "Comparison of clones productivity within the dry (left) and wet (right) locations - Height") +
    ggplot2::theme(axis.text.x = element_text(color = "grey20",
                                              size = 11),
                   axis.text.y = element_text(color = "grey20",
                                              size = 14),
                   axis.title.y = element_text(size = 18),
                   text = element_text(size = 14))

## DBH
dbh <- ggstatsplot::grouped_ggwithinstats(
    data = dplyr::filter(
        .data = phenot_tree_level,
        location %in% c("Dry", "Wet")
    ),
    x = Clone,
    y = DBH,
    k = 2,
    grouping.var = location, # grouping variable
    pairwise.comparisons = TRUE, # display significant pairwise comparisons
    pairwise.annotation = "asterisk",
    pairwise.display = "all",
    p.adjust.method = "fdr",
    results.subtitle = T,
    sample.size.label = FALSE,
    mean.ci = FALSE,
    conf.level = 0.95, # changing confidence level to 99%
    ylab = "DBH (cm)",
    xlab = "",
    title.prefix = "Clone",
    mean.color = "darkblue",
    messages = T,
    bf.message = F,
    package = "RColorBrewer",
    palette = "Paired",
    mean.label.size = 5,
    mean.size = 3,
    notch = T,
    nrow = 1,
    title.text = "Comparison of clones productivity in the wet and dry locations - DBH") +
    ggplot2::theme(axis.text.x = element_text(color = "grey20",
                                              size = 11),
                   axis.text.y = element_text(color = "grey20",
                                              size = 14),
                   axis.title.y = element_text(size = 18),
                   text = element_text(size = 14))

svg(filename = paste(snakemake@output[[2]]),
    width = 14,
    height = 24,
    pointsize = 12)
ggstatsplot::combine_plots(
    vol,height,dbh,
    nrow = 3,
    labels = c("A)", "B)", "C)"),
    title.text = "",
    #    caption.text = "Source: Gapminder Foundation",
    title.size = 14,
    caption.size = 12)
dev.off()

