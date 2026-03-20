## Clear Global Environment
rm(list = ls())

## Directories
dir_1718_H3N2 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H3N2/realigned",sep="")
dir_1819_H3N2 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H3N2/realigned",sep="")
dir_1718_H1N1 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H1N1/realigned",sep="")
dir_1819_H1N1 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H1N1/realigned",sep="")
dir_s         <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/figures/supp",sep="")

#### Session prep ####
## Install packages and load libraries as required
if(!require(tidyverse)){
  install.packages("tidyverse",dependencies = T)
  library(tidyverse)
}
if(!require(ggplot2)){
  install.packages("ggplot2",dependencies = T)
  library(ggplot2)
}
if(!require(cowplot)){
  install.packages("cowplot",dependencies = T)
  library(cowplot)
}

#### Plot standards ####
## factors
factor_by_gene <- c("PB2", "PB1", "PA", "HA", "HA_antigenic", "HA_nonantigenic", "NP", "Neuraminidase", "M1", "M2", "NS1", "NEP", "PA-X", "PB1-F2")
factor_by_segment <- c("PB2", "PB1", "PA", "HA", "NP", "Neuraminidase", "M1", "M2", "NS1", "NEP")
factor_by_gene_cleaned <- c("PB2", "PB1", "PA", "HA", "Anti. HA", "Nonanti. HA", "NP", "NA", "M1", "M2", "NS1", "NEP", "PA-X", "PB1-F2"   )

## theme
axis_formatting <- theme(axis.text.x = element_text(size = 6),
                         axis.text.y = element_text(size = 6),
                         axis.title.x = element_text(size = 6, margin = margin(t = 4)),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))

legend_formatting <- theme(legend.text = element_text(size = 6),
                           legend.key.height= unit(0.5, 'cm'),
                           legend.key.width= unit(0.5, 'cm'))

background_formatting <- theme(panel.border = element_rect(color = "grey", fill = NA, size = .5),
                               panel.grid = element_blank(),
                               strip.background = element_blank(),
                               panel.background = element_blank(),
                               legend.background = element_blank())

## plot_grid
Size_adjust = 12
LR_adjust = -0.5 # less = right
UD_adjust = 1.1 # less = down 

## palettes
# mutation types
palette_muts <- c("Nonsynonymous" = "#FF7F20",
                  "Synonymous" = "#4F7899",
                  "Nonsense" = "#7AD9C2")
# mutation types without stop
palette_muts_NS_S <- c("Nonsynonymous" = "#FF7F20",
                       "Synonymous" = "#4F7899")
# flu subtypes and B
palette_subtypes <- c("H3N2" = "#035953",
                      "H1N1" = "#017F6C")
palette_subtypes_year <- c("17-18_H3N2" = "#137C74",
                           "18-19_H3N2" = "#035953",
                           "17-18_H1N1" = "#017F6C",
                           "18-19_H1N1" = "#017F6C")
palette_subtypes_H3N2 <- c("17-18_H3N2" = "#137C74",
                           "18-19_H3N2" = "#035953")

# vax and unvax
palette_vax <- c("Vaccinated" = "Light Grey",
                 "Unvaccinated" = "Dark grey")
# pairs
palette_pairs <- c("Household" = "Red", 
                   "Random" = "Grey")


#### Functions ####
## data summary
ds <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm=TRUE),
      median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sum(sd(x[[col]]) / sqrt(length(x[[1]]))))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = paste("mean", varname, sep = "_")))
  data_sum <- rename(data_sum, c("median" = paste("median", varname, sep = "_")))
  return(data_sum)
}
# examples:
# x <- as.data.frame(ds(df, 
#                       varname = "varname_of_interest", 
#                       groupnames = c("group_of_interest")))
# x <- as.data.frame(ds(df, 
#                       varname = "varname_of_interest", 
#                       groupnames = c("group1", "group2", "...")))

## not in
'%!in%' <- function(x,y)!('%in%'(x,y))
# useful for removing column values of y from df x
# e.g., x <- filter(x,y %!in% c("unwantedValue1", "unwantedValue2"))


#### Import 1718_H3N2 ####
## create dir
setwd(dir_1718_H3N2); getwd(); dir()
percov_1718_H3N2 <- read.table("percent_coverage.tsv", sep = " ", header = FALSE)
colnames(percov_1718_H3N2) <- c("WSLH_ID", "Per_Coverage")
percov_1718_H3N2 <- percov_1718_H3N2[percov_1718_H3N2$WSLH_ID!="FileGene",]
bstats_1718_H3N2 <- as.data.frame(read_tsv("bam_stats.tsv"))
to_remove <- as.data.frame(read_tsv("to_remove.tsv"))

## create list for percov tsvs
setwd(paste(dir_1718_H3N2, "/position_coverage", sep = "")); getwd(); dir()
tsv <- dir(pattern="merged_position_coverage.tsv")
names_trunc <- gsub("_merged_position_coverage.tsv","",tsv)
n <- length(tsv)
list <- vector("list",n)
## Read all tables in tsv, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read_tsv(tsv[i], show_col_types = F, col_names = F)
  list[[i]] <- as.data.frame(list[[i]])
  if (length(list[[i]])==0) {print(paste("Empty data frame:", 
                                   names_trunc[i], sep = " "))}
  if (length(list[[i]])>0) {list[[i]]$WSLH_ID <- names_trunc[i]}
  names(list[[i]])[names(list[[i]]) == "X1"] <- "refseq"
  names(list[[i]])[names(list[[i]]) == "X2"] <- "POS"
  names(list[[i]])[names(list[[i]]) == "X3"] <- "COV"
  names(list) <- names_trunc}
## remove blank data frames
list <- Filter(function(x) dim(x)[1] > 0, list)
## separate columns
n <- length(list)
for (i in 1:n) {
  # separate the info column into its respective pieces
  list[[i]] <- separate(list[[i]],"refseq",
                        c("refseq", "EPI_ISL", "Gene_Name"),
                        sep="\\|",convert=FALSE)
  # factors
  list[[i]]$WSLH_ID <- as.factor(list[[i]]$WSLH_ID)
  list[[i]]$Gene_Name <- as.factor(list[[i]]$Gene_Name)
  # integers
  list[[i]]$POS <- as.integer(list[[i]]$POS)
  list[[i]]$COV <- as.integer(list[[i]]$COV)
}
## list to df
df_1718_H3N2 <- Reduce(full_join,list)

## remove fails from both .tsv files
percov_1718_H3N2 <- percov_1718_H3N2[percov_1718_H3N2$WSLH_ID %!in% to_remove$WSLH_ID,]
bstats_1718_H3N2 <- separate(bstats_1718_H3N2, "FileGene", 
                             c("WSLH_ID", "merged", "Gene"), 
                             sep = "_", convert = F)
bstats_1718_H3N2$FilgeGene <- paste(bstats_1718_H3N2$WSLH_ID, bstats_1718_H3N2$merged, bstats_1718_H3N2$Gene, sep = "_")
bstats_1718_H3N2 <- bstats_1718_H3N2[bstats_1718_H3N2$WSLH_ID %!in% to_remove$WSLH_ID,]
df_1718_H3N2     <- df_1718_H3N2[df_1718_H3N2$WSLH_ID %!in% to_remove$WSLH_ID,]


#### Import 1819_H3N2 ####
## create dir
setwd(dir_1819_H3N2); getwd(); dir()
percov_1819_H3N2 <- read.table("percent_coverage.tsv", sep = " ", header = FALSE)
colnames(percov_1819_H3N2) <- c("WSLH_ID", "Per_Coverage")
percov_1819_H3N2 <- percov_1819_H3N2[percov_1819_H3N2$WSLH_ID!="FileGene",]
bstats_1819_H3N2 <- as.data.frame(read_tsv("bam_stats.tsv"))
to_remove <- as.data.frame(read_tsv("to_remove.tsv"))

## create list for percov tsvs
setwd(paste(dir_1819_H3N2, "/position_coverage", sep = "")); getwd(); dir()
tsv <- dir(pattern="merged_position_coverage.tsv")
names_trunc <- gsub("_merged_position_coverage.tsv","",tsv)
n <- length(tsv)
list <- vector("list",n)
## Read all tables in tsv, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read_tsv(tsv[i], show_col_types = F, col_names = F)
  list[[i]] <- as.data.frame(list[[i]])
  if (length(list[[i]])==0) {print(paste("Empty data frame:", 
                                         names_trunc[i], sep = " "))}
  if (length(list[[i]])>0) {list[[i]]$WSLH_ID <- names_trunc[i]}
  names(list[[i]])[names(list[[i]]) == "X1"] <- "refseq"
  names(list[[i]])[names(list[[i]]) == "X2"] <- "POS"
  names(list[[i]])[names(list[[i]]) == "X3"] <- "COV"
  names(list) <- names_trunc}
## remove blank data frames
list <- Filter(function(x) dim(x)[1] > 0, list)
## separate columns
n <- length(list)
for (i in 1:n) {
  # separate the info column into its respective pieces
  list[[i]] <- separate(list[[i]],"refseq",
                        c("refseq", "EPI_ISL", "Gene_Name"),
                        sep="\\|",convert=FALSE)
  # factors
  list[[i]]$WSLH_ID <- as.factor(list[[i]]$WSLH_ID)
  list[[i]]$Gene_Name <- as.factor(list[[i]]$Gene_Name)
  # integers
  list[[i]]$POS <- as.integer(list[[i]]$POS)
  list[[i]]$COV <- as.integer(list[[i]]$COV)
}
## list to df
df_1819_H3N2 <- Reduce(full_join,list)

## remove fails from both .tsv files
percov_1819_H3N2 <- percov_1819_H3N2[percov_1819_H3N2$WSLH_ID %!in% to_remove$WSLH_ID,]
bstats_1819_H3N2 <- separate(bstats_1819_H3N2, "FileGene", 
                             c("WSLH_ID", "merged", "Gene"), 
                             sep = "_", convert = F)
bstats_1819_H3N2$FilgeGene <- paste(bstats_1819_H3N2$WSLH_ID, bstats_1819_H3N2$merged, bstats_1819_H3N2$Gene, sep = "_")
bstats_1819_H3N2 <- bstats_1819_H3N2[bstats_1819_H3N2$WSLH_ID %!in% to_remove$WSLH_ID,]
df_1819_H3N2     <- df_1819_H3N2[df_1819_H3N2$WSLH_ID %!in% to_remove$WSLH_ID,]


#### Import 1718_H1N1 ####
## create dir
setwd(dir_1718_H1N1); getwd(); dir()
percov_1718_H1N1 <- read.table("percent_coverage.tsv", sep = " ", header = FALSE)
colnames(percov_1718_H1N1) <- c("WSLH_ID", "Per_Coverage")
percov_1718_H1N1 <- percov_1718_H1N1[percov_1718_H1N1$WSLH_ID!="FileGene",]
bstats_1718_H1N1 <- as.data.frame(read_tsv("bam_stats.tsv"))
to_remove <- as.data.frame(read_tsv("to_remove.tsv"))

## create list for percov tsvs
setwd(paste(dir_1718_H1N1, "/position_coverage", sep = "")); getwd(); dir()
tsv <- dir(pattern="merged_position_coverage.tsv")
names_trunc <- gsub("_merged_position_coverage.tsv","",tsv)
n <- length(tsv)
list <- vector("list",n)
## Read all tables in tsv, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read_tsv(tsv[i], show_col_types = F, col_names = F)
  list[[i]] <- as.data.frame(list[[i]])
  if (length(list[[i]])==0) {print(paste("Empty data frame:", 
                                         names_trunc[i], sep = " "))}
  if (length(list[[i]])>0) {list[[i]]$WSLH_ID <- names_trunc[i]}
  names(list[[i]])[names(list[[i]]) == "X1"] <- "refseq"
  names(list[[i]])[names(list[[i]]) == "X2"] <- "POS"
  names(list[[i]])[names(list[[i]]) == "X3"] <- "COV"
  names(list) <- names_trunc}
## remove blank data frames
list <- Filter(function(x) dim(x)[1] > 0, list)
## separate columns
n <- length(list)
for (i in 1:n) {
  # separate the info column into its respective pieces
  list[[i]] <- separate(list[[i]],"refseq",
                        c("refseq", "EPI_ISL", "Gene_Name"),
                        sep="\\|",convert=FALSE)
  # factors
  list[[i]]$WSLH_ID <- as.factor(list[[i]]$WSLH_ID)
  list[[i]]$Gene_Name <- as.factor(list[[i]]$Gene_Name)
  # integers
  list[[i]]$POS <- as.integer(list[[i]]$POS)
  list[[i]]$COV <- as.integer(list[[i]]$COV)
}
## list to df
df_1718_H1N1 <- Reduce(full_join,list)

## remove fails from both .tsv files
percov_1718_H1N1 <- percov_1718_H1N1[percov_1718_H1N1$WSLH_ID %!in% to_remove$WSLH_ID,]
bstats_1718_H1N1 <- separate(bstats_1718_H1N1, "FileGene", 
                             c("WSLH_ID", "merged", "Gene"), 
                             sep = "_", convert = F)
bstats_1718_H1N1$FilgeGene <- paste(bstats_1718_H1N1$WSLH_ID, bstats_1718_H1N1$merged, bstats_1718_H1N1$Gene, sep = "_")
bstats_1718_H1N1 <- bstats_1718_H1N1[bstats_1718_H1N1$WSLH_ID %!in% to_remove$WSLH_ID,]
df_1718_H1N1     <- df_1718_H1N1[df_1718_H1N1$WSLH_ID %!in% to_remove$WSLH_ID,]


#### Import 1819_H1N1 ####
## create dir
setwd(dir_1819_H1N1); getwd(); dir()
percov_1819_H1N1 <- read.table("percent_coverage.tsv", sep = " ", header = FALSE)
colnames(percov_1819_H1N1) <- c("WSLH_ID", "Per_Coverage")
percov_1819_H1N1 <- percov_1819_H1N1[percov_1819_H1N1$WSLH_ID!="FileGene",]
bstats_1819_H1N1 <- as.data.frame(read_tsv("bam_stats.tsv"))
to_remove <- as.data.frame(read_tsv("to_remove.tsv"))

## create list for percov tsvs
setwd(paste(dir_1819_H1N1, "/position_coverage", sep = "")); getwd(); dir()
tsv <- dir(pattern="merged_position_coverage.tsv")
names_trunc <- gsub("_merged_position_coverage.tsv","",tsv)
n <- length(tsv)
list <- vector("list",n)
## Read all tables in tsv, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read_tsv(tsv[i], show_col_types = F, col_names = F)
  list[[i]] <- as.data.frame(list[[i]])
  if (length(list[[i]])==0) {print(paste("Empty data frame:", 
                                         names_trunc[i], sep = " "))}
  if (length(list[[i]])>0) {list[[i]]$WSLH_ID <- names_trunc[i]}
  names(list[[i]])[names(list[[i]]) == "X1"] <- "refseq"
  names(list[[i]])[names(list[[i]]) == "X2"] <- "POS"
  names(list[[i]])[names(list[[i]]) == "X3"] <- "COV"
  names(list) <- names_trunc}
## remove blank data frames
list <- Filter(function(x) dim(x)[1] > 0, list)
## separate columns
n <- length(list)
for (i in 1:n) {
  # separate the info column into its respective pieces
  list[[i]] <- separate(list[[i]],"refseq",
                        c("refseq", "EPI_ISL", "Gene_Name"),
                        sep="\\|",convert=FALSE)
  # factors
  list[[i]]$WSLH_ID <- as.factor(list[[i]]$WSLH_ID)
  list[[i]]$Gene_Name <- as.factor(list[[i]]$Gene_Name)
  # integers
  list[[i]]$POS <- as.integer(list[[i]]$POS)
  list[[i]]$COV <- as.integer(list[[i]]$COV)
}
## list to df
df_1819_H1N1 <- Reduce(full_join,list)

## remove fails from both .tsv files
percov_1819_H1N1 <- percov_1819_H1N1[percov_1819_H1N1$WSLH_ID %!in% to_remove$WSLH_ID,]
bstats_1819_H1N1 <- separate(bstats_1819_H1N1, "FileGene", 
                             c("WSLH_ID", "merged", "Gene"), 
                             sep = "_", convert = F)
bstats_1819_H1N1$FilgeGene <- paste(bstats_1819_H1N1$WSLH_ID, bstats_1819_H1N1$merged, bstats_1819_H1N1$Gene, sep = "_")
bstats_1819_H1N1 <- bstats_1819_H1N1[bstats_1819_H1N1$WSLH_ID %!in% to_remove$WSLH_ID,]
df_1819_H1N1     <- df_1819_H1N1[df_1819_H1N1$WSLH_ID %!in% to_remove$WSLH_ID,]

#### Plot 1718_H3N2 ####
## 1718_H3N2
x <- as.data.frame(ds(df_1718_H3N2, 
                      varname = "COV", 
                      groupnames = c("Gene_Name", "POS")))
genes <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")
x_axis_start <- c(0,0,0,0,0,0,0,0)
x_axis_end <- c(2295,2286,2163,1718,1530,1415,995,865)
by <- c(1000,1000,1000,750,750,500,400,500)
spec <- data.frame(genes, x_axis_start, x_axis_end, by)
plots <- list()
for (gene in genes) {
  df_plot     <- df_1718_H3N2[df_1718_H3N2$Gene_Name==gene,]
  trendline   <- x[x$Gene_Name==gene,]
  x_axis_start<- spec$x_axis_start[spec$genes==gene]
  x_axis_plot <- spec$x_axis_end[spec$genes==gene]
  by_plot     <- spec$by[spec$genes==gene]
  name_plot   <- paste(gene,"_","plot")
  name_plot   <- gsub(" ","",name_plot)
  
  ## y-axis formatting
  if (gene == "PB2"| gene == 'NP'){
    y_aesthetics = theme(axis.line.y = element_line(colour="grey"),
                         axis.text.y = element_text(hjust=0.5, size = 6),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))            
  } else {
    y_aesthetics = theme(axis.line.y=element_blank(),
                         axis.ticks.y= element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.y = element_blank())
  }
  ## legend formatting
  if (gene == "HA"){
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "right")            
  } else {
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "none")    
  }
  plot <- ggplot() + 
    geom_line(data = df_plot, 
              aes(x = POS, y = COV, 
                  group = WSLH_ID),
              alpha = 0.2, color = "black") + 
    geom_line(data = trendline, 
              aes(x = POS, y = median_COV), 
              alpha = 0.8, color = "red") + 
    labs(y = "Coverage", x = "", title=gene) + 
    coord_cartesian(ylim = c(0, 10000), xlim = c(x_axis_start, x_axis_plot)) + 
    scale_x_continuous(breaks = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot),
                       labels = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot)) + 
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 6),
          axis.line.x = element_line(colour="grey"),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 6, margin = margin(t = 4))) + 
    y_aesthetics + legend_aesthetics + legend_formatting
  plots[[name_plot]] <- plot
}
Fig1supp_top <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                          nrow = 1, rel_widths = c(0.32, 0.22, 0.21, 0.20))
Fig1supp_bot <- plot_grid(plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                          nrow = 1, rel_widths = c(0.35, 0.25, 0.25, 0.15))
Fig1asupp_1718_H3N2 <- plot_grid(Fig1supp_top, Fig1supp_bot, 
                                 nrow = 2, ncol = 1,
                                 labels = "A", label_size = Size_adjust, 
                                 hjust = LR_adjust, vjust = UD_adjust)
rm(Fig1supp_top)
rm(Fig1supp_bot)

#### Plot 1819_H3N2 ####
## 1819_H3N2
x <- as.data.frame(ds(df_1819_H3N2, 
                      varname = "COV", 
                      groupnames = c("Gene_Name", "POS")))
genes <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")
x_axis_start <- c(0,0,0,0,0,0,0,0)
x_axis_end <- c(2295,2286,2163,1718,1530,1415,995,865)
by <- c(1000,1000,1000,750,750,500,400,500)
spec <- data.frame(genes, x_axis_start, x_axis_end, by)
plots <- list()
for (gene in genes) {
  df_plot     <- df_1819_H3N2[df_1819_H3N2$Gene_Name==gene,]
  trendline   <- x[x$Gene_Name==gene,]
  x_axis_start<- spec$x_axis_start[spec$genes==gene]
  x_axis_plot <- spec$x_axis_end[spec$genes==gene]
  by_plot     <- spec$by[spec$genes==gene]
  name_plot   <- paste(gene,"_","plot")
  name_plot   <- gsub(" ","",name_plot)
  
  ## y-axis formatting
  if (gene == "PB2"| gene == 'NP'){
    y_aesthetics = theme(axis.line.y = element_line(colour="grey"),
                         axis.text.y = element_text(hjust=0.5, size = 6),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))            
  } else {
    y_aesthetics = theme(axis.line.y=element_blank(),
                         axis.ticks.y= element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.y = element_blank())
  }
  ## legend formatting
  if (gene == "HA"){
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "right")            
  } else {
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "none")    
  }
  plot <- ggplot() + 
    geom_line(data = df_plot, 
              aes(x = POS, y = COV, 
                  group = WSLH_ID),
              alpha = 0.2, color = "black") + 
    geom_line(data = trendline, 
              aes(x = POS, y = median_COV), 
              alpha = 0.8, color = "red") + 
    labs(y = "Coverage", x = "", title=gene) + 
    coord_cartesian(ylim = c(0, 10000), xlim = c(x_axis_start, x_axis_plot)) + 
    scale_x_continuous(breaks = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot),
                       labels = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot)) + 
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 6),
          axis.line.x = element_line(colour="grey"),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 6, margin = margin(t = 4))) + 
    y_aesthetics + legend_aesthetics + legend_formatting
  plots[[name_plot]] <- plot
}
Fig1supp_top <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                          nrow = 1, rel_widths = c(0.32, 0.22, 0.21, 0.20))
Fig1supp_bot <- plot_grid(plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                          nrow = 1, rel_widths = c(0.35, 0.25, 0.25, 0.15))
Fig1asupp_1819_H3N2 <- plot_grid(Fig1supp_top, Fig1supp_bot, 
                                 nrow = 2, ncol = 1,
                                 labels = "B", label_size = Size_adjust, 
                                 hjust = LR_adjust, vjust = UD_adjust)
rm(Fig1supp_top)
rm(Fig1supp_bot)

#### Plot 1718_H1N1 ####
## 1819_H3N2
x <- as.data.frame(ds(df_1718_H1N1, 
                      varname = "COV", 
                      groupnames = c("Gene_Name", "POS")))
genes <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")
x_axis_start <- c(0,0,0,0,0,0,0,0)
x_axis_end <- c(2295,2286,2163,1718,1530,1415,995,865)
by <- c(1000,1000,1000,750,750,500,400,500)
spec <- data.frame(genes, x_axis_start, x_axis_end, by)
plots <- list()
for (gene in genes) {
  df_plot     <- df_1718_H1N1[df_1718_H1N1$Gene_Name==gene,]
  trendline   <- x[x$Gene_Name==gene,]
  x_axis_start<- spec$x_axis_start[spec$genes==gene]
  x_axis_plot <- spec$x_axis_end[spec$genes==gene]
  by_plot     <- spec$by[spec$genes==gene]
  name_plot   <- paste(gene,"_","plot")
  name_plot   <- gsub(" ","",name_plot)
  
  ## y-axis formatting
  if (gene == "PB2"| gene == 'NP'){
    y_aesthetics = theme(axis.line.y = element_line(colour="grey"),
                         axis.text.y = element_text(hjust=0.5, size = 6),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))            
  } else {
    y_aesthetics = theme(axis.line.y=element_blank(),
                         axis.ticks.y= element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.y = element_blank())
  }
  ## legend formatting
  if (gene == "HA"){
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "right")            
  } else {
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "none")    
  }
  plot <- ggplot() + 
    geom_line(data = df_plot, 
              aes(x = POS, y = COV, 
                  group = WSLH_ID),
              alpha = 0.2, color = "black") + 
    geom_line(data = trendline, 
              aes(x = POS, y = median_COV), 
              alpha = 0.8, color = "red") + 
    labs(y = "Coverage", x = "", title=gene) + 
    coord_cartesian(ylim = c(0, 10000), xlim = c(x_axis_start, x_axis_plot)) + 
    scale_x_continuous(breaks = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot),
                       labels = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot)) + 
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 6),
          axis.line.x = element_line(colour="grey"),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 6, margin = margin(t = 4))) + 
    y_aesthetics + legend_aesthetics + legend_formatting
  plots[[name_plot]] <- plot
}
Fig1supp_top <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                          nrow = 1, rel_widths = c(0.32, 0.22, 0.21, 0.20))
Fig1supp_bot <- plot_grid(plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                          nrow = 1, rel_widths = c(0.35, 0.25, 0.25, 0.15))
Fig1asupp_1718_H1N1 <- plot_grid(Fig1supp_top, Fig1supp_bot, 
                                 nrow = 2, ncol = 1,
                                 labels = "C", label_size = Size_adjust, 
                                 hjust = LR_adjust, vjust = UD_adjust)
rm(Fig1supp_top)
rm(Fig1supp_bot)

#### Plot 1819_H1N1 ####
## 1819_H3N2
x <- as.data.frame(ds(df_1819_H1N1, 
                      varname = "COV", 
                      groupnames = c("Gene_Name", "POS")))
genes <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")
x_axis_start <- c(0,0,0,0,0,0,0,0)
x_axis_end <- c(2295,2286,2163,1718,1530,1415,995,865)
by <- c(1000,1000,1000,750,750,500,400,500)
spec <- data.frame(genes, x_axis_start, x_axis_end, by)
plots <- list()
for (gene in genes) {
  df_plot     <- df_1819_H1N1[df_1819_H1N1$Gene_Name==gene,]
  trendline   <- x[x$Gene_Name==gene,]
  x_axis_start<- spec$x_axis_start[spec$genes==gene]
  x_axis_plot <- spec$x_axis_end[spec$genes==gene]
  by_plot     <- spec$by[spec$genes==gene]
  name_plot   <- paste(gene,"_","plot")
  name_plot   <- gsub(" ","",name_plot)
  
  ## y-axis formatting
  if (gene == "PB2"| gene == 'NP'){
    y_aesthetics = theme(axis.line.y = element_line(colour="grey"),
                         axis.text.y = element_text(hjust=0.5, size = 6),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))            
  } else {
    y_aesthetics = theme(axis.line.y=element_blank(),
                         axis.ticks.y= element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.y = element_blank())
  }
  ## legend formatting
  if (gene == "HA"){
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "right")            
  } else {
    legend_aesthetics = theme(legend.background = element_blank(),
                              legend.title = element_blank(),
                              legend.key = element_blank(),
                              legend.position = "none")    
  }
  plot <- ggplot() + 
    geom_line(data = df_plot, 
              aes(x = POS, y = COV, 
                  group = WSLH_ID),
              alpha = 0.2, color = "black") + 
    geom_line(data = trendline, 
              aes(x = POS, y = median_COV), 
              alpha = 0.8, color = "red") + 
    labs(y = "Coverage", x = "", title=gene) + 
    coord_cartesian(ylim = c(0, 10000), xlim = c(x_axis_start, x_axis_plot)) + 
    scale_x_continuous(breaks = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot),
                       labels = seq(signif(x_axis_start,3),
                                    signif(x_axis_plot,3),by_plot)) + 
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 6),
          axis.line.x = element_line(colour="grey"),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 6, margin = margin(t = 4))) + 
    y_aesthetics + legend_aesthetics + legend_formatting
  plots[[name_plot]] <- plot
}
Fig1supp_top <- plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
                          nrow = 1, rel_widths = c(0.32, 0.22, 0.21, 0.20))
Fig1supp_bot <- plot_grid(plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
                          nrow = 1, rel_widths = c(0.35, 0.25, 0.25, 0.15))
Fig1asupp_1819_H1N1 <- plot_grid(Fig1supp_top, Fig1supp_bot, 
                                 nrow = 2, ncol = 1,
                                 labels = "D", label_size = Size_adjust, 
                                 hjust = LR_adjust, vjust = UD_adjust)
rm(Fig1supp_top)
rm(Fig1supp_bot)

#### Plot all ####
## all together now!
Fig1supp_AB <- plot_grid(Fig1asupp_1718_H3N2, Fig1asupp_1819_H3N2,
                         nrow = 1, ncol = 2)
Fig1supp_CD <- plot_grid(Fig1asupp_1718_H1N1, Fig1asupp_1819_H1N1,
                         nrow = 1, ncol = 2)

Fig1supp <- plot_grid(Fig1supp_AB, Fig1supp_CD,
                      nrow = 2, ncol = 1)

#### Save ####
setwd(dir_s); getwd(); dir()
ggsave("Fig1supp_coverage.pdf", Fig1supp,
       width = 8, height = 11, 
       units = "in", dpi = 320)
