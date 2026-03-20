## Clear Global Environment
rm(list = ls())

## Values
minAF <- .03
minCov <- 100

## README
# Run each season_subtype through this pipeline
# Input: recalled/recall_fn_ann
# Input: recalled/snpgenie/
# Output: dfs and three plots

## Directories
## 17-18 H3N2
#dir    <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H3N2/recalled/recall_fn_ann",sep="")
#dir_sg <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H3N2/recalled/snpgenie",sep="")
#dir_s  <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H3N2/dfs",sep="")

## 18-19 H3N2
#dir    <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H3N2/recalled/recall_fn_ann",sep="")
#dir_sg <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H3N2/recalled/snpgenie",sep="")
#dir_s  <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H3N2/dfs",sep="")

## 17-18 H1N1
#dir    <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H1N1/recalled/recall_fn_ann",sep="")
#dir_sg <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H1N1/recalled/snpgenie",sep="")
#dir_s  <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H1N1/dfs",sep="")

## 18-19 H1N1
dir    <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H1N1/recalled/recall_fn_ann",sep="")
dir_sg <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H1N1/recalled/snpgenie",sep="")
dir_s  <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H1N1/dfs",sep="")

# metadata
dir_md <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/summary",sep="")

# neutral_expectations
dir_ne <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/code/neutral_expectations",sep="")

#### Session prep ####
## Install packages and load libraries as required
if(!require(tidyverse)){
  install.packages("tidyverse",dependencies = T)
  library(tidyverse)
}
if(!require(lubridate)){
  install.packages("lubridate",dependencies = T)
  library(lubridate)
}
if(!require(ggplot2)){
  install.packages("ggplot2",dependencies = T)
  library(ggplot2)
}
if(!require(ggrepel)){
  install.packages("ggrepel",dependencies = T)
  library(ggrepel)
}
if(!require(grid)){
  install.packages("grid",dependencies = T)
  library(grid)
}
if(!require(gridExtra)){
  install.packages("gridExtra",dependencies = T)
  library(gridExtra)
}
if(!require(cowplot)){
  install.packages("cowplot")
  library(cowplot)
}
if(!require(reshape2)){
  install.packages("reshape2",dependencies = T)
  library(reshape2)
}
if(!require(ggsignif)){
  install.packages("ggsignif",dependencies = T)
  library(ggsignif)
}
if(!require(vcfR)){
  install.packages("vcfR")
  library(vcfR)
}

#### Functions ####
## standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}

## data summary
ds <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sum(sd(x[[col]]) / sqrt(length(x[[1]]))))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
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


#### Import ####
## create dir
setwd(dir); getwd(); dir()

## create list for vcf(s)
vcf <- dir(pattern=".vcf")
names_trunc <- gsub("_recall_fn_ann.vcf","",vcf)
n <- length(vcf)
list <- vector("list",n)

## Read all tables in vcf, apply to list, change columns
for (i in 1:n) {
  list[[i]] <- read.vcfR(vcf[i], verbose=T)
  list[[i]] <- as.data.frame(list[[i]]@fix)
  names(list) <- names_trunc}

## remove blank data frames (there should be none if QC was performed!)
list <- Filter(function(x) dim(x)[1] > 0, list)

## separate columns
n <- length(list)
for (i in 1:n) {
  # separate the info column into its respective pieces
  list[[i]] <- separate(list[[i]],"INFO",c("SN","STA","STO","TYP",
                                           "R1P","R1M","R2P","R2M",
                                           "AD","DP","MCO","PPC",
                                           "AF","RAF","LS","MQS",
                                           "MQM","BQS","BQM","EDS",
                                           "EDM","IDS","IDM","NVC",
                                           "FLG","CED","HMP","SB",
                                           "snpEff"),
                        sep=";",convert=FALSE)
  # separate snpEff column into its respective pieces
  list[[i]] <- separate(list[[i]],"snpEff",c("Allele","Ann","Ann_Impact","Gene_Name",
                                             "Gene_ID","Feature_Type","Feature_ID","Transcript_BioType",
                                             "Rank","HGVS.c","HGVS.p"),
                        sep="\\|",convert=FALSE)
  # separate the leftover info column into the label (discarded) and value
  list[[i]] <- separate(list[[i]],"AF",c("AF_label","AF"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"AD",c("AD_label","AD"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"DP",c("DP_label","DP"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"STA",c("STA_label","STA"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"STO",c("STO_label","STO"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"R1P",c("R1P_label","R1P"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"R1M",c("R1M_label","R1M"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"NVC",c("NVC_label","NVC"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"FLG",c("FLG_label","FLG"),sep="=",convert=FALSE)  
  list[[i]] <- separate(list[[i]],"HMP",c("HMP_label","HMP"),sep="=",convert=FALSE)  
  # remove columns by string in name
  list[[i]] <- list[[i]] %>% select(-contains(c("_label","TYP","MCO","PPC","RAF","LS",
                                                "MQS","MQM","BQS","BQM","EDS","EDM",
                                                "IDS","IDM","CED","SB","R2P","R2M","SN")))
  # create WSLH_ID column
  list[[i]]$WSLH_ID <- paste(names_trunc[i])
  # factors
  list[[i]]$CHROM <- as.factor(list[[i]]$CHROM)
  list[[i]]$WSLH_ID <- as.factor(list[[i]]$WSLH_ID)
  list[[i]]$Gene_Name <- as.factor(list[[i]]$Gene_Name)
  list[[i]]$Gene_ID <- as.factor(list[[i]]$Gene_ID)
  list[[i]]$Ann <- as.factor(list[[i]]$Ann)
  list[[i]]$FILTER <- as.factor(list[[i]]$FILTER)
  # integers
  list[[i]]$POS <- as.integer(list[[i]]$POS)
  list[[i]]$ID <- as.integer(list[[i]]$ID)
  list[[i]]$STA <- as.integer(list[[i]]$STA)
  list[[i]]$STO <- as.integer(list[[i]]$STO)
  list[[i]]$R1P <- as.integer(list[[i]]$R1P)
  list[[i]]$R1M <- as.integer(list[[i]]$R1M)
  list[[i]]$AD <- as.integer(list[[i]]$AD)
  list[[i]]$DP <- as.integer(list[[i]]$DP)
  list[[i]]$NVC <- as.integer(list[[i]]$NVC)
  list[[i]]$FLG <- as.integer(list[[i]]$FLG)
  list[[i]]$HMP <- as.integer(list[[i]]$HMP)
  # numerics
  list[[i]]$AF <- as.numeric(list[[i]]$AF)
  list[[i]]$QUAL <- as.numeric(list[[i]]$QUAL)
}

## import SNPGenie data for later
setwd(dir_sg); getwd(); dir()
sg <- as.data.frame(read_tsv("snpgenie_compiled.tsv"))
sg_sw <- as.data.frame(read_tsv("snpgenie_compiled_sliding_window.tsv"))
sg_to_remove <- read_tsv("to_remove.tsv")
sg_sw_antigenic <- read_tsv("antigenic_regions.tsv")
sg_sw_1 <- read_tsv("snpgenie_compiled_sliding_window_1.tsv")

## import metadata for later
setwd(dir_md); getwd(); dir()
md <- as.data.frame(read_csv("df_samplelist_metadata.csv"))

## import neutral expectations table for iSNV frequency spectra
# generated from ~/code/neutral_expectations/neutral_expectations.ipynb
setwd(dir_ne); getwd(); dir()
df_mut_bins_neutral <- as.data.frame(read_csv("neutral_df.csv"))

#### md ####
## remove index cases and donor second time point
md$WSLH_ID_trunc <- md$WSLH_ID
md_donors <- md[md$Case=="Index",] %>% 
  group_by(Study_ID) %>%
  top_n(-1, WSLH_ID_trunc)

## extract recipients and bind single-donors
md_recipients <- md[md$Case=="Contact",]
md <- rbind(md_donors, md_recipients)
md <- md[order(as.integer(md$Num)),]

## Index to donor
md$Case <- gsub("Index", "Donor", md$Case)
md$Case_Num <- gsub("Index", "Donor", md$Case_Num)

## Contact to recipient
md$Case <- gsub("Contact", "Recipient", md$Case)
md$Case_Num <- gsub("Contact", "Recipient", md$Case_Num)

## Vaccination status
md$flu_vaccine <- gsub("No", "Unvaccinated", md$flu_vaccine)
md$flu_vaccine <- gsub("Yes", "Vaccinated", md$flu_vaccine)

#### df_all_iSNVs ####
## merge data frames from list and change class
df_all_iSNVs <- Reduce(full_join,list)

## remove iSNVs in the non-coding regions, where a protein is not encoded
df_all_iSNVs <- df_all_iSNVs[df_all_iSNVs$HGVS.p!="",]

## factor
df_all_iSNVs$CHROM <- as.factor(df_all_iSNVs$CHROM)
df_all_iSNVs$FILTER <- as.factor(df_all_iSNVs$FILTER)
df_all_iSNVs$WSLH_ID <- as.factor(df_all_iSNVs$WSLH_ID)
df_all_iSNVs$Ann_Impact <- as.factor(df_all_iSNVs$Ann_Impact)
## numeric
df_all_iSNVs$POS <- as.numeric(df_all_iSNVs$POS)
df_all_iSNVs$QUAL <- as.numeric(df_all_iSNVs$QUAL)
df_all_iSNVs$STA <- as.numeric(df_all_iSNVs$STA)
df_all_iSNVs$STO <- as.numeric(df_all_iSNVs$STO)
df_all_iSNVs$R1P <- as.numeric(df_all_iSNVs$R1P)
df_all_iSNVs$R1M <- as.numeric(df_all_iSNVs$R1M)
df_all_iSNVs$AD <- as.numeric(df_all_iSNVs$AD)
df_all_iSNVs$DP <- as.numeric(df_all_iSNVs$DP)
df_all_iSNVs$AF <- as.numeric(df_all_iSNVs$AF)
df_all_iSNVs$NVC <- as.numeric(df_all_iSNVs$NVC)
df_all_iSNVs$FLG <- as.numeric(df_all_iSNVs$FLG)
df_all_iSNVs$HMP <- as.numeric(df_all_iSNVs$HMP)

## Apply antigenicity annotations
df_all_iSNVs$Antigenicity <- "Non-antigenic"
df_all_iSNVs$Antigenic_site <- ""

# Which rows of df_all_iSNVs are HA iSNVs?
str_rows_HA <- which(matrix(grepl("HA", df_all_iSNVs$Feature_ID)), arr.ind=T)[,1]
for (i in 1:length(str_rows_HA)) {
  # cycle through all HA iSNVs, row by row
  # the ith sample of str_rows_HA is not the ith row of df_al_iSNVs!
  for (j in 1:length(sg_sw_antigenic$Start)) {
    # cycle through annotations for each iSNV, row by row
    if(df_all_iSNVs$POS[str_rows_HA[i]]>=sg_sw_antigenic$Start[j] & 
       #if iSNV i is positioned at or after the START of the antigenic site j
       df_all_iSNVs$POS[str_rows_HA[i]]<=sg_sw_antigenic$Stop[j]) { 
       #if iSNV i is positioned at or before the END of antigenic site j
      df_all_iSNVs$Antigenicity[str_rows_HA[i]] <- sg_sw_antigenic$ID[j] 
      #if both of those statements are satisfied, transfer the annotation
      df_all_iSNVs$Antigenic_site[str_rows_HA[i]] <- sg_sw_antigenic$Antigenic_Site[j] 
      #also transfer the antigenic site classification
    }
  }
}

# factor
df_all_iSNVs$Antigenicity <- as.factor(df_all_iSNVs$Antigenicity)
df_all_iSNVs$Antigenic_site <- as.factor(df_all_iSNVs$Antigenic_site)

# remove NA values
df_all_iSNVs$Antigenic_site[is.na(df_all_iSNVs$Antigenic_site)] <- ""

# create Gene_ID
df_all_iSNVs$Gene_ID <- paste(df_all_iSNVs$Gene_Name,df_all_iSNVs$Antigenicity,
                              df_all_iSNVs$Antigenic_site, sep = "_")

## factor Gene_ID
df_all_iSNVs$Gene_ID <- as.factor(df_all_iSNVs$Gene_ID)

## Add GROUP, which will be used to sort data
df_all_iSNVs$WSLH_ID <- df_all_iSNVs$FILTER
df_all_iSNVs$GROUP <- paste(df_all_iSNVs$WSLH_ID, df_all_iSNVs$Gene_ID, sep = "_")
df_all_iSNVs$GROUP <- as.factor(df_all_iSNVs$GROUP)

## Remove unwanted genes (MP and NS segments, PA-X and PB1-F2 genes)
str_unwanted_genes <- c("M1", "M2", "NEP", "NS1", "PA-X", "PB1-F2")
df_all_iSNVs <- df_all_iSNVs[df_all_iSNVs$Gene_Name %!in% str_unwanted_genes,]
# re-factor Gene_Name
df_all_iSNVs$Gene_Name <- as.factor(as.character(df_all_iSNVs$Gene_Name))

## identify samples that are in the top 10% for number of within-host variants above .5% AF
df_all_iSNVs <- df_all_iSNVs[df_all_iSNVs$AF>.005,]
#plot(table(df_all_iSNVs$WSLH_ID), main = "Before removal of top 10%",
#     xlab = "Sample ID", ylab = "Number of iSNVs")
iSNV_table <- as.data.frame(table(df_all_iSNVs$WSLH_ID))
# order from highest to lowest
iSNV_table <- iSNV_table[order(iSNV_table$Freq, decreasing = T),]
n_samples <- length(iSNV_table$Var1)
top_ten_percent_iSNVs <- head(iSNV_table, round(n_samples*.1,0))
#plot(iSNV_table, main = "No removal of top 10% or AF>=.03",
#     xlab = "Sample ID", ylab = "Number of iSNVs")
# remove those samples
df_all_iSNVs <- df_all_iSNVs[df_all_iSNVs$WSLH_ID %!in% top_ten_percent_iSNVs[,1],]
# filter by minAF and minCov
df_all_iSNVs <- filter(df_all_iSNVs, AF>=minAF)
df_all_iSNVs <- filter(df_all_iSNVs, DP>=minCov)
# plot the revised counts
iSNV_table_rev <- as.data.frame(table(df_all_iSNVs$WSLH_ID))
iSNV_table_rev <- as.data.frame(iSNV_table_rev[order(iSNV_table_rev$Freq, decreasing = T),])
#plot(iSNV_table_rev, main = "After removal of top 10% and AF>=.03",
#     xlab = "Sample ID", ylab = "Number of iSNVs")

#quantile(iSNV_table$Freq, c(0.05))
#quantile(iSNV_table$Freq, c(0.95))

#plot(table(df_all_iSNVs$WSLH_ID), main = "Only AF>=.03", xlab = "Sample ID", ylab = "Number of iSNVs")

#### sg ####
## format snpgenie file
# file name (will be WSLH_ID)
sg$file <- gsub("./segments/", "", sg$file)
sg$file <- gsub("_recall_fn_ann_", "", sg$file)
sg$file <- gsub(".vcf", "", sg$file)
sg$file <- gsub(".*/", "", sg$file)
sg$product <- gsub("Neuraminidase", "NA", sg$product)
sg$file <- gsub("Neuraminidase", "NA", sg$file)
sg$file <- gsub("HA", "", sg$file)
sg$file <- gsub("M1", "", sg$file)
sg$file <- gsub("M2", "", sg$file)
sg$file <- gsub("NA", "", sg$file)
sg$file <- gsub("NEP", "", sg$file)
sg$file <- gsub("NP", "", sg$file)
sg$file <- gsub("NS1", "", sg$file)
sg$file <- gsub("PA", "", sg$file)
sg$file <- gsub("PB1", "", sg$file)
sg$file <- gsub("PB2", "", sg$file)
sg$file <- gsub("NS", "", sg$file)
sg$file <- gsub("MP", "", sg$file)

sg_sw$file <- gsub("./segments/", "", sg_sw$file)
sg_sw$file <- gsub("_recall_fn_ann_", "", sg_sw$file)
sg_sw$file <- gsub(".vcf", "", sg_sw$file)
sg_sw$file <- gsub(".*/", "", sg_sw$file)
sg_sw$file <- gsub("HA", "", sg_sw$file)

sg_sw_1$file <- gsub("./segments/", "", sg_sw_1$file)
sg_sw_1$file <- gsub("_recall_fn_ann_", "", sg_sw_1$file)
sg_sw_1$file <- gsub(".vcf", "", sg_sw_1$file)
sg_sw_1$file <- gsub(".*/", "", sg_sw_1$file)
sg_sw_1$file <- gsub("HA", "", sg_sw_1$file)

# filter out QC_fails from sg, sg_sw, and sg_sw_1
sg <- sg[sg$file %!in% sg_to_remove$WSLH_ID,]
sg_sw <- sg_sw[sg_sw$file %!in% sg_to_remove$WSLH_ID,]
sg_sw_1 <- sg_sw_1[sg_sw_1$file %!in% sg_to_remove$WSLH_ID,]

## remove highly divergent samples from list and sg analyses
top_ten_percent_iSNVs$Var1
list <- list[names(list) %!in% top_ten_percent_iSNVs[,1]]
sg <- sg[sg$file %!in% top_ten_percent_iSNVs[,1],]
sg_sw <- sg_sw[sg_sw$file %!in% top_ten_percent_iSNVs[,1],]
sg_sw_1 <- sg_sw_1[sg_sw_1$file %!in% top_ten_percent_iSNVs[,1],]

# factor
levels(as.factor(as.character(sg$product)))
sg$product <- as.factor(sg$product)
sg$file <- as.factor(sg$file)
levels(as.factor(as.character(sg_sw$product)))
sg_sw$product <- as.factor(sg_sw$product)
sg_sw$file <- as.factor(sg_sw$file)
sg_sw_1$file <- as.factor(sg_sw_1$file)

## calculate
sg$piN <- as.numeric(sg$piN)
sg$piS <- as.numeric(sg$piS)
sg <- transform(sg, 
                piN_minus_piS      = sg$piN - sg$piS, 
                #piN.confirm = sg$N_diffs / sg$N_sites, ## should equal piN in table
                #piS.confirm = sg$S_diffs / sg$S_sites, ## should equal piS in table
                pi                 = sg$piN + sg$piS)

## GROUP column (WSLH_ID and gene_name)
sg$GROUP <- paste(sg$file, sg$product, sep="-")


## remove unwanted genes from sg
sg <- sg[sg$product %!in% str_unwanted_genes,]

#### ds_sg_sw ####
## prelim sw plots
# calculate piN-piS
sg_sw <- transform(sg_sw, piN_minus_piS=piN-piS)
ds_sg_sw_piN <- ds(sg_sw[sg_sw$product=="HA",], varname="piN", groupnames=c("first_site"))
ds_sg_sw_piS <- ds(sg_sw[sg_sw$product=="HA",], varname="piS", groupnames=c("first_site"))
ds_sg_sw_piNpiS <- ds(sg_sw[sg_sw$product=="HA",], varname="piN_minus_piS", groupnames=c("first_site"))
ds_sg_sw <- merge(ds_sg_sw_piN, ds_sg_sw_piS, by = "first_site")
ds_sg_sw <- merge(ds_sg_sw,     ds_sg_sw_piNpiS, by = "first_site")

## molten means
molten_ds_sg_sw <- melt(ds_sg_sw,
                        id.vars = c("first_site"),
                        measure.vars = c("piN", "piS", "piN_minus_piS"),
                        variable.name = "Pi")
colnames(molten_ds_sg_sw)[2] <- "Pi1"
colnames(molten_ds_sg_sw)[3] <- "Mean"
molten_ds_sg_sw$Pi1 <- gsub("piN", "Nonsynonymous", molten_ds_sg_sw$Pi1)
molten_ds_sg_sw$Pi1 <- gsub("piS", "Synonymous", molten_ds_sg_sw$Pi1)

## molten SEs
molten_ds_sg_sw_se <- melt(ds_sg_sw,
                          id.vars = c("first_site"),
                          measure.vars = c("se.x", "se.y", "se"),
                          variable.name = "SE")
colnames(molten_ds_sg_sw_se)[1] <- "first_site2"
colnames(molten_ds_sg_sw_se)[2] <- "Pi"
colnames(molten_ds_sg_sw_se)[3] <- "SE"
molten_ds_sg_sw_se$Pi <- gsub("se.x", "Nonsynonymous", molten_ds_sg_sw_se$Pi)
molten_ds_sg_sw_se$Pi <- gsub("se.y", "Synonymous",    molten_ds_sg_sw_se$Pi)
molten_ds_sg_sw_se$Pi <- gsub("se",   "piN_minus_piS", molten_ds_sg_sw_se$Pi)

## cbind the two molten dfs
molten_ds_sg_sw <- cbind(molten_ds_sg_sw, molten_ds_sg_sw_se)
molten_ds_sg_sw$Pi1 <- NULL
molten_ds_sg_sw$first_site2 <- NULL

## prelim plot
axis_formatting <- theme(axis.text.x = element_text(size = 6),
                         axis.text.y = element_text(size = 6),
                         axis.title.x = element_text(size = 6, margin = margin(t = 4)),
                         axis.title.y = element_text(size = 6, margin = margin(r = 4)))

legend_formatting <- theme(legend.text = element_text(size = 6),
                           legend.title = element_text(size = 8),
                           legend.key.height= unit(0.5, 'cm'),
                           legend.key.width= unit(0.5, 'cm'))

background_formatting <- theme(panel.border = element_rect(color = "grey", fill = NA, size = .5),
                               panel.grid = element_blank(),
                               strip.background = element_blank(),
                               panel.background = element_blank(),
                               legend.background = element_blank())

title <- paste(
  gsub("_", " ", 
       gsub("/recalled/snpgenie","",
            gsub("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/", "", 
                 paste(dir_sg)))), "sliding window of HA", sep = ": ")


palette_muts_NS_S <- c("NS" = "#FF7F20",
                       "S" = "#4F7899")

palette_H3_anti_s <- c("A"   = "red",
                       "B"   = "darkgreen",
                       "C"   = "blue",
                       "D"   = "orange",
                       "E"   = "purple")

palette_H1_anti_s <- c("Cb"  = "red",
                       "Sa"  = "darkgreen",
                       "Ca1" = "blue",
                       "Ca2" = "orange",
                       "Sb"  = "purple")

# assign the appropriate palette
palette_HA_anti_s <- if (grepl("H3N2", title)==T) palette_H3_anti_s else palette_H1_anti_s

# gsub NS and S to simplify legend
molten_ds_sg_sw$Pi <- gsub("Nonsynonymous", "NS", molten_ds_sg_sw$Pi)
molten_ds_sg_sw$Pi <- gsub("Synonymous",     "S", molten_ds_sg_sw$Pi)
molten_ds_sg_sw$Pi <- as.factor(molten_ds_sg_sw$Pi)
molten_ds_sg_sw_pi <- molten_ds_sg_sw[molten_ds_sg_sw$Pi!="piN_minus_piS",]

temp1 <- ggplot() + 
  geom_errorbar(data = molten_ds_sg_sw_pi, 
                aes(x = first_site, y = Mean, ymin = Mean-SE, ymax = Mean+SE, color = Pi, group = Pi), 
                position = position_dodge(), width = .3, alpha = .2) + 
  geom_point(data = molten_ds_sg_sw_pi, 
             aes(x = first_site, y = Mean, color = Pi, group = Pi), 
             size = 1, position = position_dodge(.5), stroke = NA) + 
  coord_cartesian(ylim = c(0, c(max(molten_ds_sg_sw_pi$Mean)+max(molten_ds_sg_sw_pi$SE)))) + 
  scale_color_manual(values = palette_muts_NS_S) +   
  labs(x="Nucleotide postion", y="Nucleotide diversity (\u03c0)", 
       title = title, color = "Mutation\ntype") + 
  axis_formatting + background_formatting + legend_formatting + 
  theme(legend.key = element_blank(),
        legend.position = "right")

temp2 <- ggplot() + 
  geom_rect(data=sg_sw_antigenic[!is.na(sg_sw_antigenic$Antigenic_Site),], inherit.aes=FALSE,
            aes(xmin=Start, xmax=Stop, ymin=-Inf, ymax=Inf,
                group=Antigenic_Site, fill = factor(Antigenic_Site)), alpha=0.2)+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey") + 
  geom_errorbar(data = molten_ds_sg_sw[molten_ds_sg_sw$Pi=="piN_minus_piS",], 
                aes(x = first_site, y = Mean, ymin = Mean-SE, ymax = Mean+SE), 
                position = position_dodge(), width = .3, alpha = .2, color = "black") + 
  geom_point(data = molten_ds_sg_sw[molten_ds_sg_sw$Pi=="piN_minus_piS",], 
             aes(x = first_site, y = Mean), position = position_dodge(.5), 
             size = 1, stroke = NA, color = "black") + 
  scale_fill_manual(values = palette_HA_anti_s) + 
  coord_cartesian(ylim = c(-c(max(molten_ds_sg_sw$Mean)+max(molten_ds_sg_sw$SE)), 
                  c(max(molten_ds_sg_sw$Mean)+max(molten_ds_sg_sw$SE)))) + 
  labs(x="Nucleotide postion", y="\u03c0N - \u03c0S", fill = "Antigenic\nRegion") + 
  axis_formatting + legend_formatting + background_formatting +
  theme(legend.key = element_blank(),
        legend.position = "right")

plot_Pi_sliding_window <- plot_grid(temp1, temp2, 
                                    ncol=1, rel_heights = c(1.1, 1),
                                    align="hv", axis ="lr"); plot_Pi_sliding_window

#### df_sg_AA_HA ####
## change "*" to 0 values for PiN and PiS
sg_sw_1$piN <- as.numeric(gsub("\\*", "0", sg_sw_1$piN))
sg_sw_1$piS <- as.numeric(gsub("\\*", "0", sg_sw_1$piS))
df_sg_AA_HA <- sg_sw_1

## find difference 
df_sg_AA_HA <- transform(df_sg_AA_HA, piNminuspiS=piN-piS)

## subtract UTR from first_site, plus three to count first position as 1, not 0
df_sg_AA_HA <- transform(df_sg_AA_HA, ORF_first_pos=first_site-min(df_sg_AA_HA$first_site)+3)

title <- paste(gsub("_", " ", 
                    gsub("/recalled/snpgenie","",
                         gsub("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/", "", 
                              paste(dir_sg)))))

if (grepl("H3N2", title)==T) {
  df_sg_AA_HA <- transform(df_sg_AA_HA, adj_first_site=first_site-63)
}
if (grepl("H1N1", title)==T) {
  df_sg_AA_HA <- transform(df_sg_AA_HA, adj_first_site=first_site-70)
}

## create AA_pos column denoting the AA position in question
df_sg_AA_HA <- transform(df_sg_AA_HA, AA_pos=adj_first_site/3)

## null out unwanted columns
df_sg_AA_HA$product          <- NULL
df_sg_AA_HA$ORF_first_pos    <- NULL
df_sg_AA_HA$last_site        <- NULL
df_sg_AA_HA$last_codon       <- NULL
df_sg_AA_HA$sum_nonsyn_diffs <- NULL
df_sg_AA_HA$sum_syn_diffs    <- NULL
df_sg_AA_HA$sum_nonsyn_sites <- NULL
df_sg_AA_HA$sum_syn_sites    <- NULL
df_sg_AA_HA$`piN.piS`        <- NULL
df_sg_AA_HA$`piN.piS_vs_ref` <- NULL
df_sg_AA_HA$piN_vs_ref       <- NULL
df_sg_AA_HA$piS_vs_ref       <- NULL
df_sg_AA_HA$adj_first_site   <- NULL

## Apply antigenicity annotations
df_sg_AA_HA <- merge(df_sg_AA_HA, sg_sw_antigenic,
                     by.x = "AA_pos",
                     by.y = "HA_number",
                     all.x = TRUE)
df_sg_AA_HA$Start <- NULL
df_sg_AA_HA$Stop  <- NULL

df_sg_AA_HA$ID[is.na(df_sg_AA_HA$ID)] <- "Non-antigenic"
df_sg_AA_HA$ID <- as.factor(df_sg_AA_HA$ID)
df_sg_AA_HA$piN <- as.numeric(df_sg_AA_HA$piN)
df_sg_AA_HA$piS <- as.numeric(df_sg_AA_HA$piS)
df_sg_AA_HA$piNminuspiS <- as.numeric(df_sg_AA_HA$piNminuspiS)

#### df_AA_HA ####
## Add "1" value in columns for counts
df_sg_AA_HA$One <- 1

### piN-piS
## group by antigenicity (yes/no) and site (if applicable)
df_sg_AA_HA$GROUP <- paste(df_sg_AA_HA$file, paste(df_sg_AA_HA$ID, df_sg_AA_HA$Antigenic_Site, sep = "_"), sep = "/")
df_AA_HA_piN_piS <- aggregate(df_sg_AA_HA$piNminuspiS, by=list(GROUP=df_sg_AA_HA$GROUP), FUN=sum)
# count how many occurance of ONE there are for every GROUP
df_AA_HA_1 <- aggregate(df_sg_AA_HA$One, by=list(GROUP=df_sg_AA_HA$GROUP), FUN=sum)
# name change from default
names(df_AA_HA_1)[names(df_AA_HA_1) == "x"] <- "Count"
names(df_AA_HA_piN_piS)[names(df_AA_HA_piN_piS) == "x"] <- "piN_minus_piS"
# merge the two dfs together
df_AA_HA_1 <- merge(df_AA_HA_piN_piS,
                    df_AA_HA_1,
                    by = "GROUP")
# calculate piN-piS per antigenic/non-ant. site
df_AA_HA_1 <- transform(df_AA_HA_1, perSite_piN_minus_piS=piN_minus_piS/Count)
## group by antigenic (yes/no) only so we can get all of the antigenic values
df_sg_AA_HA$GROUP_2 <- paste(df_sg_AA_HA$file, paste(df_sg_AA_HA$ID, "all", sep = "_"), sep = "/")
df_AA_HA_piN_piS_2 <- aggregate(df_sg_AA_HA$piNminuspiS, by=list(GROUP=df_sg_AA_HA$GROUP_2), FUN=sum)
# count how many occurance of ONE there are for every GROUP
df_AA_HA_1_2 <- aggregate(df_sg_AA_HA$One, by=list(GROUP=df_sg_AA_HA$GROUP_2), FUN=sum)
# name change from default
names(df_AA_HA_1_2)[names(df_AA_HA_1_2) == "x"] <- "Count"
names(df_AA_HA_piN_piS_2)[names(df_AA_HA_piN_piS_2) == "x"] <- "piN_minus_piS"
# merge the two dfs together
df_AA_HA_2 <- merge(df_AA_HA_piN_piS_2,
                    df_AA_HA_1_2,
                    by = "GROUP")
# calculate piN-piS per antigenic/non-ant. site
df_AA_HA_2 <- transform(df_AA_HA_2, perSite_piN_minus_piS=piN_minus_piS/Count)
## rbind df_AA_HA_1 and df_AA_HA_2 together
df_AA_HA_piNminus_piS <- rbind(df_AA_HA_1,
                  df_AA_HA_2,
                  by = "GROUP")
# weird bug where GROUP is a row. Remove it.
df_AA_HA_piNminus_piS <- df_AA_HA_piNminus_piS[df_AA_HA_piNminus_piS$GROUP!="GROUP",]




### piN
## group by antigenicity (yes/no) and site (if applicable)
df_sg_AA_HA$GROUP <- paste(df_sg_AA_HA$file, paste(df_sg_AA_HA$ID, df_sg_AA_HA$Antigenic_Site, sep = "_"), sep = "/")
df_AA_HA_piN <- aggregate(df_sg_AA_HA$piN, by=list(GROUP=df_sg_AA_HA$GROUP), FUN=sum)
# count how many occurance of ONE there are for every GROUP
df_AA_HA_1 <- aggregate(df_sg_AA_HA$One, by=list(GROUP=df_sg_AA_HA$GROUP), FUN=sum)
# name change from default
names(df_AA_HA_1)[names(df_AA_HA_1) == "x"] <- "Count"
names(df_AA_HA_piN)[names(df_AA_HA_piN) == "x"] <- "piN"
# merge the two dfs together
df_AA_HA_1 <- merge(df_AA_HA_piN,
                    df_AA_HA_1,
                    by = "GROUP")
# calculate piN per antigenic/non-ant. site
df_AA_HA_1 <- transform(df_AA_HA_1, perSite_piN=piN/Count)
## group by antigenic (yes/no) only so we can get all of the antigenic values
df_sg_AA_HA$GROUP_2 <- paste(df_sg_AA_HA$file, paste(df_sg_AA_HA$ID, "all", sep = "_"), sep = "/")
df_AA_HA_piN_2 <- aggregate(df_sg_AA_HA$piN, by=list(GROUP=df_sg_AA_HA$GROUP_2), FUN=sum)
# count how many occurance of ONE there are for every GROUP
df_AA_HA_1_2 <- aggregate(df_sg_AA_HA$One, by=list(GROUP=df_sg_AA_HA$GROUP_2), FUN=sum)
# name change from default
names(df_AA_HA_1_2)[names(df_AA_HA_1_2) == "x"] <- "Count"
names(df_AA_HA_piN_2)[names(df_AA_HA_piN_2) == "x"] <- "piN"
# merge the two dfs together
df_AA_HA_2 <- merge(df_AA_HA_piN_2,
                    df_AA_HA_1_2,
                    by = "GROUP")
# calculate piN per antigenic/non-ant. site
df_AA_HA_2 <- transform(df_AA_HA_2, perSite_piN=piN/Count)
## rbind df_AA_HA_1 and df_AA_HA_2 together
df_AA_HA_piN <- rbind(df_AA_HA_1,
                  df_AA_HA_2,
                  by = "GROUP")
# weird bug where GROUP is a row. Remove it.
df_AA_HA_piN <- df_AA_HA_piN[df_AA_HA_piN$GROUP!="GROUP",]




### piS
## group by antigenicity (yes/no) and site (if applicable)
df_sg_AA_HA$GROUP <- paste(df_sg_AA_HA$file, paste(df_sg_AA_HA$ID, df_sg_AA_HA$Antigenic_Site, sep = "_"), sep = "/")
df_AA_HA_piS <- aggregate(df_sg_AA_HA$piS, by=list(GROUP=df_sg_AA_HA$GROUP), FUN=sum)
# count how many occurance of ONE there are for every GROUP
df_AA_HA_1 <- aggregate(df_sg_AA_HA$One, by=list(GROUP=df_sg_AA_HA$GROUP), FUN=sum)
# name change from default
names(df_AA_HA_1)[names(df_AA_HA_1) == "x"] <- "Count"
names(df_AA_HA_piS)[names(df_AA_HA_piS) == "x"] <- "piS"
# merge the two dfs together
df_AA_HA_1 <- merge(df_AA_HA_piS,
                    df_AA_HA_1,
                    by = "GROUP")
# calculate piS per antigenic/non-ant. site
df_AA_HA_1 <- transform(df_AA_HA_1, perSite_piS=piS/Count)
## group by antigenic (yes/no) only so we can get all of the antigenic values
df_sg_AA_HA$GROUP_2 <- paste(df_sg_AA_HA$file, paste(df_sg_AA_HA$ID, "all", sep = "_"), sep = "/")
df_AA_HA_piS_2 <- aggregate(df_sg_AA_HA$piS, by=list(GROUP=df_sg_AA_HA$GROUP_2), FUN=sum)
# count how many occurance of ONE there are for every GROUP
df_AA_HA_1_2 <- aggregate(df_sg_AA_HA$One, by=list(GROUP=df_sg_AA_HA$GROUP_2), FUN=sum)
# name change from default
names(df_AA_HA_1_2)[names(df_AA_HA_1_2) == "x"] <- "Count"
names(df_AA_HA_piS_2)[names(df_AA_HA_piS_2) == "x"] <- "piS"
# merge the two dfs together
df_AA_HA_2 <- merge(df_AA_HA_piS_2,
                    df_AA_HA_1_2,
                    by = "GROUP")
# calculate piS per antigenic/non-ant. site
df_AA_HA_2 <- transform(df_AA_HA_2, perSite_piS=piS/Count)
## rbind df_AA_HA_1 and df_AA_HA_2 together
df_AA_HA_piS <- rbind(df_AA_HA_1,
                  df_AA_HA_2,
                  by = "GROUP")
# weird bug where GROUP is a row. Remove it.
df_AA_HA_piS <- df_AA_HA_piS[df_AA_HA_piS$GROUP!="GROUP",]








df_AA_HA_piN$Count <- NULL
df_AA_HA_piS$Count <- NULL
df_AA_HA_piNminus_piS$Count <- NULL
df_AA_HA <- merge(df_AA_HA_piN,
                  df_AA_HA_piS,
                  by = "GROUP")

df_AA_HA <- merge(df_AA_HA,
                  df_AA_HA_piNminus_piS,
                  by = "GROUP")

# separate GROUP column (smart, eh!?)
df_AA_HA <- separate(df_AA_HA,"GROUP",
                         c("file", "ID"),
                         sep="/",convert=FALSE)
df_AA_HA <- separate(df_AA_HA,"ID",
                         c("Antigenicity", "Antigenic_Site"),
                         sep="_",convert=FALSE)
# remove duplicate Antigenic_NA rows
df_AA_HA <- df_AA_HA[df_AA_HA$Antigenic_Site!="NA",]


#### df_gene_counts ####
## Add "1" value in columns for counts
df_all_iSNVs$One <- 1

## Calculate per-sample-per-gene counts
# count how many occurance of ONE there are for every GROUP
df_gene_counts <- aggregate(df_all_iSNVs$One, by=list(GROUP=df_all_iSNVs$GROUP), FUN=sum)
# separate GROUP column (smart, eh!?)
df_gene_counts <- separate(df_gene_counts,"GROUP",
                           c("WSLH_ID", "Gene", "Antigenicity", "Antigenic_site", "RBD"),
                           sep="_",convert=FALSE)
# name change from default
names(df_gene_counts)[names(df_gene_counts) == "x"] <- "Count"
# factor
df_gene_counts$WSLH_ID <- as.factor(df_gene_counts$WSLH_ID)
df_gene_counts$Gene <- as.factor(df_gene_counts$Gene)
df_gene_counts$Antigenicity <- as.factor(df_gene_counts$Antigenicity)

## summary
summary(df_gene_counts)

#### df_sub50AF ####
## filter to only include iSNVs
df_sub50AF <- filter(df_all_iSNVs, AF<=.5)

## only NS, S, and Stop
df_sub50AF$Ann <- as.factor(as.character(df_sub50AF$Ann))
levels(df_sub50AF$Ann)
df_sub50AF_NS <- filter(df_sub50AF, Ann=="missense_variant")
df_sub50AF_S <- filter(df_sub50AF, Ann=="synonymous_variant")
df_sub50AF_Stop <- filter(df_sub50AF, Ann=="stop_gained")
df_sub50AF <- rbind(df_sub50AF_NS, df_sub50AF_S, df_sub50AF_Stop)
df_sub50AF$Ann <- gsub("missense_variant","Nonsynonymous",df_sub50AF$Ann)
df_sub50AF$Ann <- gsub("synonymous_variant","Synonymous",df_sub50AF$Ann)
df_sub50AF$Ann <- gsub("stop_gained","Stop gained",df_sub50AF$Ann)

## add AF bins
df_sub50AF$Bin[df_sub50AF$AF>minAF & df_sub50AF$AF<=0.1] <- "3-10%"
df_sub50AF$Bin[df_sub50AF$AF>0.1 & df_sub50AF$AF<=0.2] <- "10-20%"
df_sub50AF$Bin[df_sub50AF$AF>0.2 & df_sub50AF$AF<=0.3] <- "20-30%"
df_sub50AF$Bin[df_sub50AF$AF>0.3 & df_sub50AF$AF<=0.4] <- "30-40%"
df_sub50AF$Bin[df_sub50AF$AF>0.4 & df_sub50AF$AF<=0.5] <- "40-50%"

## freq to percent
df_sub50AF$AF <- df_sub50AF$AF*100

#### df_sub50AF_md ####
## merge with md
df_sub50AF_md <- merge(df_sub50AF, md,
                       all.x = TRUE,
                       by.x = "WSLH_ID",
                       by.y = "WSLH_ID")
df_sub50AF_md$Subtype <- as.factor(df_sub50AF_md$Subtype)
df_sub50AF_md$flu_vaccine <- as.factor(df_sub50AF_md$flu_vaccine)
df_sub50AF_md$WSLH_ID <- as.factor(as.character(df_sub50AF_md$WSLH_ID))

#### df_sub50AF_md_n_SNVs_by_subtype ####
## count subtypes by vax status
df_sub50AF_md_n_SNVs_by_subtype <- as.data.frame(table(df_sub50AF_md$Subtype,df_sub50AF_md$WSLH_ID))
names(df_sub50AF_md_n_SNVs_by_subtype)[names(df_sub50AF_md_n_SNVs_by_subtype) == "Var1"] <- "Subtype"
names(df_sub50AF_md_n_SNVs_by_subtype)[names(df_sub50AF_md_n_SNVs_by_subtype) == "Var2"] <- "WSLH_ID"
names(df_sub50AF_md_n_SNVs_by_subtype)[names(df_sub50AF_md_n_SNVs_by_subtype) == "Freq"] <- "Count"

## factors
df_sub50AF_md_n_SNVs_by_subtype$Subtype <- factor(df_sub50AF_md_n_SNVs_by_subtype$Subtype,
                                                  levels = c("H3N2", "H1N1", "B"))

#### df_sub50AF_md_n_SNVs_by_vax ####
df_sub50AF_md_n_SNVs_by_vax <- as.data.frame(table(df_sub50AF_md$flu_vaccine,df_sub50AF_md$WSLH_ID))
names(df_sub50AF_md_n_SNVs_by_vax)[names(df_sub50AF_md_n_SNVs_by_vax) == "Var1"] <- "flu_vaccine"
names(df_sub50AF_md_n_SNVs_by_vax)[names(df_sub50AF_md_n_SNVs_by_vax) == "Var2"] <- "WSLH_ID"
names(df_sub50AF_md_n_SNVs_by_vax)[names(df_sub50AF_md_n_SNVs_by_vax) == "Freq"] <- "Count"
df_sub50AF_md_n_SNVs_by_vax <- df_sub50AF_md_n_SNVs_by_vax[df_sub50AF_md_n_SNVs_by_vax$Count>0,]

## factors
df_sub50AF_md_n_SNVs_by_vax$Subtype <- factor(df_sub50AF_md_n_SNVs_by_vax$flu_vaccine,
                                              levels = c("Vaccinated", "Unvaccinated"))

#### df_sub50AF_md_n_SNVs_by_bin_and_mut ####
## create df_sub50AF_md_n_SNVs_by_bin_and_mut
df_sub50AF_md_n_SNVs_by_bin_and_mut <- as.data.frame(table(df_sub50AF$Bin, df_sub50AF$Ann, df_sub50AF$WSLH_ID))
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var1"] <- "Bins"
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var2"] <- "Mutation_type"
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Var3"] <- "WSLH_ID"
names(df_sub50AF_md_n_SNVs_by_bin_and_mut)[names(df_sub50AF_md_n_SNVs_by_bin_and_mut) == "Freq"] <- "Count"
df_sub50AF_md_n_SNVs_by_bin_and_mut$ID <- paste(df_sub50AF_md_n_SNVs_by_bin_and_mut$Bins,
                                                df_sub50AF_md_n_SNVs_by_bin_and_mut$Mutation_type,
                                                df_sub50AF_md_n_SNVs_by_bin_and_mut$WSLH_ID, sep = "/")

#### list_mut_bins_prop ####
### df_sub50AF_md_n_SNVs_by_bin_and_mut
## split df by WSLH_ID
list_mut_bins_prop <- split(df_sub50AF_md_n_SNVs_by_bin_and_mut, 
                            with(df_sub50AF_md_n_SNVs_by_bin_and_mut, 
                                 WSLH_ID, drop=T))
## split Mutation_type lists by Mutation_type
for (i in 1:length(list_mut_bins_prop)) {
  list_mut_bins_prop[[i]] <- split(list_mut_bins_prop[[i]], 
                                   with(list_mut_bins_prop[[i]], 
                                        Mutation_type, drop=T))
}
## calculate total number of iSNVs per Mutation_type per WSLH_ID, the proportional column
# run the calculation
for (i in 1:length(list_mut_bins_prop)) {
  # will loop through every WSLH_ID
  for (j in 1:length(list_mut_bins_prop[[1]])) {
    # will loop through every Mutation_type
    for (k in 1:length(list_mut_bins_prop[[i]][[j]]$Bins)) {
      # will loop through every bin (some samples don't have iSNVs in every bin)
      list_mut_bins_prop[[i]][[j]]$Prop[k] <- sum(list_mut_bins_prop[[i]][[j]][k,"Count"]/sum(list_mut_bins_prop[[i]][[j]][,"Count"]))
    }
  }
}

## reduce from twice-nested list to once-nested list
for (i in 1:length(list_mut_bins_prop)) {
  list_mut_bins_prop[[i]] <- Reduce(full_join, list_mut_bins_prop[[i]])
}

#### df_mut_bins_prop ####
## list to df
df_mut_bins_prop_all <- Reduce(full_join,list_mut_bins_prop)

## only include some of the mutation types
df_mut_bins_prop_S <- filter(df_mut_bins_prop_all, Mutation_type=="Synonymous")
df_mut_bins_prop_NS <- filter(df_mut_bins_prop_all, Mutation_type=="Nonsynonymous")
df_mut_bins_prop_Stop_g <- filter(df_mut_bins_prop_all, Mutation_type=="Stop gained")
df_mut_bins_prop <- rbind(df_mut_bins_prop_S, df_mut_bins_prop_NS, df_mut_bins_prop_Stop_g)

## convert NaN values (0/0) to zero
df_mut_bins_prop$Prop[df_mut_bins_prop$Prop=="NaN"] <- 0

## generate proportion mean and 95% confidence intervals for each bin and mutation_type
# split by Mutation_type_bin
df_mut_bins_prop$Mutation_type_bin <- paste(df_mut_bins_prop$Mutation_type, df_mut_bins_prop$Bins, sep = "_") 
list_mut_bins_prop_summary <- split(df_mut_bins_prop, with(df_mut_bins_prop, Mutation_type_bin, drop=T))
for (i in 1:length(list_mut_bins_prop_summary)) {
  # where i is a level of a list (e.g., `Synonymous_40-50%`)
  list_mut_bins_prop_summary[[i]] <- data.frame("Bins"          = list_mut_bins_prop_summary[[i]]$Bins[1],
                                                "Mutation_type" = list_mut_bins_prop_summary[[i]]$Mutation_type[1],
                                                "Mean_prop"     = mean(list_mut_bins_prop_summary[[i]]$Prop),
                                                "StdError_prop" = se(list_mut_bins_prop_summary[[i]]$Prop),
                                                "Total"         = sum(list_mut_bins_prop_summary[[i]]$Count))
  }

## reduce back to df_mut_bins_prop
df_mut_bins_prop <- Reduce(full_join, list_mut_bins_prop_summary)

## save a copy for the no-neutral df
df_mut_bins_prop_noNeut <- df_mut_bins_prop

## df for neutral expectation
df_mut_bins_neutral$...1 <- NULL
colnames(df_mut_bins_neutral) <- c("Bins", "Mean_prop")
df_mut_bins_neutral$StdError_prop <- 0
df_mut_bins_neutral$Mutation_type <- "Neutral expectation"
df_mut_bins_neutral$Total <- ""

## factor order
df_mut_bins_prop_noNeut$Mutation_type <- factor(df_mut_bins_prop_noNeut$Mutation_type,
                                                levels = c("Nonsynonymous", "Synonymous", "Stop gained"))

## add neutral expectation
df_mut_bins_prop <- rbind(df_mut_bins_prop, df_mut_bins_neutral)

## factor order
df_mut_bins_prop$Mutation_type <- factor(df_mut_bins_prop$Mutation_type,
                                         levels = c("Nonsynonymous", "Synonymous",
                                                    "Stop gained","Neutral expectation"))

## preliminary plot
# palette
palette_muts <- c("Nonsynonymous" = "#FF7F20",
                  "Synonymous" = "#4F7899",
                  "Stop gained" = "#7AD9C2")
# plot
plot_mut_bins_prop <- ggplot() + 
  geom_bar(data = df_mut_bins_prop_noNeut, 
           aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")), 
               y = Mean_prop, fill = Mutation_type),
           stat = "identity", position = position_dodge()) +
  geom_text(data = df_mut_bins_prop_noNeut, 
            aes(x = Bins, y = 0, group = Mutation_type, label = Total),
            position=position_dodge(width=0.9), vjust=1.5, size = 1.5) + 
  geom_errorbar(data = df_mut_bins_prop_noNeut, 
                aes(x = Bins, y = Mean_prop, 
                    ymin = Mean_prop-StdError_prop, ymax = Mean_prop+StdError_prop, 
                    group = Mutation_type, width = .3), color = "grey", position=position_dodge(.9)) +   
  scale_x_discrete(breaks = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")) + 
  geom_point(data = df_mut_bins_neutral, 
             aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")),
                 y = Mean_prop), size = 1) + 
  geom_line(data = df_mut_bins_neutral, 
            aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")), 
                y = Mean_prop, group = 1), linewidth = .5) +
  scale_fill_manual(values = palette_muts) + 
  labs(x = "Within-host iSNV frequency bin", 
       y = "Proportion of variants per mutation type",
       title = title) + 
  ylim(0,1) + 
  axis_formatting + 
  legend_formatting + 
  theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = c(.8, .8)) + 
  background_formatting

#### df_divergence_antigenic_Ann ####
## calculate divergence by totaling AF for each sample
df_divergence_antigenic_Ann <- aggregate(df_all_iSNVs$AF, by=list(GROUP=df_all_iSNVs$GROUP, df_all_iSNVs$Ann), FUN=sum)
df_divergence_antigenic_Ann <- separate(df_divergence_antigenic_Ann,"GROUP", c("WSLH_ID","Gene_Name", "Antigenicity", "Antigenic_Site"), sep="_",convert=FALSE)

## remove "Non-antigenic" from non-HA genes
df_divergence_antigenic_Ann$Antigenicity[df_divergence_antigenic_Ann$Gene_Name!="HA"] <- ""

## change names
names(df_divergence_antigenic_Ann)[names(df_divergence_antigenic_Ann) == "x"] <- "Divergence"
names(df_divergence_antigenic_Ann)[names(df_divergence_antigenic_Ann) == "Group.2"] <- "Ann"

## remove Antigenic_Site cols
df_divergence_antigenic_Ann$Antigenic_Site <- NULL

## add GROUP column
df_divergence_antigenic_Ann$GROUP <- paste(df_divergence_antigenic_Ann$WSLH_ID, df_divergence_antigenic_Ann$Gene_Name)
df_divergence_antigenic_Ann$GROUP <- gsub(" ", "-", df_divergence_antigenic_Ann$GROUP)

## do it again for just HA in general (combining antigenic and nonantigenic)
df_divergence_antigenic_Ann_HA <- aggregate(df_divergence_antigenic_Ann$Divergence, 
                                            by=list(GROUP=df_divergence_antigenic_Ann$GROUP, 
                                                    df_divergence_antigenic_Ann$Ann), FUN=sum)
df_divergence_antigenic_Ann_HA <- separate(df_divergence_antigenic_Ann_HA,"GROUP", c("WSLH_ID","Gene_Name"), sep="-",convert=FALSE)
names(df_divergence_antigenic_Ann_HA)[names(df_divergence_antigenic_Ann_HA) == "x"] <- "Divergence"
names(df_divergence_antigenic_Ann_HA)[names(df_divergence_antigenic_Ann_HA) == "Group.2"] <- "Ann"
df_divergence_antigenic_Ann_HA <- filter(df_divergence_antigenic_Ann_HA, Gene_Name=="HA")
df_divergence_antigenic_Ann_HA$GROUP <- paste(df_divergence_antigenic_Ann_HA$WSLH_ID, df_divergence_antigenic_Ann_HA$Gene_Name)
df_divergence_antigenic_Ann_HA$GROUP <- gsub(" ", "-", df_divergence_antigenic_Ann_HA$GROUP)
df_divergence_antigenic_Ann_HA$Antigenicity <- "Combined"

## add compiled HA to split HA
df_divergence_antigenic_Ann <- rbind(df_divergence_antigenic_Ann_HA, df_divergence_antigenic_Ann)
df_divergence_antigenic_Ann_HA <- NULL

## factor
df_divergence_antigenic_Ann$WSLH_ID <- as.factor(df_divergence_antigenic_Ann$WSLH_ID)
df_divergence_antigenic_Ann$Gene_Name <- as.factor(df_divergence_antigenic_Ann$Gene_Name)
df_divergence_antigenic_Ann$Antigenicity <- as.factor(df_divergence_antigenic_Ann$Antigenicity)

## just NS and S
levels(as.factor(as.character(df_divergence_antigenic_Ann$Ann)))
df_divergence_NS <- filter(df_divergence_antigenic_Ann, Ann=="missense_variant")
df_divergence_S <- filter(df_divergence_antigenic_Ann, Ann=="synonymous_variant")
df_divergence_antigenic_Ann <- rbind(df_divergence_NS, df_divergence_S)

## rename NS and S
df_divergence_antigenic_Ann$Ann <- gsub("missense_variant","Nonsynonymous",df_divergence_antigenic_Ann$Ann)
df_divergence_antigenic_Ann$Ann <- gsub("synonymous_variant","Synonymous",df_divergence_antigenic_Ann$Ann)
df_divergence_antigenic_Ann$Ann <- as.factor(df_divergence_antigenic_Ann$Ann)

## Create x axis factors
df_divergence_antigenic_Ann$x_axis <- paste(df_divergence_antigenic_Ann$Gene_Name,df_divergence_antigenic_Ann$Antigenicity)
df_divergence_antigenic_Ann$x_axis <- gsub(" NA", "", df_divergence_antigenic_Ann$x_axis)
df_divergence_antigenic_Ann$x_axis <- gsub(" ", "-", df_divergence_antigenic_Ann$x_axis)

## merge df with sg
df_divergence_antigenic_Ann <- merge(df_divergence_antigenic_Ann, sg, 
                                     by.x = "GROUP", by.y = "GROUP") ## sg does not have same GROUP as df

## two dfs (back together soon!)
df_divergence_S <- filter(df_divergence_antigenic_Ann, Ann == "Synonymous")
df_divergence_NS <- filter(df_divergence_antigenic_Ann, Ann == "Nonsynonymous")

## divergence per site
df_divergence_S <- transform(df_divergence_S, DivergencePerSite = df_divergence_S$Divergence / df_divergence_S$S_sites)
df_divergence_NS <- transform(df_divergence_NS, DivergencePerSite = df_divergence_NS$Divergence / df_divergence_NS$N_sites)

## rbind
df_divergence <- rbind(df_divergence_S, df_divergence_NS)

## data summary
ds_divergence <- as.data.frame(ds(df_divergence, 
                                  varname="DivergencePerSite", 
                                  groupnames=c("x_axis", "Ann")))

## factor x axis
levels(as.factor(as.character(df_divergence$x_axis)))
df_divergence$x_axis <- gsub("[[:punct:]]$", "", df_divergence$x_axis)
df_divergence$x_axis <- factor(df_divergence$x_axis,
                               levels = c("PB2", "PB1", "PA", "HA-Combined", 
                                          "HA-Antigenic", "HA-Non-antigenic",
                                          "NP", "NA", "M1", "M2", 
                                          "NS1", "NEP"))
ds_divergence$x_axis <- gsub("[[:punct:]]$", "", ds_divergence$x_axis)
ds_divergence$x_axis <- factor(ds_divergence$x_axis,
                               levels = c("PB2", "PB1", "PA", "HA-Combined", 
                                          "HA-Antigenic", "HA-Non-antigenic",
                                          "NP", "NA", "M1", "M2", 
                                          "NS1", "NEP"))

#### df_AA_HA_paired ####
## merge with md
df_AA_HA_paired <- merge(df_AA_HA, md, by.x="file", by.y="WSLH_ID")

## factor
df_AA_HA_paired$Case <- as.factor(df_AA_HA_paired$Case)
df_AA_HA_paired$Case_Num <- as.factor(df_AA_HA_paired$Case_Num)

## find samples with pairs by pulling by group name and excluding other group, reciprocally
# recips
string_df_AA_HA_recips <- as.character(gsub("Recipient ", "", levels(as.factor(df_AA_HA_paired$Case_Num[df_AA_HA_paired$Case=="Recipient"]))))
string_df_AA_HA_recips <- as.numeric(gsub("Donor .*", "0", string_df_AA_HA_recips))
string_df_AA_HA_recips <- string_df_AA_HA_recips[string_df_AA_HA_recips>0]
# donors
string_df_AA_HA_donors <- as.character(gsub("Donor ", "", levels(as.factor(df_AA_HA_paired$Case_Num[df_AA_HA_paired$Case=="Donor"]))))
string_df_AA_HA_donors <- as.numeric(gsub("Recipient .*", "0", string_df_AA_HA_donors))
string_df_AA_HA_donors <- string_df_AA_HA_donors[string_df_AA_HA_donors>0]
# intersection
string_df_AA_HA_pairs <- intersect(string_df_AA_HA_recips, string_df_AA_HA_donors); string_df_AA_HA_pairs

## remove Case 1544, which has the same WSLH_ID as two unique 
# demographic observations, making it look like a donor-recipient pair
string_df_AA_HA_pairs <- string_df_AA_HA_pairs[string_df_AA_HA_pairs %!in% c("1544")]

## only df_AA_HA_paired data that are pairs
df_AA_HA_paired <- filter(df_AA_HA_paired, Num %in% string_df_AA_HA_pairs)

## factor
df_AA_HA_paired$file <- as.factor(as.character(df_AA_HA_paired$file))
df_AA_HA_paired$Case <- as.factor(as.character(df_AA_HA_paired$Case))
df_AA_HA_paired$Case_Num <- as.factor(as.character(df_AA_HA_paired$Case_Num))

## create identifier and list of df_AA_HA_paired by Identifier
df_AA_HA_paired$Identifier <- paste(df_AA_HA_paired$Num, df_AA_HA_paired$Case, sep="-")
df_AA_HA_list_paired <- split(df_AA_HA_paired, with(df_AA_HA_paired, Identifier, drop=T))

## clean up the list!
n <- length(df_AA_HA_list_paired)
for (i in 1:n) {
  df_AA_HA_list_paired[[i]]$WSLH_ID <- as.factor(as.character(df_AA_HA_list_paired[[i]]$file))
  df_AA_HA_list_paired[[i]]$ID <- as.integer(df_AA_HA_list_paired[[i]]$WSLH_ID)
  df_AA_HA_list_paired[[i]]$ID[df_AA_HA_list_paired[[i]]$Case=="Donor"] <- 0
  df_AA_HA_list_paired[[i]]$ID <- paste(df_AA_HA_list_paired[[i]]$Case, df_AA_HA_list_paired[[i]]$ID, sep=" ")
}

## back to df
df_AA_HA_paired <- Reduce(full_join, df_AA_HA_list_paired)

## Gene column
df_AA_HA_paired$Gene <- as.factor(paste(df_AA_HA_paired$Antigenicity, df_AA_HA_paired$Antigenic_Site, sep = "_"))

## all_genes
all_genes <- data.frame("Gene"=levels(df_AA_HA_paired$Gene), 
                        "piN"=0, "perSite_piN"=0, 
                        "piS"=0, "perSite_piS"=0, 
                        "piN_minus_piS"=0, "perSite_piN_minus_piS"=0)

## back to list
df_AA_HA_list_paired <- split(df_AA_HA_paired, with(df_AA_HA_paired, WSLH_ID, drop=T))

## add Gene where piN_minus_piS=0 to df_AA_HA_list_paired
n <- length(df_AA_HA_list_paired)
for (i in 1:n) {
  temp <- all_genes
  temp1 <- merge(temp, df_AA_HA_list_paired[[i]][,c("Gene", 
                                                    "piN", "perSite_piN", 
                                                    "piS", "perSite_piS", 
                                                    "piN_minus_piS", "perSite_piN_minus_piS")], by="Gene", all.x=TRUE)
  temp1[is.na(temp1)] <- 0
  
  temp1$piN <- temp1$piN.y
  temp1$piN.x <- NULL
  temp1$piN.y <- NULL
  temp1$perSite_piN <- temp1$perSite_piN.y
  temp1$perSite_piN.x <- NULL
  temp1$perSite_piN.y <- NULL
  temp1$piS <- temp1$piS.y
  temp1$piS.x <- NULL
  temp1$piS.y <- NULL
  temp1$perSite_piS <- temp1$perSite_piS.y
  temp1$perSite_piS.x <- NULL
  temp1$perSite_piS.y <- NULL
  temp1$piN_minus_piS <- temp1$piN_minus_piS.y
  temp1$piN_minus_piS.x <- NULL
  temp1$piN_minus_piS.y <- NULL
  temp1$perSite_piN_minus_piS <- temp1$perSite_piN_minus_piS.y
  temp1$perSite_piN_minus_piS.x <- NULL
  temp1$perSite_piN_minus_piS.y <- NULL
  
  temp1$WSLH_ID <- df_AA_HA_list_paired[[i]]$WSLH_ID[1]
  temp1$Case <- df_AA_HA_list_paired[[i]]$Case[1]
  temp1$Num <- df_AA_HA_list_paired[[i]]$Num[1]
  temp1$Case_Num <- paste(temp1$Case, temp1$Num, sep = "-")
  temp1$ID <- df_AA_HA_list_paired[[i]]$ID[1]
  df_AA_HA_list_paired[[i]] <- temp1
}
temp <- NULL
temp1 <- NULL

## back to df
df_AA_HA_paired <- Reduce(full_join, df_AA_HA_list_paired)

## back to list
df_AA_HA_list_paired <- split(df_AA_HA_paired, with(df_AA_HA_paired, Num, drop=T))

## split the list dfs into twice-nested dfs within list
n <- length(df_AA_HA_list_paired)
for (i in 1:n) {
  df_AA_HA_list_paired[[i]] <- split(df_AA_HA_list_paired[[i]], with(df_AA_HA_list_paired[[i]], ID, drop=T))
}

## only one WSLH_ID per donor or recipient in pairing
n <- length(df_AA_HA_list_paired)
for (i in 1:n) {
  ## donor dups
  df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID <- as.factor(as.character(df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID))
  if (length(levels(df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID)) >= 2) {
    levels <- levels(df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID)
    df_AA_HA_list_paired[[i]]$`Donor 0` <- df_AA_HA_list_paired[[i]]$`Donor 0`[df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID==levels[1],]
  }
  ## recip 1 dups
  df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID <- as.factor(as.character(df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID))
  if (length(levels(df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID)) >= 2) {
    levels <- levels(df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID)
    df_AA_HA_list_paired[[i]]$`Recipient 1` <- df_AA_HA_list_paired[[i]]$`Recipient 1`[df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID==levels[1],]
  }
  ## recip 2 dups
  df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID <- as.factor(as.character(df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID))
  if (length(levels(df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID)) >= 2) {
    levels <- levels(df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID)
    df_AA_HA_list_paired[[i]]$`Recipient 2` <- df_AA_HA_list_paired[[i]]$`Recipient 2`[df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID==levels[1],]
  }
  ## recip 3 dups
  df_AA_HA_list_paired[[i]]$`Recipient 3`$WSLH_ID <- as.factor(as.character(df_AA_HA_list_paired[[i]]$`Recipient 3`$WSLH_ID))
  if (length(levels(df_AA_HA_list_paired[[i]]$`Recipient 3`$WSLH_ID)) >= 2) {
    levels <- levels(df_AA_HA_list_paired[[i]]$`Recipient 3`$WSLH_ID)
    df_AA_HA_list_paired[[i]]$`Recipient 3` <- df_AA_HA_list_paired[[i]]$`Recipient 3`[df_AA_HA_list_paired[[i]]$`Recipient 3`$WSLH_ID==levels[1],]
  }
}

## remove empty dfs and single donors/recips from nested list
# empty
n <- length(df_AA_HA_list_paired)
for (i in 1:n) {
  df_AA_HA_list_paired[[i]] <- Filter(function(x) dim(x)[1] > 0, df_AA_HA_list_paired[[i]])
}
# donor/recips-only by identifying which list of lists is problematic
to_remove <- c(NA)
n <- length(df_AA_HA_list_paired)
for (i in 1:n) {
  if (length(df_AA_HA_list_paired[[i]])==1) {
    to_remove <- c(i, to_remove)
  }
}
# removes identified donor/recip-only cases
if (length(to_remove[!is.na(to_remove)])>0) {
  df_AA_HA_list_paired <- df_AA_HA_list_paired[-to_remove[!is.na(to_remove)]]
}



## create df in pairing format
n <- length(df_AA_HA_list_paired)
df_AA_HA_paired <- data.frame("Gene"=levels(as.factor(df_AA_HA_paired$Gene)), 
                              "Num" = "Remove",
                              "Donor_WSLH_ID" = 0,
                              "Recipient_WSLH_ID" = 0,
                              "Donor_piN_minus_piS" = 0,
                              "Recipient_piN_minus_piS" = 0,
                              "Donor_perSite_piN_minus_piS" = 0,
                              "Recipient_perSite_piN_minus_piS" = 0,
                              "Donor_piN" = 0,
                              "Recipient_piN" = 0,
                              "Donor_perSite_piN" = 0,
                              "Recipient_perSite_piN" = 0,
                              "Donor_piS" = 0,
                              "Recipient_piS" = 0,
                              "Donor_perSite_piS" = 0,
                              "Recipient_perSite_piS" = 0)
for (i in 1:n) {
  if (length(df_AA_HA_list_paired[[i]]) == 2) {
    temp <- data.frame("Num"                                 = df_AA_HA_list_paired[[i]]$`Donor 0`$Num,
                       "Gene"                                = df_AA_HA_list_paired[[i]]$`Donor 0`$Gene,
                       "Donor_WSLH_ID"                       = df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID,
                       "Recipient_WSLH_ID"                   = df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID,
                       "Donor_piN_minus_piS"                 = df_AA_HA_list_paired[[i]]$`Donor 0`$piN_minus_piS,
                       "Recipient_piN_minus_piS"             = df_AA_HA_list_paired[[i]]$`Recipient 1`$piN_minus_piS,
                       "Donor_perSite_piN_minus_piS"         = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN_minus_piS,
                       "Recipient_perSite_piN_minus_piS"     = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piN_minus_piS,
                       "Donor_piN"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piN,
                       "Recipient_piN"                       = df_AA_HA_list_paired[[i]]$`Recipient 1`$piN,
                       "Donor_perSite_piN"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN,
                       "Recipient_perSite_piN"               = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piN,
                       "Donor_piS"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piS,
                       "Recipient_piS"                       = df_AA_HA_list_paired[[i]]$`Recipient 1`$piS,
                       "Donor_perSite_piS"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piS,
                       "Recipient_perSite_piS"               = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piS)
  }
  if (length(df_AA_HA_list_paired[[i]]) == 3) {
    temp1 <- data.frame("Num"                                 = df_AA_HA_list_paired[[i]]$`Donor 0`$Num,
                        "Gene"                                = df_AA_HA_list_paired[[i]]$`Donor 0`$Gene,
                        "Donor_WSLH_ID"                       = df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID,
                        "Recipient_WSLH_ID"                   = df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID,
                        "Donor_piN_minus_piS"                 = df_AA_HA_list_paired[[i]]$`Donor 0`$piN_minus_piS,
                        "Recipient_piN_minus_piS"             = df_AA_HA_list_paired[[i]]$`Recipient 1`$piN_minus_piS,
                        "Donor_perSite_piN_minus_piS"         = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN_minus_piS,
                        "Recipient_perSite_piN_minus_piS"     = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piN_minus_piS,
                        "Donor_piN"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piN,
                        "Recipient_piN"                       = df_AA_HA_list_paired[[i]]$`Recipient 1`$piN,
                        "Donor_perSite_piN"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN,
                        "Recipient_perSite_piN"               = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piN,
                        "Donor_piS"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piS,
                        "Recipient_piS"                       = df_AA_HA_list_paired[[i]]$`Recipient 1`$piS,
                        "Donor_perSite_piS"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piS,
                        "Recipient_perSite_piS"               = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piS)
    temp2 <- data.frame("Num"                                 = df_AA_HA_list_paired[[i]]$`Donor 0`$Num,
                        "Gene"                                = df_AA_HA_list_paired[[i]]$`Donor 0`$Gene,
                        "Donor_WSLH_ID"                       = df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID,
                        "Recipient_WSLH_ID"                   = df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID,
                        "Donor_piN_minus_piS"                 = df_AA_HA_list_paired[[i]]$`Donor 0`$piN_minus_piS,
                        "Recipient_piN_minus_piS"             = df_AA_HA_list_paired[[i]]$`Recipient 2`$piN_minus_piS,
                        "Donor_perSite_piN_minus_piS"         = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN_minus_piS,
                        "Recipient_perSite_piN_minus_piS"     = df_AA_HA_list_paired[[i]]$`Recipient 2`$perSite_piN_minus_piS,
                        "Donor_piN"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piN,
                        "Recipient_piN"                       = df_AA_HA_list_paired[[i]]$`Recipient 2`$piN,
                        "Donor_perSite_piN"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN,
                        "Recipient_perSite_piN"               = df_AA_HA_list_paired[[i]]$`Recipient 2`$perSite_piN,
                        "Donor_piS"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piS,
                        "Recipient_piS"                       = df_AA_HA_list_paired[[i]]$`Recipient 2`$piS,
                        "Donor_perSite_piS"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piS,
                        "Recipient_perSite_piS"               = df_AA_HA_list_paired[[i]]$`Recipient 2`$perSite_piS)
    temp <- rbind(temp1, temp2)
  }
  if (length(df_AA_HA_list_paired[[i]]) == 4) {
    temp1 <- data.frame("Num"                                 = df_AA_HA_list_paired[[i]]$`Donor 0`$Num,
                        "Gene"                                = df_AA_HA_list_paired[[i]]$`Donor 0`$Gene,
                        "Donor_WSLH_ID"                       = df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID,
                        "Recipient_WSLH_ID"                   = df_AA_HA_list_paired[[i]]$`Recipient 1`$WSLH_ID,
                        "Donor_piN_minus_piS"                 = df_AA_HA_list_paired[[i]]$`Donor 0`$piN_minus_piS,
                        "Recipient_piN_minus_piS"             = df_AA_HA_list_paired[[i]]$`Recipient 1`$piN_minus_piS,
                        "Donor_perSite_piN_minus_piS"         = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN_minus_piS,
                        "Recipient_perSite_piN_minus_piS"     = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piN_minus_piS,
                        "Donor_piN"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piN,
                        "Recipient_piN"                       = df_AA_HA_list_paired[[i]]$`Recipient 1`$piN,
                        "Donor_perSite_piN"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN,
                        "Recipient_perSite_piN"               = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piN,
                        "Donor_piS"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piS,
                        "Recipient_piS"                       = df_AA_HA_list_paired[[i]]$`Recipient 1`$piS,
                        "Donor_perSite_piS"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piS,
                        "Recipient_perSite_piS"               = df_AA_HA_list_paired[[i]]$`Recipient 1`$perSite_piS)
    temp2 <- data.frame("Num"                                 = df_AA_HA_list_paired[[i]]$`Donor 0`$Num,
                        "Gene"                                = df_AA_HA_list_paired[[i]]$`Donor 0`$Gene,
                        "Donor_WSLH_ID"                       = df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID,
                        "Recipient_WSLH_ID"                   = df_AA_HA_list_paired[[i]]$`Recipient 2`$WSLH_ID,
                        "Donor_piN_minus_piS"                 = df_AA_HA_list_paired[[i]]$`Donor 0`$piN_minus_piS,
                        "Recipient_piN_minus_piS"             = df_AA_HA_list_paired[[i]]$`Recipient 2`$piN_minus_piS,
                        "Donor_perSite_piN_minus_piS"         = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN_minus_piS,
                        "Recipient_perSite_piN_minus_piS"     = df_AA_HA_list_paired[[i]]$`Recipient 2`$perSite_piN_minus_piS,
                        "Donor_piN"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piN,
                        "Recipient_piN"                       = df_AA_HA_list_paired[[i]]$`Recipient 2`$piN,
                        "Donor_perSite_piN"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN,
                        "Recipient_perSite_piN"               = df_AA_HA_list_paired[[i]]$`Recipient 2`$perSite_piN,
                        "Donor_piS"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piS,
                        "Recipient_piS"                       = df_AA_HA_list_paired[[i]]$`Recipient 2`$piS,
                        "Donor_perSite_piS"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piS,
                        "Recipient_perSite_piS"               = df_AA_HA_list_paired[[i]]$`Recipient 2`$perSite_piS)
    temp3 <- data.frame("Num"                                 = df_AA_HA_list_paired[[i]]$`Donor 0`$Num,
                        "Gene"                                = df_AA_HA_list_paired[[i]]$`Donor 0`$Gene,
                        "Donor_WSLH_ID"                       = df_AA_HA_list_paired[[i]]$`Donor 0`$WSLH_ID,
                        "Recipient_WSLH_ID"                   = df_AA_HA_list_paired[[i]]$`Recipient 3`$WSLH_ID,
                        "Donor_piN_minus_piS"                 = df_AA_HA_list_paired[[i]]$`Donor 0`$piN_minus_piS,
                        "Recipient_piN_minus_piS"             = df_AA_HA_list_paired[[i]]$`Recipient 3`$piN_minus_piS,
                        "Donor_perSite_piN_minus_piS"         = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN_minus_piS,
                        "Recipient_perSite_piN_minus_piS"     = df_AA_HA_list_paired[[i]]$`Recipient 3`$perSite_piN_minus_piS,
                        "Donor_piN"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piN,
                        "Recipient_piN"                       = df_AA_HA_list_paired[[i]]$`Recipient 3`$piN,
                        "Donor_perSite_piN"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piN,
                        "Recipient_perSite_piN"               = df_AA_HA_list_paired[[i]]$`Recipient 3`$perSite_piN,
                        "Donor_piS"                           = df_AA_HA_list_paired[[i]]$`Donor 0`$piS,
                        "Recipient_piS"                       = df_AA_HA_list_paired[[i]]$`Recipient 3`$piS,
                        "Donor_perSite_piS"                   = df_AA_HA_list_paired[[i]]$`Donor 0`$perSite_piS,
                        "Recipient_perSite_piS"               = df_AA_HA_list_paired[[i]]$`Recipient 3`$perSite_piS)
    temp <- rbind(temp1, temp2)
    temp <- rbind(temp, temp3)
  }
  df_AA_HA_paired <- rbind(df_AA_HA_paired, temp)
  df_AA_HA_paired <- na.omit(df_AA_HA_paired)
}

## remove the dummy data frame
df_AA_HA_paired <- df_AA_HA_paired[df_AA_HA_paired$Num!="Remove",]

## change to numeric
df_AA_HA_paired$Recipient_piN                    <- as.numeric(df_AA_HA_paired$Recipient_piN)
df_AA_HA_paired$Recipient_perSite_piN            <- as.numeric(df_AA_HA_paired$Recipient_perSite_piN)
df_AA_HA_paired$Recipient_piS                    <- as.numeric(df_AA_HA_paired$Recipient_piS)
df_AA_HA_paired$Recipient_perSite_piS            <- as.numeric(df_AA_HA_paired$Recipient_perSite_piS)
df_AA_HA_paired$Recipient_piN_minus_piS          <- as.numeric(df_AA_HA_paired$Recipient_piN_minus_piS)
df_AA_HA_paired$Recipient_perSite_piN_minus_piS  <- as.numeric(df_AA_HA_paired$Recipient_perSite_piN_minus_piS)

df_AA_HA_paired$Donor_piN                    <- as.numeric(df_AA_HA_paired$Donor_piN)
df_AA_HA_paired$Donor_perSite_piN            <- as.numeric(df_AA_HA_paired$Donor_perSite_piN)
df_AA_HA_paired$Donor_piS                    <- as.numeric(df_AA_HA_paired$Donor_piS)
df_AA_HA_paired$Donor_perSite_piS            <- as.numeric(df_AA_HA_paired$Donor_perSite_piS)
df_AA_HA_paired$Donor_piN_minus_piS          <- as.numeric(df_AA_HA_paired$Donor_piN_minus_piS)
df_AA_HA_paired$Donor_perSite_piN_minus_piS  <- as.numeric(df_AA_HA_paired$Donor_perSite_piN_minus_piS)

## calculate change between donor and recipient
df_AA_HA_paired <- transform(df_AA_HA_paired, 
                             delta.piN                     = df_AA_HA_paired$Recipient_piN                   - df_AA_HA_paired$Donor_piN,
                             delta.perSite_piN             = df_AA_HA_paired$Recipient_perSite_piN           - df_AA_HA_paired$Donor_perSite_piN,
                             delta.piS                     = df_AA_HA_paired$Recipient_piS                   - df_AA_HA_paired$Donor_piS,
                             delta.perSite_piS             = df_AA_HA_paired$Recipient_perSite_piS           - df_AA_HA_paired$Donor_perSite_piS,
                             delta.piN_minus_piS           = df_AA_HA_paired$Recipient_piN_minus_piS         - df_AA_HA_paired$Donor_piN_minus_piS,
                             delta.perSite_piN_minus_piS   = df_AA_HA_paired$Recipient_perSite_piN_minus_piS - df_AA_HA_paired$Donor_perSite_piN_minus_piS)

## add ID and factor
df_AA_HA_paired$ID <- paste(df_AA_HA_paired$Donor_WSLH_ID, df_AA_HA_paired$Recipient_WSLH_ID, sep="-")
df_AA_HA_paired$ID <- as.factor(df_AA_HA_paired$ID)

## remove duplicate rows
df_AA_HA_paired <- df_AA_HA_paired[!duplicated(df_AA_HA_paired), ]


#### df_shared_iSNVs ####
## merge df_all_iSNVs with md so we can see which samples are donor and recipient
df_shared_iSNVs <- merge(df_all_iSNVs, md,
                         by.x = "WSLH_ID",
                         by.y = "WSLH_ID")

## create iSNV column with REF, POS, ALT, and Gene_Name so we can create iSNV list
df_shared_iSNVs$iSNV <- paste(df_shared_iSNVs$REF, df_shared_iSNVs$POS, df_shared_iSNVs$ALT, sep = "")
df_shared_iSNVs$iSNV <- paste(df_shared_iSNVs$iSNV, df_shared_iSNVs$Gene_Name, sep = "-")

## create identifier and list of df_shared_iSNVs by Identifier
df_shared_iSNVs$Identifier <- paste(df_shared_iSNVs$Num, df_shared_iSNVs$Case, sep="-")
list_case <- split(df_shared_iSNVs, with(df_shared_iSNVs, Identifier, drop=T))

## clean up the list!
n <- length(list_case)
for (i in 1:n) {
  list_case[[i]]$WSLH_ID <- as.factor(as.character(list_case[[i]]$WSLH_ID))
  list_case[[i]]$ID <- as.integer(list_case[[i]]$WSLH_ID)
  list_case[[i]]$ID[list_case[[i]]$Case=="Donor"] <- 0
  list_case[[i]]$ID <- paste(list_case[[i]]$Case, list_case[[i]]$ID, sep=" ")
}

## reduce list to df_shared_iSNVs
df_shared_iSNVs <- Reduce(full_join, list_case)

## create iSNV string of all iSNVs in df_shared_iSNVs
df_shared_iSNVs$iSNV <- as.factor(df_shared_iSNVs$iSNV)
all_iSNVs <- data.frame("iSNV"=levels(df_shared_iSNVs$iSNV), "AF"=0)

## back to list!
list_case <- split(df_shared_iSNVs, with(df_shared_iSNVs, WSLH_ID, drop=T))

## add iSNVs where AF=0 to list_case
n <- length(list_case)
for (i in 1:n) {
  temp <- all_iSNVs
  temp1 <- merge(temp, list_case[[i]][, c("iSNV", "AF")], by="iSNV", all.x=TRUE)
  temp1[is.na(temp1)] <- 0
  temp1$AF <- temp1$AF.y
  temp1$AF.x <- NULL
  temp1$AF.y <- NULL
  temp1$WSLH_ID <- list_case[[i]]$WSLH_ID[1]
  temp1$Case <- list_case[[i]]$Case[1]
  temp1$Num <- list_case[[i]]$Num[1]
  temp1$Case_Num <- paste(temp1$Case, temp1$Num, sep = "-")
  temp1$ID <- list_case[[i]]$ID[1]
  list_case[[i]] <- temp1
}
temp <- NULL
temp1 <- NULL

## back to df
df_shared_iSNVs <- Reduce(full_join, list_case)

## string of all recipients; string of all donors
string_recips <- as.numeric(gsub("Recipient-", "", levels(as.factor(df_shared_iSNVs$Case_Num[df_shared_iSNVs$Case=="Recipient"]))))
string_donors <- as.numeric(gsub("Donor-", "", levels(as.factor(df_shared_iSNVs$Case_Num[df_shared_iSNVs$Case=="Donor"]))))
string_paired <- intersect(string_recips, string_donors); string_paired

## remove Case 1544, which has the same WSLH_ID as two unique 
# demographic observations, making it look like a donor-recipient pair
string_paired <- string_paired[string_paired %!in% c("1544")]

## only keep recipients if they have a donor and vise-versa
df_shared_iSNVs <- filter(df_shared_iSNVs, Num %in% string_paired)

## back to list
list_shared <- split(df_shared_iSNVs, with(df_shared_iSNVs, Num, drop=T))

## split the list dfs into twice-nested dfs within list
n <- length(list_shared)
for (i in 1:n) {
  list_shared[[i]] <- split(list_shared[[i]], with(list_shared[[i]], ID, drop=T))
}

## only one WSLH_ID per donor or recipient in pairing
n <- length(list_shared)
for (i in 1:n) {
  ## donor dups
  list_shared[[i]]$`Donor 0`$WSLH_ID <- as.factor(as.character(list_shared[[i]]$`Donor 0`$WSLH_ID))
  if (length(levels(list_shared[[i]]$`Donor 0`$WSLH_ID)) >= 2) {
    levels <- levels(list_shared[[i]]$`Donor 0`$WSLH_ID)
    list_shared[[i]]$`Donor 0` <- list_shared[[i]]$`Donor 0`[list_shared[[i]]$`Donor 0`$WSLH_ID==levels[1],]
  }
  ## recip 1 dups
  list_shared[[i]]$`Recipient 1`$WSLH_ID <- as.factor(as.character(list_shared[[i]]$`Recipient 1`$WSLH_ID))
  if (length(levels(list_shared[[i]]$`Recipient 1`$WSLH_ID)) >= 2) {
    levels <- levels(list_shared[[i]]$`Recipient 1`$WSLH_ID)
    list_shared[[i]]$`Recipient 1` <- list_shared[[i]]$`Recipient 1`[list_shared[[i]]$`Recipient 1`$WSLH_ID==levels[1],]
  }
  ## recip 2 dups
  list_shared[[i]]$`Recipient 2`$WSLH_ID <- as.factor(as.character(list_shared[[i]]$`Recipient 2`$WSLH_ID))
  if (length(levels(list_shared[[i]]$`Recipient 2`$WSLH_ID)) >= 2) {
    levels <- levels(list_shared[[i]]$`Recipient 2`$WSLH_ID)
    list_shared[[i]]$`Recipient 2` <- list_shared[[i]]$`Recipient 2`[list_shared[[i]]$`Recipient 2`$WSLH_ID==levels[1],]
  }
  ## recip 3 dups
  list_shared[[i]]$`Recipient 3`$WSLH_ID <- as.factor(as.character(list_shared[[i]]$`Recipient 3`$WSLH_ID))
  if (length(levels(list_shared[[i]]$`Recipient 3`$WSLH_ID)) >= 2) {
    levels <- levels(list_shared[[i]]$`Recipient 3`$WSLH_ID)
    list_shared[[i]]$`Recipient 3` <- list_shared[[i]]$`Recipient 3`[list_shared[[i]]$`Recipient 3`$WSLH_ID==levels[1],]
  }
}

## remove empty dfs from nested list
n <- length(list_shared)
for (i in 1:n) {
  list_shared[[i]] <- Filter(function(x) dim(x)[1] > 0, list_shared[[i]])
}

## create df in format for JT plot
df_list <- list()

for (i in seq_along(list_shared)) {
  donor <- list_shared[[i]]$`Donor 0`
  recipient_names <- names(list_shared[[i]])[grepl("Recipient", names(list_shared[[i]]))]
  
  for (r in recipient_names) {
    recipient <- list_shared[[i]][[r]]
    
    shared <- inner_join(donor, recipient, by = "iSNV", suffix = c("_Donor", "_Recipient"))
    
    if (nrow(shared) > 0) {
      df <- data.frame(
        Num = shared$Num_Donor,
        iSNV = shared$iSNV,
        Donor_WSLH_ID = shared$WSLH_ID_Donor,
        Recipient_WSLH_ID = shared$WSLH_ID_Recipient,
        Donor_AF = shared$AF_Donor,
        Recipient_AF = shared$AF_Recipient
      )
      df_list[[length(df_list) + 1]] <- df
    }
  }
}

df_JTplot <- bind_rows(df_list)
df_JTplot <- na.omit(df_JTplot)
df_JTplot$ID <- paste(df_JTplot$Donor_WSLH_ID, df_JTplot$Recipient_WSLH_ID, sep="-")
df_JTplot$ID <- as.factor(df_JTplot$ID)


plot_jtplot <- ggplot() + 
  geom_point(data = df_JTplot, color = "Black", alpha = 0.5, size = 1.5, stroke = NA,
             aes(x = Donor_AF, y = Recipient_AF)) + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) + 
  scale_x_continuous(breaks = seq(0, 1, .10), labels = seq(0, 1, .10), expand = c(0.01, 0.01)) + 
  scale_y_continuous(breaks = seq(0, 1, .10), labels = seq(0, 1, .10), expand = c(0.01, 0.01)) + 
  labs(x="iSNV frequency in Donor", y="iSNV frequency in Recipient", title = title) + 
  theme(legend.title = element_blank(), legend.key = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 6, margin = margin(t = 4)),
        axis.title.y = element_text(size = 6, margin = margin(r = 4)),
        legend.text = element_text(size = 6),
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        panel.border = element_rect(color = "grey", fill = NA, size = .5),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank())

df_trunc <- filter(df_shared_iSNVs, AF < 0.95)
df_trunc <- filter(df_trunc, AF > 0)
temp_df1 <- as.data.frame(table(df_trunc$iSNV, df_trunc$Num))

#### df_prop_shared_random ####
df_prop_shared_random <- data.frame(df_all_iSNVs[4], df_all_iSNVs[2], df_all_iSNVs[5], df_all_iSNVs[21], df_all_iSNVs[14], df_all_iSNVs[27])

## create iSNV column with REF, POS, ALT, and Gene_Name so we can create iSNV list
df_prop_shared_random$iSNV <- paste(df_prop_shared_random$REF, df_prop_shared_random$POS, df_prop_shared_random$ALT, sep = "")
df_prop_shared_random$iSNV <- paste(df_prop_shared_random$iSNV, df_prop_shared_random$Gene_Name, sep = "-")

## create iSNV string of all iSNVs in df_prop_shared_random
df_prop_shared_random$iSNV <- as.factor(df_prop_shared_random$iSNV)
all_iSNVs_random <- data.frame("iSNV"=levels(df_prop_shared_random$iSNV), "AF"=0)

## to list
list_prop_shared_random <- split(df_prop_shared_random, with(df_prop_shared_random, WSLH_ID, drop=T))

## remove blank data frames (there should be none if QC was performed!)
list_prop_shared_random <- Filter(function(x) dim(x)[1] > 0, list_prop_shared_random)

## add iSNVs where AF=0 to list_prop_shared_random
n <- length(list_prop_shared_random)
for (i in 1:n) {
  temp <- all_iSNVs_random
  temp1 <- merge(temp, list_prop_shared_random[[i]][, c("iSNV", "AF")], by="iSNV", all.x=TRUE)
  temp1[is.na(temp1)] <- 0
  temp1$AF <- temp1$AF.y
  temp1$AF.x <- NULL
  temp1$AF.y <- NULL
  temp1$WSLH_ID <- list_prop_shared_random[[i]]$WSLH_ID[1]
  temp1$ID <- list_prop_shared_random[[i]]$ID[1]
  list_prop_shared_random[[i]] <- temp1
}
temp <- NULL
temp1 <- NULL

## set up environment for random pairings
n <- length(list_prop_shared_random)
df_random <- data.frame("iSNV"="", "Donor_AF"=0, "Recipient_AF"=0, "ID"="")
df_random <- df_random[0,]

df_compiled <- data.frame("ID"="", "Prop_Var_Shared"=0)
df_compiled <- df_compiled[0,]

list_random <- vector("list", n)
list_random_shared <- vector("list", n)

## change Donor_AF column to just AF
for (i in 1:n) {names(list_prop_shared_random[[i]])[names(list_prop_shared_random[[i]]) == "Donor_AF"] <- "AF"}

## random pairings (pairwise, not self-pair, only one-way pairs)
x = 1
start = Sys.time()
for (i in 2:n) {
  for (j in 1:x) {
    df_temp <- data.frame(
      "ID"               = paste(names(list_prop_shared_random[i]), names(list_prop_shared_random[j]), sep="-"),
      "iSNV"             = list_prop_shared_random[[i]]$iSNV,
      "Donor_AF"         = list_prop_shared_random[[i]]$AF,
      "Recipient_AF"     = list_prop_shared_random[[j]]$AF,
      "Donor_Binary"     = ifelse(list_prop_shared_random[[i]]$AF > 0, list_prop_shared_random[[i]]$iSNV, NA),
      "Recipient_Binary" = ifelse(list_prop_shared_random[[j]]$AF > 0, list_prop_shared_random[[j]]$iSNV, NA)
    )
    df_temp2 <- data.frame(
      "Shared_Binary"    = ifelse(df_temp$Donor_Binary == df_temp$Recipient_Binary, df_temp$iSNV, 0),
      "Donor_Count"      = length(which(df_temp$Donor_AF != 0)),
      "Recipient_Count"  = length(which(df_temp$Recipient_AF != 0))
    )
    df_temp3 <- data.frame("Binary_Count"     = length(which(df_temp2$Shared_Binary != 0)))
    df_temp4 <- data.frame("Total_iSNVs"      = sum(df_temp2$Donor_Count[1] + df_temp2$Recipient_Count[1]))
    df_temp5 <- data.frame("Prop_Var_Shared"  = sum(sum(2 * df_temp3$Binary_Count[1]) / df_temp4$Total_iSNVs[1]))
    df_temp <- cbind(df_temp, df_temp2, df_temp3, df_temp4, df_temp5)
    list_random[[i-1]][[j]] <- rbind(df_random, df_temp)
    df_temp6 <- data.frame(
      "ID" = df_temp$ID[1], 
      "Prop_Var_Shared" = df_temp$Prop_Var_Shared[1]
    )
    df_compiled <- rbind(df_compiled, df_temp6)
  }
  out <- paste(i, "/", n, ": pairing pairs..."); print(out)
  x = x + 1
}
finish <- Sys.time()
out1 <- paste("Start: ", start)
out2 <- paste("Finish:", finish)
print(out1); print(out2)
df_compiled$Pairing <- "Random"

## create bin range
bin_range <- data.frame("Lower" = seq(0,0.95,0.05), "Upper" = seq(0.05,1,0.05))
bin_range$Bins <- paste(bin_range[,1], bin_range[,2], sep="-")

## add Prop_Var_Shared bins
# this isn't the fastest method, but it's adjustable (bin_range above)
# these nested loops parse every iSNV and bin it based on its AF
n <- length(df_compiled[,1])
m <- length(bin_range[,1])
df_compiled$Bins <- ""
start = Sys.time()
for (i in 1:n) {
  for (j in 1:m) {
    df_compiled[i,]$Bins[df_compiled[i,]$Prop_Var_Shared >= bin_range[j,1] && 
                           df_compiled[i,]$Prop_Var_Shared < bin_range[j,2]] <- bin_range[j,3]
    df_compiled[i,]$Bins[df_compiled[i,]$Prop_Var_Shared == 0.3] <- "0.3-0.35"
    df_compiled[i,]$Bins[df_compiled[i,]$Prop_Var_Shared == 1.0] <- "0.95-1"
    # not sure why 0.3 isn't mapping
  }
  out <- paste(i, "/", n, ": binning iSNVs..."); print(out)
}
finish <- Sys.time()
out1 <- paste("Start: ", start)
out2 <- paste("Finish:", finish)
print(out1); print(out2)

df_compiled$Pairing <- "Random"

## split df into bin and bin frequency
df_compiled <- as.data.frame(table(df_compiled$Bins, df_compiled$Pairing))
list_compiled <- split(df_compiled, with(df_compiled, Var2, drop=T))
n <- length(list_compiled)
m <- length(list_compiled[[1]][,1])
for (i in 1:n) {
  for (j in 1:m) {
    list_compiled[[i]]$Prop[j] <- sum(list_compiled[[i]][j,3]/sum(list_compiled[[i]][,3]))
  }
}

## back to df
df_compiled <- Reduce(full_join, list_compiled)

#### df_prop_shared ####
## split df to list by ID (each pair)
df_prop_shared <- data.frame(df_JTplot[2],df_JTplot[5],df_JTplot[6],df_JTplot[7])
df_prop_shared$Pairing <- "Household"
list_prop_shared <- split(df_prop_shared, with(df_prop_shared, ID, drop=T))

## counts
n <- length(list_prop_shared)
for (i in 1:n) {
  list_prop_shared[[i]] <- list_prop_shared[[i]] %>%
    mutate(Donor_Binary = ifelse(Donor_AF > 0, iSNV, 0)) %>% 
    mutate(Recip_Binary = ifelse(Recipient_AF > 0, iSNV, NA)) %>% 
    mutate(Shared_Binary = ifelse(Donor_Binary == Recip_Binary, iSNV, 0))
  list_prop_shared[[i]]$Donor_Count <- length(which(list_prop_shared[[i]]$Donor_AF != 0))
  list_prop_shared[[i]]$Recip_Count <- length(which(list_prop_shared[[i]]$Recipient_AF != 0))
  list_prop_shared[[i]]$Binary_Count <- length(which(list_prop_shared[[i]]$Shared_Binary != 0))
  list_prop_shared[[i]]$Total_iSNVs <- sum(list_prop_shared[[i]]$Donor_Count[1] + list_prop_shared[[i]]$Recip_Count[1])
  list_prop_shared[[i]]$Prop_Var_Shared <- sum(sum(2 * list_prop_shared[[i]]$Binary_Count[1]) / list_prop_shared[[i]]$Total_iSNVs[1])
}

df_prop_shared <-  list_prop_shared[[1]][0,]
for (i in 1:n) {
  df_prop_shared <- rbind(df_prop_shared, list_prop_shared[[i]][1,])
}
df_prop_shared$iSNV <- NULL
df_prop_shared$Donor_AF <- NULL
df_prop_shared$Recipient_AF <- NULL
df_prop_shared$Donor_Binary <- NULL
df_prop_shared$Recip_Binary <- NULL
df_prop_shared$Shared_Binary <- NULL

df_household <- df_prop_shared

## add Prop_Var_Shared bins
n <- length(df_prop_shared[,1])
m <- length(bin_range[,1])
df_prop_shared$Bins <- ""
for (i in 1:n) {
  for (j in 1:m)
    if(df_prop_shared[i,]$Prop_Var_Shared > bin_range[j,1] &
       df_prop_shared[i,]$Prop_Var_Shared <= bin_range[j,2])
      df_prop_shared[i,]$Bins <- bin_range[j,3]
}
df_prop_shared$Bins[df_prop_shared$Prop_Var_Shared==0] <- "0-0.05"

## split df by 
df_prop_shared <- as.data.frame(table(df_prop_shared$Bins, df_prop_shared$Pairing))
list_prop_shared <- split(df_prop_shared, with(df_prop_shared, Var2, drop=T))
n <- length(list_prop_shared)
m <- length(list_prop_shared[[1]][,1])
for (i in 1:n) {
  for (j in 1:m) {
    list_prop_shared[[i]]$Prop[j] <- sum(list_prop_shared[[i]][j,3]/sum(list_prop_shared[[i]][,3]))
  }
}

## back to df
df_prop_shared_household <- Reduce(full_join, list_prop_shared)

df_prop_shared <- rbind(df_compiled, df_prop_shared_household)

Fig5alabels <- c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5",
                 "0.6", "0.7", "0.8", "0.9", "1.0")

df_prop_shared$Var1 <- factor(df_prop_shared$Var1, levels = levels(as.factor(bin_range$Bins)))

#### sg_paired ####
## merge with md
sg_paired <- merge(sg, md, by.x="file", by.y="WSLH_ID")

# remove duplicates if present
sg_paired <- sg_paired[!duplicated(sg_paired), ]

## factor
sg_paired$Case <- as.factor(sg_paired$Case)
sg_paired$Case_Num <- as.factor(sg_paired$Case_Num)

## find samples with pairs by pulling by group name and excluding other group, reciprocally
# recips
string_sg_recips <- as.character(gsub("Recipient ", "", levels(as.factor(sg_paired$Case_Num[sg_paired$Case=="Recipient"]))))
string_sg_recips <- as.numeric(gsub("Donor .*", "0", string_sg_recips))
string_sg_recips <- string_sg_recips[string_sg_recips>0]
# donors
string_sg_donors <- as.character(gsub("Donor ", "", levels(as.factor(sg_paired$Case_Num[sg_paired$Case=="Donor"]))))
string_sg_donors <- as.numeric(gsub("Recipient .*", "0", string_sg_donors))
string_sg_donors <- string_sg_donors[string_sg_donors>0]
# intersection
string_sg_pairs <- intersect(string_sg_recips, string_sg_donors); string_sg_pairs

## remove Case 1544, which has the same WSLH_ID as two unique 
# demographic observations, making it look like a donor-recipient pair
string_sg_pairs <- string_sg_pairs[string_sg_pairs %!in% c("1544")]

## only sg_paired data that are pairs
sg_paired <- filter(sg_paired, Num %in% string_sg_pairs)

## factor
sg_paired$file <- as.factor(as.character(sg_paired$file))
sg_paired$Case <- as.factor(as.character(sg_paired$Case))
sg_paired$Case_Num <- as.factor(as.character(sg_paired$Case_Num))

## create identifier and list of sg_paired by Identifier
sg_paired$Identifier <- paste(sg_paired$Num, sg_paired$Case, sep="-")
sg_list_paired <- split(sg_paired, with(sg_paired, Identifier, drop=T))

## clean up the list!
n <- length(sg_list_paired)
for (i in 1:n) {
  sg_list_paired[[i]]$WSLH_ID <- as.factor(as.character(sg_list_paired[[i]]$file))
  sg_list_paired[[i]]$ID <- as.integer(sg_list_paired[[i]]$WSLH_ID)
  sg_list_paired[[i]]$ID[sg_list_paired[[i]]$Case=="Donor"] <- 0
  sg_list_paired[[i]]$ID <- paste(sg_list_paired[[i]]$Case, sg_list_paired[[i]]$ID, sep=" ")
}

## back to df
sg_paired <- Reduce(full_join, sg_list_paired)

## Gene column
sg_paired$Gene <- sg_paired$product

## all_genes
all_genes <- data.frame("Gene"=levels(as.factor(as.character(sg_paired$product))), 
                        "piN_minus_piS"=0, "pi"=0)

## back to list
sg_list_paired <- split(sg_paired, with(sg_paired, WSLH_ID, drop=T))

## add Gene where piN_minus_piS=0 to sg_list_paired
n <- length(sg_list_paired)
for (i in 1:n) {
  temp <- all_genes
  temp1 <- merge(temp, sg_list_paired[[i]][, c("Gene", "piN_minus_piS", "pi", "piN", "piS")], by="Gene", all.x=TRUE)
  temp1[is.na(temp1)] <- 0
  temp1$piN_minus_piS <- temp1$piN_minus_piS.y
  temp1$piN_minus_piS.x <- NULL
  temp1$piN_minus_piS.y <- NULL
  temp1$pi <- temp1$pi.y
  temp1$pi.x <- NULL
  temp1$pi.y <- NULL
  temp1$WSLH_ID <- sg_list_paired[[i]]$WSLH_ID[1]
  temp1$Case <- sg_list_paired[[i]]$Case[1]
  temp1$Num <- sg_list_paired[[i]]$Num[1]
  temp1$Case_Num <- paste(temp1$Case, temp1$Num, sep = "-")
  temp1$ID <- sg_list_paired[[i]]$ID[1]
  sg_list_paired[[i]] <- temp1
}
temp <- NULL
temp1 <- NULL

## back to df
sg_paired <- Reduce(full_join, sg_list_paired)

## back to list
sg_list_paired <- split(sg_paired, with(sg_paired, Num, drop=T))

## split the list dfs into twice-nested dfs within list
n <- length(sg_list_paired)
for (i in 1:n) {
  sg_list_paired[[i]] <- split(sg_list_paired[[i]], with(sg_list_paired[[i]], ID, drop=T))
}

## only one WSLH_ID per donor or recipient in pairing
n <- length(sg_list_paired)
for (i in 1:n) {
  ## donor dups
  sg_list_paired[[i]]$`Donor 0`$WSLH_ID <- as.factor(as.character(sg_list_paired[[i]]$`Donor 0`$WSLH_ID))
  if (length(levels(sg_list_paired[[i]]$`Donor 0`$WSLH_ID)) >= 2) {
    levels <- levels(sg_list_paired[[i]]$`Donor 0`$WSLH_ID)
    sg_list_paired[[i]]$`Donor 0` <- sg_list_paired[[i]]$`Donor 0`[sg_list_paired[[i]]$`Donor 0`$WSLH_ID==levels[1],]
  }
  ## recip 1 dups
  sg_list_paired[[i]]$`Recipient 1`$WSLH_ID <- as.factor(as.character(sg_list_paired[[i]]$`Recipient 1`$WSLH_ID))
  if (length(levels(sg_list_paired[[i]]$`Recipient 1`$WSLH_ID)) >= 2) {
    levels <- levels(sg_list_paired[[i]]$`Recipient 1`$WSLH_ID)
    sg_list_paired[[i]]$`Recipient 1` <- sg_list_paired[[i]]$`Recipient 1`[sg_list_paired[[i]]$`Recipient 1`$WSLH_ID==levels[1],]
  }
  ## recip 2 dups
  sg_list_paired[[i]]$`Recipient 2`$WSLH_ID <- as.factor(as.character(sg_list_paired[[i]]$`Recipient 2`$WSLH_ID))
  if (length(levels(sg_list_paired[[i]]$`Recipient 2`$WSLH_ID)) >= 2) {
    levels <- levels(sg_list_paired[[i]]$`Recipient 2`$WSLH_ID)
    sg_list_paired[[i]]$`Recipient 2` <- sg_list_paired[[i]]$`Recipient 2`[sg_list_paired[[i]]$`Recipient 2`$WSLH_ID==levels[1],]
  }
  ## recip 3 dups
  sg_list_paired[[i]]$`Recipient 3`$WSLH_ID <- as.factor(as.character(sg_list_paired[[i]]$`Recipient 3`$WSLH_ID))
  if (length(levels(sg_list_paired[[i]]$`Recipient 3`$WSLH_ID)) >= 2) {
    levels <- levels(sg_list_paired[[i]]$`Recipient 3`$WSLH_ID)
    sg_list_paired[[i]]$`Recipient 3` <- sg_list_paired[[i]]$`Recipient 3`[sg_list_paired[[i]]$`Recipient 3`$WSLH_ID==levels[1],]
  }
}

## remove empty dfs and single donors/recips from nested list
# empty
n <- length(sg_list_paired)
for (i in 1:n) {
  sg_list_paired[[i]] <- Filter(function(x) dim(x)[1] > 0, sg_list_paired[[i]])
}
# donor/recips-only by identifying which list of lists is problematic
to_remove <- c(NA)
n <- length(sg_list_paired)
for (i in 1:n) {
  if (length(sg_list_paired[[i]])==1) {
    to_remove <- c(i, to_remove)
  }
}
# removes identified donor/recip-only cases
if (length(to_remove[!is.na(to_remove)])>0) {
  sg_list_paired <- sg_list_paired[-to_remove[!is.na(to_remove)]]
}

## create df in pairing format
# Initialize a list to collect rows
paired_list <- list()

# Loop over households
for (i in seq_along(sg_list_paired)) {
  donor <- sg_list_paired[[i]]$`Donor 0`
  recipient_names <- names(sg_list_paired[[i]])[grepl("Recipient", names(sg_list_paired[[i]]))]
  
  # Loop over all recipients in the current household
  for (r in recipient_names) {
    recipient <- sg_list_paired[[i]][[r]]
    
    # Create a paired dataframe for this donor-recipient pair
    paired_df <- data.frame(
      Num                   = donor$Num,
      Gene                  = donor$Gene,
      Donor_WSLH_ID         = donor$WSLH_ID,
      Recipient_WSLH_ID     = recipient$WSLH_ID,
      Donor_piN_minus_piS   = donor$piN_minus_piS,
      Recipient_piN_minus_piS = recipient$piN_minus_piS,
      Donor_piN             = donor$piN,
      Recipient_piN         = recipient$piN,
      Donor_piS             = donor$piS,
      Recipient_piS         = recipient$piS,
      Donor_pi              = donor$pi,
      Recipient_pi          = recipient$pi
    )
    
    paired_list[[length(paired_list) + 1]] <- paired_df
  }
}

# Combine all rows and drop any with NAs
sg_paired <- do.call(rbind, paired_list)
sg_paired <- na.omit(sg_paired)

## calculate change between donor and recipient
sg_paired <- transform(sg_paired, 
                       delta.piN_minus_piS     = sg_paired$Recipient_piN_minus_piS - sg_paired$Donor_piN_minus_piS,
                       delta.piN               = sg_paired$Recipient_piN           - sg_paired$Donor_piN,
                       delta.piS               = sg_paired$Recipient_piS           - sg_paired$Donor_piS,
                       delta.pi                = sg_paired$Recipient_pi            - sg_paired$Donor_pi)

## add ID and factor
sg_paired$ID <- paste(sg_paired$Donor_WSLH_ID, sg_paired$Recipient_WSLH_ID, sep="-")
sg_paired$ID <- as.factor(sg_paired$ID)

## remove duplicate rows
sg_paired <- sg_paired[!duplicated(sg_paired), ]

## ds_sg_paired
ds_sg_paired_piN_minus_piS     <- as.data.frame(ds(sg_paired, varname="delta.piN_minus_piS", groupnames=c("Gene")))
ds_sg_paired_pi         <- as.data.frame(ds(sg_paired, varname="delta.pi", groupnames=c("Gene")))
ds_sg_paired_piN         <- as.data.frame(ds(sg_paired, varname="delta.piN", groupnames=c("Gene")))
ds_sg_paired_piS         <- as.data.frame(ds(sg_paired, varname="delta.piS", groupnames=c("Gene")))
colnames(ds_sg_paired_piN_minus_piS)     <- c("Gene", "delta.piN_minus_piS", "sd.piN_minus_piS", "se.piN_minus_piS")
colnames(ds_sg_paired_pi)         <- c("Gene", "delta.pi", "sd.pi", "se.pi")
colnames(ds_sg_paired_piN)         <- c("Gene", "delta.piN", "sd.piN", "se.piN")
colnames(ds_sg_paired_piS)         <- c("Gene", "delta.piS", "sd.piS", "se.piS")
ds_sg_paired <- merge(ds_sg_paired_piN_minus_piS, ds_sg_paired_pi, by="Gene")
ds_sg_paired <- merge(ds_sg_paired, ds_sg_paired_piN, by="Gene")
ds_sg_paired <- merge(ds_sg_paired, ds_sg_paired_piS, by="Gene")

### logfile ####
## header
str_header_dir    <- paste("VCF location:", dir, sep = "\t")
str_header_dir_sg <- paste("SNPGenie location:", dir_sg, sep = "\t")
str_header_dir_s  <- paste("Save location:", dir_s, sep = "\t")
str_header_dir_md <- paste("Metadata location:", dir_md, sep = "\t")
str_header_minAF  <- paste("Minimum Allele Freqncy:", minAF, sep = "\t")
## number of samples passing qc
str_n_samples <- paste("Number of samples:",length(list), sep = "\t")
## number of cases with pairs
str_n_cases_with_pairs <- paste("Number of cases w/pairs:", length(string_paired), sep = "\t")
## number of pairing events
str_n_pairing_events <- paste("Number of pairing events:", length(df_household$ID), sep = "\t")
## cases with pairs
str_cases_with_pairs <- paste("Cases w/pairs:",paste(string_paired, collapse = ", "), sep = "\t")


## logfile
logfile <- c(str_header_dir,str_header_dir_sg,str_header_dir_s ,str_header_dir_md, 
             str_n_samples, str_n_cases_with_pairs, str_n_pairing_events, str_cases_with_pairs); logfile

#### SAVE ####
setwd(dir_s); getwd(); dir()

## iSNV plot
tiff(filename = "0-plot_iSNVs.tiff", 
     width = 5, height = 5, units = "in", res = 750)
plot(table(df_all_iSNVs$WSLH_ID), 
     xlab = "Sample ID", ylab = "Number of iSNVs", main = title)
dev.off()

ggsave("0-plot_JTplot.tiff", plot_jtplot,
       width = 3, height = 3, 
       units = "in", device='tiff', dpi=750)
write.csv(df_all_iSNVs,                        "df_all_iSNVs.csv")
write.csv(df_sub50AF,                          "df_sub50AF.csv")
write.csv(df_sub50AF_md_n_SNVs_by_subtype,     "df_sub50AF_md_n_SNVs_by_subtype.csv")
write.csv(df_sub50AF_md_n_SNVs_by_vax,         "df_sub50AF_md_n_SNVs_by_vax.csv")
write.csv(df_mut_bins_prop_noNeut,             "df_mut_bins_prop_noNeut.csv")
write.csv(df_compiled,                         "df_compiled.csv")
write.csv(df_divergence,                       "df_divergence.csv")
write.csv(ds_divergence,                       "ds_divergence.csv")
write.csv(df_household,                        "df_household.csv")
write.csv(df_prop_shared,                      "df_prop_shared.csv")
write.csv(df_JTplot,                           "df_JTplot.csv")
write.csv(sg,                                  "sg.csv")
write.csv(ds_sg_paired,                        "ds_sg_paired.csv")
write.csv(sg_paired,                           "sg_paired.csv")
write.csv(df_mut_bins_prop,                    "df_mut_bins_prop.csv")
write_tsv(as.data.frame(logfile),              "logfile.tsv")
write.csv(ds_sg_sw,                            "ds_sg_sw.csv")
write.csv(sg_sw,                               "sg_sw.csv")
write.csv(df_sg_AA_HA,                         "df_sg_AA_HA.csv")
write.csv(df_AA_HA,                            "df_AA_HA.csv")
write.csv(df_AA_HA_paired,                     "df_AA_HA_paired.csv")
write.csv(df_sub50AF_md_n_SNVs_by_bin_and_mut, "df_sub50AF_md_n_SNVs_by_bin_and_mut.csv")
write.csv(df_mut_bins_prop_all,                "df_mut_bins_prop_all.csv")


#### ####