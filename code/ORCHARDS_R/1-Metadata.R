## Clear Global Environment
rm(list = ls())

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
if(!require(ggpubr)){
  install.packages("ggpubr",dependencies = T)
  library(ggpubr)
}
if(!require(cowplot)){
  install.packages("cowplot",dependencies = T)
  library(cowplot)
}
if(!require(grid)){
  install.packages("grid",dependencies = T)
  library(grid)
}
if(!require(gridExtra)){
  install.packages("gridExtra",dependencies = T)
  library(gridExtra)
}
if(!require(ggpubr)){
  install.packages("ggpubr",dependencies = T)
  library(ggpubr)
}

#### Functions ####
## User-defined functions
'%!in%' <- function(x,y)!('%in%'(x,y))
# useful for removing column values of y from df x
# e.g., x <- filter(x,y %!in% c("unwantedValue1", "unwantedValue2"))

#### Plot standards ####
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
# axis_formatting + legend_formatting + background_formatting

## plot_grid
Size_adjust = 12
LR_adjust = -0.5 # less = right
UD_adjust = 1.1 # less = down 

#### Import ####
setwd("/Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/summary"); dir()
dir_s <- paste("/Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs")

## read csv
sampleList <- as.data.frame(read_csv("orchards_sample_list.csv"))
metadata <- as.data.frame(read_csv("0-ORCHARDS_20172019_data_tall.csv"))
samples_QC_pass <- as.data.frame(read_tsv("samples_QC_pass.tsv"))

#### sampleList - Clean ####
## Rename columns
sampleList <- dplyr::rename(sampleList,Study_ID=`Study ID`, 
                            DoC=`Date of Collection`,
                            Type_Subtype=`Influenza PCR Result`, 
                            WSLH_ID=`WSLH ID`, 
                            Ct=`FLU PCR Ct`)

## only keep samples that passed QC
sampleList <- sampleList[sampleList$WSLH_ID %in% samples_QC_pass$WSLH_ID,]

## TIME
# Remove time from DoC
sampleList <- separate(sampleList, "DoC", c("DoC", "H"), sep = " ")
print("NAs are okay someitmes")
sampleList$H <- NULL

## convert DoC to time format
sampleList$DoC_mdy <- mdy(sampleList$DoC)
sampleList$DoC <- sampleList$DoC_mdy
sampleList$DoC_mdy <- NULL
sampleList$DoC <- as.Date(sampleList$DoC)

## Split duplicate rows
# split rows with two Ct values and Types
sampleList <- separate_rows(sampleList,"Type_Subtype",sep=",\\s+",convert=FALSE)
# Split Ct by A and B
sampleList$Ct_backup <- sampleList$Ct
sampleList <- sampleList %>%
  separate("Ct",c("Ct_A","Ct_B"),sep=", ")
sampleList$Ct <- sampleList$Ct_backup
sampleList$Ct_backup <- NULL
# Move A values to A and B values to B
sampleList$Ct <- ifelse(sampleList$Type_Subtype %in% c("Influenza A", "Influenza A H1N1", "Influenza A H3N2"),sampleList$Ct_A,sampleList$Ct)
sampleList$Ct <- ifelse(sampleList$WSLH_ID %in% c("REDACTED", "REDACTED", "REDACTED","REDACTED") & 
                          sampleList$Type_Subtype %in% c("Influenza B", "Influenza B (Victoria)", "Influenza B (Yamagata)"),sampleList$Ct_B,sampleList$Ct)
sampleList$Ct_A <- NULL
sampleList$Ct_B <- NULL
# Remove some leftover characters
sampleList$Ct <- gsub(" \\(A)","",sampleList$Ct)
sampleList$Ct <- gsub(" \\(B)","",sampleList$Ct)
sampleList$Ct <- as.numeric(sampleList$Ct)

## NOMECLATURE
# Factor
sampleList$Type_Subtype <- as.factor(sampleList$Type_Subtype)
# Full name
sampleList$Type_Subtype <- gsub("Influenza B (Y)", "Influenza B (Yamagata)",sampleList$Type_Subtype,fixed=TRUE)
sampleList$Type_Subtype <- gsub("Influenza A H1N1", "Influenza A (H1N1)",sampleList$Type_Subtype,fixed=TRUE)
sampleList$Type_Subtype <- gsub("Influenza A H3N2", "Influenza A (H3N2)",sampleList$Type_Subtype,fixed=TRUE)
# Type
sampleList$Type <- sampleList$Type_Subtype
sampleList$Type <- gsub("Influenza A \\(H3N2)","Influenza A",sampleList$Type)
sampleList$Type <- gsub("Influenza A \\(H1N1)","Influenza A",sampleList$Type)
sampleList$Type <- gsub("Influenza B \\(Yamagata)","Influenza B",sampleList$Type)
sampleList$Type <- gsub("Influenza B \\(Victoria)","Influenza B",sampleList$Type)
# Subtype
sampleList$Subtype <- sampleList$Type_Subtype
sampleList$Subtype <- gsub("Influenza A \\(H3N2)","H3N2",sampleList$Subtype)
sampleList$Subtype <- gsub("Influenza A \\(H1N1)","H1N1",sampleList$Subtype)
sampleList$Subtype <- gsub("Influenza B \\(Yamagata)","Yamagata",sampleList$Subtype)
sampleList$Subtype <- gsub("Influenza B \\(Victoria)","Victoria",sampleList$Subtype)
# Reorder df
sampleList <- sampleList[,c(1,4,2,5,6,7,3)]

## DUPLICATES
# see samples with same WSLH_ID and Ct
dups_removed <- sampleList[duplicated(sampleList[c("WSLH_ID","Ct")]),]
# remove one of two duplicates with same WSLH_ID and Ct 
sampleList <- sampleList[!duplicated(sampleList[c("WSLH_ID","Ct")]),]
# identify Study_ID dups—samples run twice
run_twice <- data.frame("Study_ID" = sampleList$Study_ID[duplicated(sampleList["Study_ID"])],
                        "Longitudinal" = "Yes")
# Mark longitudinal samples
sampleList$Longitudinal <- "No"
sampleList$Longitudinal[sampleList$Study_ID %in% run_twice$Study_ID] <- "Yes"

# create co-infected df
coinf <- unique(sampleList[duplicated(sampleList$WSLH_ID), "WSLH_ID"]); coinf <- coinf$WSLH_ID
coinfList <- filter(sampleList,WSLH_ID %in% coinf)
# filter out co-infected from sampleList
sampleList <- filter(sampleList,WSLH_ID %!in% coinf)

### df for sorting
df_sorting <- sampleList
# convert DoC to year and month
df_sorting$DoC <- ymd(df_sorting$DoC)
df_sorting$DoC_Y <- year(df_sorting$DoC)
df_sorting$DoC_M <- month(df_sorting$DoC)
df_sorting$DoC_D <- day(df_sorting$DoC)
df_sorting$DoC_MY <- paste(df_sorting$DoC_M, "-",df_sorting$DoC_Y)
df_sorting$DoC_MY <- gsub(" ", "", df_sorting$DoC_MY)
# Season column
df_sorting$Season <- ""
df_sorting$Season[df_sorting$DoC_Y==2016 & df_sorting$DoC_M>6] <- "16-17"
df_sorting$Season[df_sorting$DoC_Y==2017 & df_sorting$DoC_M<6] <- "16-17"
df_sorting$Season[df_sorting$DoC_Y==2017 & df_sorting$DoC_M>6] <- "17-18"
df_sorting$Season[df_sorting$DoC_Y==2018 & df_sorting$DoC_M<6] <- "17-18"
df_sorting$Season[df_sorting$DoC_Y==2018 & df_sorting$DoC_M>6] <- "18-19"
df_sorting$Season[df_sorting$DoC_Y==2019 & df_sorting$DoC_M<6] <- "18-19"
# clean
df_sorting <- df_sorting[,c(2,6,13)]
df_sorting$ID <- as.factor(paste(df_sorting$Season, df_sorting$Subtype, sep = "_"))

#### sampleList - samples to sort ####
string_1718_H3N2 <- df_sorting$WSLH_ID[df_sorting$ID=="17-18_H3N2"]; paste(length(string_1718_H3N2), "17-18_H3N2")
string_1819_H3N2 <- df_sorting$WSLH_ID[df_sorting$ID=="18-19_H3N2"]; paste(length(string_1819_H3N2), "18-19_H3N2")
string_1718_H1N1 <- df_sorting$WSLH_ID[df_sorting$ID=="17-18_H1N1"]; paste(length(string_1718_H1N1), "17-18_H1N1")
string_1819_H1N1 <- df_sorting$WSLH_ID[df_sorting$ID=="18-19_H1N1"]; paste(length(string_1819_H1N1), "18-19_H1N1")


#### sampleList - Map pairs ####
pairList <- sampleList
# create Recipient number column value
pairList$Pair_num <- pairList$Study_ID
pairList$Pair_num <- gsub(".*-","Recipient ",pairList$Pair_num)
# truncate Pair_num to be shared Study_ID
pairList$Pair <- pairList$Study_ID
pairList$Pair <- gsub("-.*","",pairList$Pair)
# create df of rows sharing Study_ID
pairList$Pair <- as.integer(pairList$Pair)
n_occur <- data.frame(table(pairList$Pair))
n_occur <- filter(n_occur, Freq>1)
Pairs <- as.integer(as.character(n_occur$Var1))
pair_df <- filter(pairList, Pair %in% Pairs)
n_occur <- NULL
Pairs <- NULL
pairList <- NULL
# create Donor number column value
exceptions <- c("Recipient 1","Recipient 2","Recipient 3","Recipient 4","Recipient 5","Recipient 6")
pair_df <- pair_df %>%
  mutate(Case=if_else(!(Pair_num %in% exceptions),"Donor",Pair_num))
pair_df$Pair_num <- NULL
pair_df$Case_num <- pair_df$Case
pair_df$Case_Pair <- pair_df$Case
pair_df$Case_Pair <- gsub(" 1","",pair_df$Case_Pair)
pair_df$Case_Pair <- gsub(" 2","",pair_df$Case_Pair)
pair_df$Case_Pair <- gsub(" 3","",pair_df$Case_Pair)
pair_df$Case_Pair <- gsub(" 4","",pair_df$Case_Pair)
pair_df$Case_Pair <- gsub(" 5","",pair_df$Case_Pair)
pair_df$Case_Pair <- gsub(" 6","",pair_df$Case_Pair)
pair_df$Case <- pair_df$Case_Pair
pair_df$Case_Pair <- paste(pair_df$Case_Pair,pair_df$Pair)
pair_df <- arrange(pair_df,Pair)
# count pairs
n_Donors <- length(which(pair_df$Case=="Donor")); print(paste(n_Donors, "donors"))
n_Recipients <- length(which(pair_df$Case=="Recipient")); print(paste(n_Recipients, "recipients"))
# factor
pair_df$Type <- as.factor(pair_df$Type)
pair_df$Subtype <- as.factor(pair_df$Subtype)
pair_df$Type_Subtype <- as.factor(pair_df$Type_Subtype)
pair_df$Pair <- as.factor(pair_df$Pair)
pair_df$Case <- as.factor(pair_df$Case)
pair_df$Case_num <- as.factor(pair_df$Case_num)
pair_df$Case_Pair <- as.factor(pair_df$Case_Pair)
# metrics for plot
pair_plot_df <- pair_df
pair_plot_df$value <- pair_plot_df$Case_num
pair_plot_df$value <- gsub("Recipient ","",pair_plot_df$value)
pair_plot_df$value <- gsub("Donor",0,pair_plot_df$value)
pair_plot_df$value <- as.integer(pair_plot_df$value)

## pairs
test <- pair_plot_df
test$WSLH_ID_trunc <- ""
test$WSLH_ID_trunc <- as.integer(gsub(".*VR", "", test$WSLH_ID))
test <- test[!duplicated(test$Study_ID),]
test$WSLH_ID_trunc <- NULL
pair_plot_df <- test

# alter values for plot 
pair_plot_df[pair_plot_df$Pair==665 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==720 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==726 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==730 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==738 & pair_plot_df$Case_num=="Recipient 5", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==940 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==941 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==942 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==948 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==966 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==977 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1002 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1002 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1003 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1025 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1031 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1033 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1036 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1048 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1171 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1179 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1186 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1194 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1205 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1208 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1220 & pair_plot_df$Case_num=="Recipient 6", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1228 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1236 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1251 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1271 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1272 & pair_plot_df$Case_num=="Recipient 5", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1536 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1573 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1574 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1585 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1621 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1636 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1641 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1641 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1656 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1673 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1673 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1696 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1696 & pair_plot_df$Case_num=="Recipient 4", "value"] <- 2
pair_plot_df[pair_plot_df$Pair==1705 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1707 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1736 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1750 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1761 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1769 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1772 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1773 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1787 & pair_plot_df$Case_num=="Recipient 5", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1808 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1817 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1836 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1845 & pair_plot_df$Case_num=="Recipient 2", "value"] <- 1
pair_plot_df[pair_plot_df$Pair==1850 & pair_plot_df$Case_num=="Recipient 3", "value"] <- 1

# alter Case_num to reflect revised Recipient numbers
pair_plot_df$Case_num <- paste(pair_plot_df$Case, pair_plot_df$value)
pair_plot_df$Case_num <- gsub(" 0","",pair_plot_df$Case_num)

# remove rows with no Donor or Recipient (double-timepoint Donors/Recipients)
#pair_plot_df <- pair_plot_df[duplicated(pair_plot_df[c("Case_Pair", "Case_num")]),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="924"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="944"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="953"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="967"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="995"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1056"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1066"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1082"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1106"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1132"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1168"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1192"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1203"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1246"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1252"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1254"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1255"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1270"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1274"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1285"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1286"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1544"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1584"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1590"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1607"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1622"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1626"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1645"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1655"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1665"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1685"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1687"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1699"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1716"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1756"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1762"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1818"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1831"),]
#pair_plot_df <- pair_plot_df[!(pair_plot_df$Pair=="1190"),] #B Donor, A Recipient
pair_df <- pair_plot_df

# duplicate Donors for Recipients > 1
pair_plot_df$line <- pair_plot_df$value
pair_plot_df[pair_plot_df$line==0,"line"] <- 1
Donor_643    <- filter(pair_plot_df,Pair=="643" & Case_num=="Donor")
Donor_941_1  <- filter(pair_plot_df,Pair=="941" & Case_num=="Donor")
Donor_941_2  <- filter(pair_plot_df,Pair=="941" & Case_num=="Donor")
Donor_942_1  <- filter(pair_plot_df,Pair=="942" & Case_num=="Donor")
Donor_942_2  <- filter(pair_plot_df,Pair=="942" & Case_num=="Donor")
Donor_948    <- filter(pair_plot_df,Pair=="948" & Case_num=="Donor")
Donor_951    <- filter(pair_plot_df,Pair=="951" & Case_num=="Donor")
Donor_1002_1 <- filter(pair_plot_df,Pair=="1002" & Case_num=="Donor")
Donor_1002_2 <- filter(pair_plot_df,Pair=="1002" & Case_num=="Donor")
Donor_1031   <- filter(pair_plot_df,Pair=="1031" & Case_num=="Donor")
Donor_1033_1 <- filter(pair_plot_df,Pair=="1033" & Case_num=="Donor")
Donor_1033_2 <- filter(pair_plot_df,Pair=="1033" & Case_num=="Donor")
Donor_1048   <- filter(pair_plot_df,Pair=="1048" & Case_num=="Donor")
Donor_1162_1 <- filter(pair_plot_df,Pair=="1162" & Case_num=="Donor")
Donor_1162_2 <- filter(pair_plot_df,Pair=="1162" & Case_num=="Donor")
Donor_1186   <- filter(pair_plot_df,Pair=="1186" & Case_num=="Donor")
Donor_1208   <- filter(pair_plot_df,Pair=="1208" & Case_num=="Donor")
Donor_1220   <- filter(pair_plot_df,Pair=="1220" & Case_num=="Donor")
Donor_1251   <- filter(pair_plot_df,Pair=="1251" & Case_num=="Donor")
Donor_1585   <- filter(pair_plot_df,Pair=="1585" & Case_num=="Donor")
Donor_1621   <- filter(pair_plot_df,Pair=="1621" & Case_num=="Donor")
Donor_1641   <- filter(pair_plot_df,Pair=="1641" & Case_num=="Donor")
Donor_1673   <- filter(pair_plot_df,Pair=="1673" & Case_num=="Donor")
Donor_1696   <- filter(pair_plot_df,Pair=="1696" & Case_num=="Donor")
Donor_1722_1 <- filter(pair_plot_df,Pair=="1722" & Case_num=="Donor")
Donor_1722_2 <- filter(pair_plot_df,Pair=="1722" & Case_num=="Donor")
Donor_1747   <- filter(pair_plot_df,Pair=="1747" & Case_num=="Donor")
Donor_1785   <- filter(pair_plot_df,Pair=="1785" & Case_num=="Donor")

Donor_643$line    <- 2
Donor_941_1$line  <- 1
Donor_942_1$line  <- 1
Donor_948$line    <- 2
Donor_951$line    <- 2
Donor_1002_1$line <- 2
Donor_1002_2$line <- 3
Donor_1031$line   <- 2
Donor_1033_1$line <- 1
Donor_1033_2$line <- 2
Donor_1048$line   <- 2
Donor_1162_1$line <- 2
Donor_1162_2$line <- 3
Donor_1186$line   <- 2
Donor_1208$line   <- 2
Donor_1220$line   <- 2
Donor_1251$line   <- 2
Donor_1585$line   <- 1
Donor_1621$line   <- 2
Donor_1641$line   <- 2
Donor_1673$line   <- 2
Donor_1696$line   <- 2
Donor_1722_1$line <- 2
Donor_1722_2$line <- 3
Donor_1747$line   <- 2
Donor_1785$line   <- 2

pair_plot_df <- rbind(Donor_643,
                      Donor_941_1,
                      Donor_941_2,
                      Donor_942_1,
                      Donor_942_2,
                      Donor_948,
                      Donor_951,
                      Donor_1002_1,
                      Donor_1002_2,
                      Donor_1031,
                      Donor_1033_1,
                      Donor_1033_2,
                      Donor_1048,
                      Donor_1162_1,
                      Donor_1162_2,
                      Donor_1186,
                      Donor_1208,
                      Donor_1220,
                      Donor_1251,
                      Donor_1585,
                      Donor_1621,
                      Donor_1641,
                      Donor_1673,
                      Donor_1696,
                      Donor_1722_1,
                      Donor_1722_2,
                      Donor_1747,
                      Donor_1785,
                      pair_plot_df)

Donor_643    <- NULL
Donor_941_1  <- NULL
Donor_941_2  <- NULL
Donor_942_1  <- NULL
Donor_942_2  <- NULL
Donor_948    <- NULL
Donor_951    <- NULL
Donor_1002_1 <- NULL
Donor_1002_2 <- NULL
Donor_1031   <- NULL
Donor_1033_1 <- NULL
Donor_1033_2 <- NULL
Donor_1048   <- NULL
Donor_1162_1 <- NULL
Donor_1162_2 <- NULL
Donor_1186   <- NULL
Donor_1208   <- NULL
Donor_1220   <- NULL
Donor_1251   <- NULL
Donor_1585   <- NULL
Donor_1621   <- NULL
Donor_1641   <- NULL
Donor_1673   <- NULL
Donor_1696   <- NULL
Donor_1722_1 <- NULL
Donor_1722_2 <- NULL
Donor_1747   <- NULL
Donor_1785   <- NULL

# change facet order
pair_plot_df$Pair <- as.factor(as.character(pair_plot_df$Pair))
pair_plot_df$Pair <- factor(pair_plot_df$Pair, 
                            levels=c("1162","1722",
                                     "1641","1696",
                                     "1179","1205","1536","1585","1636","1656","1705","1760","1773","1780","1787","1825","1844",
                                     "1002","1031","1033","1048","1186","1208","1220",
                                     
                                     "1021","1025",
                                     "1081","1137","1147",
                                     "1574",
                                     "1619","1621",
                                     "1743",
                                     "1761","1769",
                                     "1808","1836","1850",
                                     "922","941","942","970","977","998",
                                     "944","1747","1699","1673","1655"))


## rename
df_samplelist_cleaned <- sampleList
df_coinfection <- coinfList
df_pair <- pair_df
df_pair_plot <- pair_plot_df
sampleList <- NULL
coinfList <- NULL
pair_df <- NULL
pair_plot_df <- NULL

#### samplelist - Clean time columns ####
## convert DoC to year and month
df_samplelist_cleaned$DoC <- ymd(df_samplelist_cleaned$DoC)
df_samplelist_cleaned$DoC_Y <- year(df_samplelist_cleaned$DoC)
df_samplelist_cleaned$DoC_M <- month(df_samplelist_cleaned$DoC)
df_samplelist_cleaned$DoC_D <- day(df_samplelist_cleaned$DoC)
df_samplelist_cleaned$DoC_MY <- paste(df_samplelist_cleaned$DoC_M, "-",df_samplelist_cleaned$DoC_Y)
df_samplelist_cleaned$DoC_MY <- gsub(" ", "", df_samplelist_cleaned$DoC_MY)

## Season column
df_samplelist_cleaned$Season <- ""
df_samplelist_cleaned$Season[df_samplelist_cleaned$DoC_Y==2016 & df_samplelist_cleaned$DoC_M>6] <- "16-17"
df_samplelist_cleaned$Season[df_samplelist_cleaned$DoC_Y==2017 & df_samplelist_cleaned$DoC_M<6] <- "16-17"

df_samplelist_cleaned$Season[df_samplelist_cleaned$DoC_Y==2017 & df_samplelist_cleaned$DoC_M>6] <- "17-18"
df_samplelist_cleaned$Season[df_samplelist_cleaned$DoC_Y==2018 & df_samplelist_cleaned$DoC_M<6] <- "17-18"

df_samplelist_cleaned$Season[df_samplelist_cleaned$DoC_Y==2018 & df_samplelist_cleaned$DoC_M>6] <- "18-19"
df_samplelist_cleaned$Season[df_samplelist_cleaned$DoC_Y==2019 & df_samplelist_cleaned$DoC_M<6] <- "18-19"

## reaorder data
df_samplelist_cleaned <- df_samplelist_cleaned[,c(8,1,2,3,13,4,5,6,7,9,10,11,12)]


#### sampleList - Plot df_pair_plot ####
# plot
plot_pairs <- ggplot(df_pair_plot[!is.na(df_pair_plot$Pair),], 
                     aes(x=factor(Case, level=c("Donor","Recipient")), 
                                       y=value, 
                                       group=Pair, color=Type_Subtype)) + 
  geom_point(size = 0.75, aes(shape=Case)) + 
  geom_line(aes(group=line)) +
  coord_cartesian(ylim = c(-0.5,3.5)) +
  scale_shape_manual(values=c(16,21)) +
  scale_color_manual(
    values = c(`Influenza A`               = "#A26102", # Dark Orange
               `Influenza A (H1N1)`        = "#D28E2B", # Mid Orange
               `Influenza A (H3N2)`        = "#EEC588", # Light Orange
               `Influenza B`               = "#004689", # Dark Blue
               `Influenza B (Victoria)`    = "#1970C3", # Mid Blue
               `Influenza B (Yamagata)`    = "#5DA7EC")) + # Light Blue
  guides(color=guide_legend(override.aes=list(shape=15, size=6)),
         shape=guide_legend(override.aes=list(size=4))) + 
  facet_wrap(~Pair,ncol=10) + 
  labs(x = "", y = "") + 
  axis_formatting + legend_formatting + background_formatting + 
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size=4, vjust=-4),
        panel.border = element_rect(color="grey",fill=NA,size=.5),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(-.5, "lines"),
        panel.background = element_blank())
plot_pairs

setwd("/Users/hries/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/figures")
ggsave("TransmissionPairs.pdf", plot_pairs,
       width = 5, height = 2.5,
       units = "in", dpi = 320)

# legend
plot_pairs_legend <- ggplot(df_pair_plot, aes(x=factor(Case, level=c("Donor","Recipient")), 
                                              y=value, 
                                              group=Pair, color=Type_Subtype)) + 
  geom_point(size = 1, aes(shape=Case)) + 
  geom_line(aes(group=line)) +
  coord_cartesian(ylim = c(-0.5,3.5)) +   
  scale_shape_manual(values=c(16,21)) +
  scale_color_manual(
    values = c(`Influenza A`               = "#A26102", # Dark Orange
               `Influenza A (H1N1)`        = "#D28E2B", # Mid Orange
               `Influenza A (H3N2)`        = "#EEC588", # Light Orange
               `Influenza B`               = "#004689", # Dark Blue
               `Influenza B (Victoria)`    = "#1970C3", # Mid Blue
               `Influenza B (Yamagata)`    = "#5DA7EC")) + # Light Blue
  guides(color=guide_legend(override.aes=list(shape=15, size=6)),
         shape=guide_legend(override.aes=list(size=4))) + 
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text= element_text(size=4, vjust=-4),
        panel.border = element_rect(color="grey",fill=NA,size=.5),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(-.5, "lines"),
        panel.background = element_blank()) +
  facet_wrap(~Pair,ncol=12)

#### sampleList - Plot df_samplelist_cleaned ####
ct_time <- ggplot(df_samplelist_cleaned, aes(x=DoC, y=Ct, group=Type, color=Type_Subtype)) + 
  geom_point(size = 0.5) + 
  coord_cartesian(
    ylim = c(15, 40)) + 
  scale_y_continuous(
    breaks = c(15, 20, 25, 30, 35, 40),
    labels = c(15, 20, 25, 30, 35, 40),
    expand= c(0,0)) +
  scale_x_date(
    date_labels="%b %Y", 
    breaks=as.Date(c("2017-07-01","2017-10-01",
                     "2018-01-01","2018-04-01","2018-07-01","2018-10-01",
                     "2019-01-01","2019-04-01","2019-07-01")),
    limits=as.Date(c("2017-07-01","2019-07-01"))) +
  scale_color_manual(
    values = c(`Influenza A`               = "#A26102", # Dark Orange
               `Influenza A (H1N1)`        = "#D28E2B", # Mid Orange
               `Influenza A (H3N2)`        = "#EEC588", # Light Orange
               `Influenza B`               = "#004689", # Dark Blue
               `Influenza B (Victoria)`    = "#1970C3", # Mid Blue
               `Influenza B (Yamagata)`    = "#5DA7EC")) + # Light Blue
  labs(y = "Ct value", title = "") + 
  theme_bw() + #removes grey background
  #geom_hline(yintercept = 33, linetype = "dashed", color = "grey") + 
  #annotate("text", x = ymd("2018-07-01"), y = 38, 
  #         label = paste(length(df_samplelist_cleaned[df_samplelist_cleaned$Ct>33,]), 
  #                       " of 284 samples\nwhere Ct>33", sep = "")) + 
  theme(
    panel.grid.major = element_blank(),       #removes ver. grid lines
    #panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.title = element_blank(),           #remove legend title
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 5),     #give us some room!
    axis.title.x = element_blank(),           #remove x-axis title
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust=1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size

hist <- ggplot(df_samplelist_cleaned, aes(x=DoC, fill=Type_Subtype)) + 
  geom_histogram(binwidth=14) + 
  #geom_hline(yintercept=0.3) + 
  coord_cartesian(
    ylim = c(0, 50)) + 
  scale_y_continuous(
    breaks = seq(0, 50, by=10),
    labels = seq(0, 50, by=10),
    expand= c(0,0)) +
  scale_x_date(date_labels="%b %Y", 
               breaks=as.Date(c("2017-07-01","2017-10-01",
                                "2018-01-01","2018-04-01","2018-07-01","2018-10-01",
                                "2019-01-01","2019-04-01","2019-07-01")),
               limits=as.Date(c("2017-07-01","2019-07-01"))) +
  scale_fill_manual(
    values = c(`Influenza A`               = "#A26102", # Dark Orange
               `Influenza A (H1N1)`        = "#D28E2B", # Mid Orange
               `Influenza A (H3N2)`        = "#EEC588", # Light Orange
               `Influenza B`               = "#004689", # Dark Blue
               `Influenza B (Victoria)`    = "#1970C3", # Mid Blue
               `Influenza B (Yamagata)`    = "#5DA7EC")) + # Light Blue
  labs(y = "Count") + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),       #removes ver. grid lines
    #panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.title = element_blank(),           #remove legend title
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 5),     #give us some room!
    axis.title.x = element_blank(),           #remove x-axis title
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust=1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size

#### sampleList - Supplemental 1 ####
plot_dot_legend <- ct_time + theme(legend.position="right")

legend_df_dot_plot <- get_legend(plot_dot_legend)
legend_df_dot_plot <- as_ggplot(legend_df_dot_plot)

plot_hist_title <- hist + labs(title="Samples that passed QC") + 
  theme(plot.title = element_text(vjust=-.5))
plot_dot_title <- ct_time + labs(title="Samples that passed QC") + 
  theme(plot.title = element_text(vjust=-.5))

plot_supp1 <- plot_grid(plot_hist_title,
                        plot_dot_title,
                        legend_df_dot_plot,
                        ncol=3, nrow=1, 
                        rel_heights = c(1,1,1),
                        rel_widths = c(1,1,.75),
                        labels = c("A","B",""), 
                        label_size=12)

#### Metadata - clean ####
#metadata_backup <- metadata
#metadata <- metadata_backup

## only keep samples that passed QC
sampleList <- sampleList[sampleList$WSLH_ID %in% samples_QC_pass$WSLH_ID,]

## Num column
metadata$Num <- paste(metadata$record_id)
metadata$Num <- gsub("-1","",metadata$Num)
metadata$Num <- gsub("-2","",metadata$Num)
metadata$Num <- gsub("-3","",metadata$Num)
metadata$Num <- gsub("-4","",metadata$Num)
metadata$Num <- gsub("-5","",metadata$Num)
metadata$Num <- gsub("-6","",metadata$Num)
metadata$Num <- gsub("-7","",metadata$Num)
metadata$Num <- as.factor(metadata$Num)

## Add Case_Num column
metadata$Case_Num <- paste(metadata$Case, metadata$Num)
metadata$Case_Num <- as.factor(metadata$Case_Num)
metadata$Case <- as.factor(metadata$Case)

## Clean age column
metadata$age <- gsub("left blank",NA,metadata$age)
metadata$age <- gsub("0 \\(4 weeks)",".08",metadata$age)
metadata$age <- gsub("15 mo.","1.25",metadata$age)
metadata$age <- gsub("5 months",".42",metadata$age)
metadata$age <- gsub("20 months","1.67",metadata$age)
metadata$age <- gsub("6 days","0.02",metadata$age)
metadata$age <- gsub("1 month",".08",metadata$age)
metadata$age <- gsub("16 months","1.33",metadata$age)
metadata$age <- gsub("21 months","1.75",metadata$age)
metadata$age <- gsub("8 months",".67",metadata$age)
metadata$age <- gsub("10 months",".83",metadata$age)
metadata$age <- gsub("18 mo.","1.5",metadata$age)
metadata$age <- gsub("9 months",".75",metadata$age)
metadata$age <- gsub("18 months","1.5",metadata$age)
metadata$age <- gsub("18 Months","1.5",metadata$age)
metadata$age <- gsub("5 mo.",".42",metadata$age)
metadata$age <- gsub("7 months",".58",metadata$age)
metadata$age <- gsub("2.08s","2.08",metadata$age)
metadata$age <- as.numeric(metadata$age)

## Clean gender golumn
metadata$gender <- gsub("1","Male",metadata$gender)
metadata$gender <- gsub("2","Female",metadata$gender)
metadata$gender <- as.factor(metadata$gender)

## Clean flu_vaccine
metadata$flu_vaccine <- gsub("1","Yes", metadata$flu_vaccine)
metadata$flu_vaccine <- gsub("0","No", metadata$flu_vaccine)
metadata$flu_vaccine <- as.factor(metadata$flu_vaccine); levels(metadata$flu_vaccine)

## Clean days onset
metadata$days_between_onset_visit <- as.factor(metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("\\?","NA",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("n/a","NA",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("constant - allergies","NA",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("ongoing","NA",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("0 \\(today)","NA",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("0.21","0",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("0.5","0",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("21ish","21",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("a few days ago","3",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("<12 hours ago","0",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("<1","0",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("7+","7",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("14-21 days","14",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("2-3 weeks","14",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("7-14","7",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("7-14+","7",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("long time \\(6 months)","28",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("NA",NA,metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- as.factor(metadata$days_between_onset_visit)

## days to timeframes
metadata$days_between_onset_visit <- as.numeric(metadata$days_between_onset_visit)
metadata$days_between_onset_visit[metadata$days_between_onset_visit<7] <- 0
metadata$days_between_onset_visit[metadata$days_between_onset_visit<14 & metadata$days_between_onset_visit>0] <- 7
metadata$days_between_onset_visit[metadata$days_between_onset_visit<21 & metadata$days_between_onset_visit>=14] <- 14
metadata$days_between_onset_visit[metadata$days_between_onset_visit<28 & metadata$days_between_onset_visit>=21] <- 21
metadata$days_between_onset_visit[metadata$days_between_onset_visit>=28] <- 28
metadata$days_between_onset_visit <- as.factor(metadata$days_between_onset_visit)

## timeframes to days again
metadata$days_between_onset_visit <- as.character(metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("0","0-1 week",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("7","1-2 week",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("14","2-3 week",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("21","3-4 week",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("28","4 plus",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub(" week","",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("0-1","0-7 days",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("1-2","7-14 days",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("2-3","14-21 days",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("3-4","21-28 days",metadata$days_between_onset_visit)
metadata$days_between_onset_visit <- gsub("4 plus","28+ days",metadata$days_between_onset_visit)
metadata <- dplyr::rename(metadata,days_symptom_onset=`days_between_onset_visit`)

## clean home_visit_appointment
metadata$home_visit_appointment <- mdy_hm(metadata$home_visit_appointment)
metadata$day0_date_household <- NULL
metadata$day7_date_household <- NULL

## rearrange columns
metadata <- metadata[,c(1,8,2,9,3,4,5,6,7)]

#### Merge data frames ####
## merge by shared Study_ID or record_id
df_samplelist_metadata <- merge(df_samplelist_cleaned, metadata,
                                all.x = TRUE,
                                by.x = "Study_ID",
                                by.y = "record_id")
df_samplelist_metadata$Sample <- paste("HA","|",df_samplelist_metadata$Num,"|",df_samplelist_metadata$WSLH_ID)
df_samplelist_metadata$Sample <- gsub(" ","",df_samplelist_metadata$Sample)

df_pair_metadata <- merge(df_pair, metadata,
                          all.x = TRUE,
                          by.x = "Study_ID",
                          by.y = "record_id")
df_pair_metadata$Sample <- paste("HA","|",df_pair_metadata$Num,"|",df_pair_metadata$WSLH_ID)
df_pair_metadata$Sample <- gsub(" ","",df_pair_metadata$Sample)

df_coinfection_metadata <- merge(df_coinfection, metadata,
                                 all.x = TRUE,
                                 by.x = "Study_ID",
                                 by.y = "record_id")
df_coinfection_metadata$Sample <- paste("HA","|",df_coinfection_metadata$Num,"|",df_coinfection_metadata$WSLH_ID)
df_coinfection_metadata$Sample <- gsub(" ","",df_coinfection_metadata$Sample)

#### df_samplelist_metadata - Plot age pyramid (all) ####
df_pyramid_plot <- df_samplelist_metadata
df_pyramid_plot$age <- as.integer(df_pyramid_plot$age)
df_pyramid_plot$gender <- factor(df_pyramid_plot$gender, levels=c("Male", "Female"), labels=c("Male", "Female"))

plot_pyramid_all <- ggplot(data=df_pyramid_plot,aes(x=age,fill=gender)) + 
  geom_bar(data=subset(df_pyramid_plot,gender=="Female"),aes(y=..count..*(1))) + 
  geom_bar(data=subset(df_pyramid_plot,gender=="Male"),aes(y=..count..*(-1))) + 
  scale_y_continuous(breaks=seq(-45,45,10),labels=abs(seq(-45,45,10)),limits=c(-45,45)) + 
  scale_x_continuous(breaks=seq(-60,60,5),labels=abs(seq(-60,60,5)),limits=c(0,60)) + 
  scale_fill_manual("", values = c("Female" = "#004D40", "Male" = "#FFC107")) + 
  labs(
    y = "Count",
    x = "Age",
    title = "ORCHARDS Age Pyramid - All Samples") +  
  coord_flip() + 
  theme_bw() + #removes grey background
  theme(
    #panel.grid.major = element_blank(),       #removes ver. grid lines
    panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 20),     #give us some room!
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size

#### df_samplelist_metadata - Plot age pyramid (pair) ####
df_pyramid_plot <- df_pair_metadata
df_pyramid_plot$age <- as.integer(df_pyramid_plot$age)
df_pyramid_plot$gender <- factor(df_pyramid_plot$gender, levels=c("Male", "Female"), labels=c("Male", "Female"))

plot_pyramid_pair <- ggplot(data=df_pyramid_plot,aes(x=age,fill=gender)) + 
  geom_bar(data=subset(df_pyramid_plot,gender=="Female"),aes(y=..count..*(1))) + 
  geom_bar(data=subset(df_pyramid_plot,gender=="Male"),aes(y=..count..*(-1))) + 
  scale_y_continuous(breaks=seq(-45,45,10),labels=abs(seq(-45,45,10)),limits=c(-45,45)) + 
  scale_x_continuous(breaks=seq(-60,60,5),labels=abs(seq(-60,60,5)),limits=c(0,60)) + 
  scale_fill_manual("", values = c("Female" = "#004D40", "Male" = "#FFC107")) + 
  labs(
    y = "Count",
    x = "Age",
    title = "ORCHARDS Age Pyramid - Paired Samples") +  
  coord_flip() + 
  theme_bw() + #removes grey background
  theme(
    #panel.grid.major = element_blank(),       #removes ver. grid lines
    panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 20),     #give us some room!
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size

#### df_samplelist_metadata - Supplemental 3 ####
plot_pyramid_legend <- plot_pyramid_all
plot_pyramid_legend <- plot_pyramid_legend + theme(legend.position="right")

legend_df_samplelist_metadata <- get_legend(plot_pyramid_legend)
legend_df_samplelist_metadata <- as_ggplot(legend_df_samplelist_metadata)

plot_pyramid_all_title <- plot_pyramid_all + labs(title="All Community Samples") + 
  theme(plot.title = element_text(vjust=-.5))
plot_pyramid_pair_title <- plot_pyramid_pair + labs(title="Paired Samples") + 
  theme(plot.title = element_text(vjust=-.5))

plot_supp3 <- plot_grid(plot_pyramid_all_title,
                        plot_pyramid_pair_title,
                        legend_df_samplelist_metadata,
                        ncol=3, nrow=1, 
                        rel_heights = c(1,1,1),
                        rel_widths = c(1,1,.33),
                        labels = c("A","B",""), 
                        label_size=12)

ggsave("Supplemental_3.pdf", plot_supp3,
       dpi = 300, unit = c("in"),
       width = 7.5, height = 3)

#### df_samplelist_metadata - Plot days since symptom onset ####
df_symptom_plot <- df_samplelist_metadata
df_symptom_plot$days_symptom_onset <- as.factor(df_symptom_plot$days_symptom_onset)

plot_symptom_all <- ggplot(data=df_symptom_plot,aes(x=days_symptom_onset,fill=days_symptom_onset)) + 
  geom_bar(data=df_symptom_plot,aes(y=..count..)) + 
  scale_y_continuous(breaks=seq(0,150,25),labels=abs(seq(0,150,25)),limits=c(0,150)) + 
  scale_fill_manual("", values = c("0-7 days" = "#1E88E5",
                                   "7-14 days" = "#004D40", 
                                   "14-21 days" = "#FFC107",
                                   "21-28 days" = "#610882",
                                   "28+ days" = "#D81B60")) + 
  scale_x_discrete(limits = c("0-7 days","7-14 days",
                              "14-21 days","21-28 days",
                              "28+ days")) + 
  labs(
    y = "Count",
    x = "Days since symptom onset") +  
  coord_flip() + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),       #removes ver. grid lines
    panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 20),     #give us some room!
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size


df_symptom_plot <- df_pair_metadata
df_symptom_plot$days_symptom_onset <- as.factor(df_symptom_plot$days_symptom_onset)

plot_vax_pair <- ggplot(data=df_symptom_plot,aes(x=days_symptom_onset,fill=days_symptom_onset)) + 
  geom_bar(data=df_symptom_plot,aes(y=..count..)) + 
  scale_y_continuous(breaks=seq(0,150,25),labels=abs(seq(0,150,25)),limits=c(0,150)) + 
  scale_fill_manual("", values = c("0-7 days" = "#1E88E5",
                                   "7-14 days" = "#004D40", 
                                   "14-21 days" = "#FFC107",
                                   "21-28 days" = "#610882",
                                   "28+ days" = "#D81B60")) + 
  scale_x_discrete(limits = c("0-7 days","7-14 days",
                              "14-21 days","21-28 days",
                              "28+ days")) + 
  labs(
    y = "Count",
    x = "Days since symptom onset") +  
  coord_flip() + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),       #removes ver. grid lines
    panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 20),     #give us some room!
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size

#### df_samplelist_metadata - Supplemental 4 ####
plot_vax_all_title <- plot_vax_all + labs(title="All Community Samples") + 
  theme(plot.title = element_text(vjust=-.5))
plot_vax_pair_title <- plot_vax_pair + labs(title="Paired Samples") + 
  theme(plot.title = element_text(vjust=-.5))

plot_supp4 <- plot_grid(plot_vax_all_title,
                        plot_vax_pair_title,
                        ncol=2, nrow=1, 
                        rel_heights = c(1,1),
                        rel_widths = c(1,1),
                        labels = c("A","B"), 
                        label_size=12)

ggsave("Supplemental_4.pdf", plot_supp4,
       dpi = 300, unit = c("in"),
       width = 7.5, height = 3)

#### df_samplelist_metadata - Plot vaccination status ####
df_vax_plot <- df_samplelist_metadata
df_vax_plot$flu_vaccine <- as.character(df_vax_plot$flu_vaccine)

plot_vax_all <- ggplot(data=df_vax_plot,aes(x=flu_vaccine,fill=flu_vaccine)) + 
  geom_bar(data=df_vax_plot,aes(y=..count..)) + 
  #scale_y_continuous(breaks=seq(0,150,25),labels=abs(seq(0,150,25)),limits=c(0,150)) + 
  #scale_fill_manual("", values = c("0-7 days" = "#1E88E5",
  #                                 "7-14 days" = "#004D40", 
  #                                 "14-21 days" = "#FFC107",
  #                                 "21-28 days" = "#610882",
  #                                 "28+ days" = "#D81B60")) + 
  #scale_x_discrete(limits = c("0-7 days","7-14 days",
  #                            "14-21 days","21-28 days",
  #                            "28+ days")) + 
  labs(
    y = "Count",
    x = "Vaccination status") +  
  coord_flip() + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),       #removes ver. grid lines
    panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 20),     #give us some room!
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size


df_vax_plot <- df_pair_metadata
df_vax_plot$days_symptom_onset <- as.factor(df_vax_plot$days_symptom_onset)

plot_vax_pair <- ggplot(data=df_vax_plot,aes(x=days_symptom_onset,fill=days_symptom_onset)) + 
  geom_bar(data=df_vax_plot,aes(y=..count..)) + 
  scale_y_continuous(breaks=seq(0,150,25),labels=abs(seq(0,150,25)),limits=c(0,150)) + 
  scale_fill_manual("", values = c("0-7 days" = "#1E88E5",
                                   "7-14 days" = "#004D40", 
                                   "14-21 days" = "#FFC107",
                                   "21-28 days" = "#610882",
                                   "28+ days" = "#D81B60")) + 
  scale_x_discrete(limits = c("0-7 days","7-14 days",
                              "14-21 days","21-28 days",
                              "28+ days")) + 
  labs(
    y = "Count",
    x = "Days since symptom onset") +  
  coord_flip() + 
  theme_bw() + #removes grey background
  theme(
    panel.grid.major = element_blank(),       #removes ver. grid lines
    panel.grid.minor = element_blank(),       #removes hor. grid lines
    legend.position = "none",                 #place legend on right
    legend.text = element_text(size = 10),    #change text size
    plot.title = element_text(hjust = 0),     #move title left
    plot.margin = margin(10, 10, 10, 20),     #give us some room!
    strip.text.x = element_text(size = 6),    #change text size
    strip.text.y = element_text(size = 6))    #change text size


#### df_samplelist_metadata - Supplemental 5 ####
plot_pyramid_legend <- plot_pyramid_all
plot_pyramid_legend <- plot_pyramid_legend + theme(legend.position="right")

legend_df_samplelist_metadata <- get_legend(plot_pyramid_legend)
legend_df_samplelist_metadata <- as_ggplot(legend_df_samplelist_metadata)

plot_pyramid_all_title <- plot_pyramid_all + labs(title="All Community Samples") + 
  theme(plot.title = element_text(vjust=-.5))
plot_pyramid_pair_title <- plot_pyramid_pair + labs(title="Paired Samples") + 
  theme(plot.title = element_text(vjust=-.5))

plot_supp5 <- plot_grid(plot_pyramid_all_title,
                        plot_pyramid_pair_title,
                        legend_df_samplelist_metadata,
                        ncol=3, nrow=1, 
                        rel_heights = c(1,1,1),
                        rel_widths = c(1,1,.33),
                        labels = c("A","B",""), 
                        label_size=12)

ggsave("Supplemental_5.pdf", plot_supp5,
       dpi = 300, unit = c("in"),
       width = 7.5, height = 3)

#### Save data frames ####
#write.csv(x=df_coinfection_metadata, file="df_coinfection_metadata.csv")
#write.csv(x=df_pair_metadata, file="df_pair_metadata.csv")
#write.csv(x=df_samplelist_metadata, file="df_samplelist_metadata.csv")
#write.csv(x=metadata, file="df_metadata.csv")

#### Save plots ####
setwd(dir_s); dir()

## Supp 1: 
ggsave("Supplemental_1.pdf", plot_supp1, 
       dpi = 300, unit = c("in"),
       width = 7.5, height = 3.5)

## Supp 2: Transmission pairs
ggsave("Supplemental_2.pdf", plot_pairs, 
       dpi = 300, unit = c("in"),
       width = 5, height = 4)
