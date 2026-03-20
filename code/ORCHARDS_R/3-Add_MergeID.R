## Clear Global Environment
rm(list = ls())

## Directories

dir_1718_H3N2 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H3N2/dfs",sep="")
dir_1819_H3N2 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H3N2/dfs",sep="")
dir_1718_H1N1 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/17-18_H1N1/dfs",sep="")
dir_1819_H1N1 <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/reads_3/220208_ORCHARDS_Share/reads/18-19_H1N1/dfs",sep="")
dir_s         <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/dfs",sep="")

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

#### Import 17-18_H3N2 ####
## 17-18_H3N2
setwd(dir_1718_H3N2); getwd(); dir()

df_all_iSNVs_1                        <- read.csv("df_all_iSNVs.csv")
df_sub50AF_1                          <- read.csv("df_sub50AF.csv")
df_sub50AF_md_n_SNVs_by_subtype_1     <- read.csv("df_sub50AF_md_n_SNVs_by_subtype.csv")
df_sub50AF_md_n_SNVs_by_vax_1         <- read.csv("df_sub50AF_md_n_SNVs_by_vax.csv")
df_mut_bins_prop_noNeut_1             <- read.csv("df_mut_bins_prop_noNeut.csv")
df_compiled_1                         <- read.csv("df_compiled.csv")
df_divergence_1                       <- read.csv("df_divergence.csv")
ds_divergence_1                       <- read.csv("ds_divergence.csv")
df_prop_shared_1                      <- read.csv("df_prop_shared.csv")
df_JTplot_1                           <- read.csv("df_JTplot.csv")
df_household_1                        <- read.csv("df_household.csv")
sg_1                                  <- read.csv("sg.csv")
sg_sw_1                               <- read.csv("sg_sw.csv")
sg_paired_1                           <- read.csv("sg_paired.csv")
df_mut_bins_prop_1                    <- read.csv("df_mut_bins_prop.csv")
ds_sg_sw_1                            <- read.csv("ds_sg_sw.csv")
df_sub50AF_md_n_SNVs_by_bin_and_mut_1 <- read.csv("df_sub50AF_md_n_SNVs_by_bin_and_mut.csv")
df_sg_AA_HA_1                         <- read.csv("df_sg_AA_HA.csv")
df_AA_HA_1                            <- read.csv("df_AA_HA.csv")
df_AA_HA_paired_1                     <- read.csv("df_AA_HA_paired.csv")
df_mut_bins_prop_all_1                <- read.csv("df_mut_bins_prop_all.csv")

## Merge_ID function
Add_Merge_ID_17_18_H3N2 <- function(data) {
  data$Merge_ID <- "17-18_H3N2"
  return(data)
}

## Execute
df_all_iSNVs_1                        <- Add_Merge_ID_17_18_H3N2(df_all_iSNVs_1)
df_sub50AF_1                          <- Add_Merge_ID_17_18_H3N2(df_sub50AF_1)
df_sub50AF_md_n_SNVs_by_subtype_1     <- Add_Merge_ID_17_18_H3N2(df_sub50AF_md_n_SNVs_by_subtype_1)
df_sub50AF_md_n_SNVs_by_vax_1         <- Add_Merge_ID_17_18_H3N2(df_sub50AF_md_n_SNVs_by_vax_1)
df_mut_bins_prop_noNeut_1             <- Add_Merge_ID_17_18_H3N2(df_mut_bins_prop_noNeut_1)
df_compiled_1                         <- Add_Merge_ID_17_18_H3N2(df_compiled_1)
df_divergence_1                       <- Add_Merge_ID_17_18_H3N2(df_divergence_1)
ds_divergence_1                       <- Add_Merge_ID_17_18_H3N2(ds_divergence_1)
df_prop_shared_1                      <- Add_Merge_ID_17_18_H3N2(df_prop_shared_1)
df_JTplot_1                           <- Add_Merge_ID_17_18_H3N2(df_JTplot_1)
df_household_1                        <- Add_Merge_ID_17_18_H3N2(df_household_1)
sg_1                                  <- Add_Merge_ID_17_18_H3N2(sg_1)
sg_sw_1                               <- Add_Merge_ID_17_18_H3N2(sg_sw_1)
sg_paired_1                           <- Add_Merge_ID_17_18_H3N2(sg_paired_1)
df_mut_bins_prop_1                    <- Add_Merge_ID_17_18_H3N2(df_mut_bins_prop_1)
ds_sg_sw_1                            <- Add_Merge_ID_17_18_H3N2(ds_sg_sw_1)
df_sub50AF_md_n_SNVs_by_bin_and_mut_1 <- Add_Merge_ID_17_18_H3N2(df_sub50AF_md_n_SNVs_by_bin_and_mut_1)
df_sg_AA_HA_1                         <- Add_Merge_ID_17_18_H3N2(df_sg_AA_HA_1)
df_AA_HA_1                            <- Add_Merge_ID_17_18_H3N2(df_AA_HA_1)
df_AA_HA_paired_1                     <- Add_Merge_ID_17_18_H3N2(df_AA_HA_paired_1)
df_mut_bins_prop_all_1                <- Add_Merge_ID_17_18_H3N2(df_mut_bins_prop_all_1)

#### Import 18-19_H3N2 ####
## 18-19_H3N2
setwd(dir_1819_H3N2); getwd(); dir()

df_all_iSNVs_2                         <- read.csv("df_all_iSNVs.csv")
df_sub50AF_2                           <- read.csv("df_sub50AF.csv")
df_sub50AF_md_n_SNVs_by_subtype_2      <- read.csv("df_sub50AF_md_n_SNVs_by_subtype.csv")
df_sub50AF_md_n_SNVs_by_vax_2          <- read.csv("df_sub50AF_md_n_SNVs_by_vax.csv")
df_mut_bins_prop_noNeut_2              <- read.csv("df_mut_bins_prop_noNeut.csv")
df_compiled_2                          <- read.csv("df_compiled.csv")
df_divergence_2                        <- read.csv("df_divergence.csv")
ds_divergence_2                        <- read.csv("ds_divergence.csv")
df_prop_shared_2                       <- read.csv("df_prop_shared.csv")
df_JTplot_2                            <- read.csv("df_JTplot.csv")
df_household_2                         <- read.csv("df_household.csv")
sg_2                                   <- read.csv("sg.csv")
sg_sw_2                                <- read.csv("sg_sw.csv")
sg_paired_2                            <- read.csv("sg_paired.csv")
df_mut_bins_prop_2                     <- read.csv("df_mut_bins_prop.csv")
ds_sg_sw_2                             <- read.csv("ds_sg_sw.csv")
df_sub50AF_md_n_SNVs_by_bin_and_mut_2  <- read.csv("df_sub50AF_md_n_SNVs_by_bin_and_mut.csv")
df_sg_AA_HA_2                          <- read.csv("df_sg_AA_HA.csv")
df_AA_HA_2                             <- read.csv("df_AA_HA.csv")
df_AA_HA_paired_2                      <- read.csv("df_AA_HA_paired.csv")
df_mut_bins_prop_all_2                 <- read.csv("df_mut_bins_prop_all.csv")

## Merge_ID
Add_Merge_ID_18_19_H3N2 <- function(data) {
  data$Merge_ID <- "18-19_H3N2"
  return(data)
}

## Execute 
df_all_iSNVs_2                         <- Add_Merge_ID_18_19_H3N2(df_all_iSNVs_2)
df_sub50AF_2                           <- Add_Merge_ID_18_19_H3N2(df_sub50AF_2)
df_sub50AF_md_n_SNVs_by_subtype_2      <- Add_Merge_ID_18_19_H3N2(df_sub50AF_md_n_SNVs_by_subtype_2)
df_sub50AF_md_n_SNVs_by_vax_2          <- Add_Merge_ID_18_19_H3N2(df_sub50AF_md_n_SNVs_by_vax_2)
df_mut_bins_prop_noNeut_2              <- Add_Merge_ID_18_19_H3N2(df_mut_bins_prop_noNeut_2)
df_compiled_2                          <- Add_Merge_ID_18_19_H3N2(df_compiled_2)
df_divergence_2                        <- Add_Merge_ID_18_19_H3N2(df_divergence_2)
ds_divergence_2                        <- Add_Merge_ID_18_19_H3N2(ds_divergence_2)
df_prop_shared_2                       <- Add_Merge_ID_18_19_H3N2(df_prop_shared_2)
df_JTplot_2                            <- Add_Merge_ID_18_19_H3N2(df_JTplot_2)
df_household_2                         <- Add_Merge_ID_18_19_H3N2(df_household_2)
sg_2                                   <- Add_Merge_ID_18_19_H3N2(sg_2)
sg_sw_2                                <- Add_Merge_ID_18_19_H3N2(sg_sw_2)
sg_paired_2                            <- Add_Merge_ID_18_19_H3N2(sg_paired_2)
df_mut_bins_prop_2                     <- Add_Merge_ID_18_19_H3N2(df_mut_bins_prop_2)
ds_sg_sw_2                             <- Add_Merge_ID_18_19_H3N2(ds_sg_sw_2)
df_sub50AF_md_n_SNVs_by_bin_and_mut_2  <- Add_Merge_ID_18_19_H3N2(df_sub50AF_md_n_SNVs_by_bin_and_mut_2)
df_sg_AA_HA_2                          <- Add_Merge_ID_18_19_H3N2(df_sg_AA_HA_2)
df_AA_HA_2                             <- Add_Merge_ID_18_19_H3N2(df_AA_HA_2)
df_AA_HA_paired_2                      <- Add_Merge_ID_18_19_H3N2(df_AA_HA_paired_2)
df_mut_bins_prop_all_2                 <- Add_Merge_ID_18_19_H3N2(df_mut_bins_prop_all_2)

#### Merge H3N2 ####
df_all_iSNVs_H3N2                          <- rbind(df_all_iSNVs_1, df_all_iSNVs_2)
df_sub50AF_H3N2                            <- rbind(df_sub50AF_1, df_sub50AF_2)
df_sub50AF_md_n_SNVs_by_subtype_H3N2       <- rbind(df_sub50AF_md_n_SNVs_by_subtype_1, df_sub50AF_md_n_SNVs_by_subtype_2)
df_sub50AF_md_n_SNVs_by_vax_H3N2           <- rbind(df_sub50AF_md_n_SNVs_by_vax_1, df_sub50AF_md_n_SNVs_by_vax_2)
df_mut_bins_prop_noNeut_H3N2               <- rbind(df_mut_bins_prop_noNeut_1, df_mut_bins_prop_noNeut_2)
df_compiled_H3N2                           <- rbind(df_compiled_1, df_compiled_2)
df_divergence_H3N2                         <- rbind(df_divergence_1, df_divergence_2)
ds_divergence_H3N2                         <- rbind(ds_divergence_1, ds_divergence_2)
df_prop_shared_H3N2                        <- rbind(df_prop_shared_1, df_prop_shared_2)
df_JTplot_H3N2                             <- rbind(df_JTplot_1, df_JTplot_2)
df_household_H3N2                          <- rbind(df_household_1, df_household_2)
sg_H3N2                                    <- rbind(sg_1, sg_2)
sg_sw_H3N2                                 <- rbind(sg_sw_1, sg_sw_2)
sg_paired_H3N2                             <- rbind(sg_paired_1, sg_paired_2)
df_mut_bins_prop_H3N2                      <- rbind(df_mut_bins_prop_1, df_mut_bins_prop_2)
ds_sg_sw_H3N2                              <- rbind(ds_sg_sw_1, ds_sg_sw_2)
df_sub50AF_md_n_SNVs_by_bin_and_mut_H3N2   <- rbind(df_sub50AF_md_n_SNVs_by_bin_and_mut_1, df_sub50AF_md_n_SNVs_by_bin_and_mut_2)
df_sg_AA_HA_H3N2                           <- rbind(df_sg_AA_HA_1, df_sg_AA_HA_2)
df_AA_HA_H3N2                              <- rbind(df_AA_HA_1, df_AA_HA_2)
df_AA_HA_paired_H3N2                       <- rbind(df_AA_HA_paired_1, df_AA_HA_paired_2)
df_mut_bins_prop_all_H3N2                  <- rbind(df_mut_bins_prop_all_1, df_mut_bins_prop_all_2)

#### Import 17-18_H1N1 ####
## 17-18_H3N2
setwd(dir_1718_H1N1); getwd(); dir()

df_all_iSNVs_1                         <- read.csv("df_all_iSNVs.csv")
df_sub50AF_1                           <- read.csv("df_sub50AF.csv")
df_sub50AF_md_n_SNVs_by_subtype_1      <- read.csv("df_sub50AF_md_n_SNVs_by_subtype.csv")
df_sub50AF_md_n_SNVs_by_vax_1          <- read.csv("df_sub50AF_md_n_SNVs_by_vax.csv")
df_mut_bins_prop_noNeut_1              <- read.csv("df_mut_bins_prop_noNeut.csv")
df_compiled_1                          <- read.csv("df_compiled.csv")
df_divergence_1                        <- read.csv("df_divergence.csv")
ds_divergence_1                        <- read.csv("ds_divergence.csv")
df_prop_shared_1                       <- read.csv("df_prop_shared.csv")
df_JTplot_1                            <- read.csv("df_JTplot.csv")
df_household_1                         <- read.csv("df_household.csv")
sg_1                                   <- read.csv("sg.csv")
sg_sw_1                                <- read.csv("sg_sw.csv")
sg_paired_1                            <- read.csv("sg_paired.csv")
df_mut_bins_prop_1                     <- read.csv("df_mut_bins_prop.csv")
ds_sg_sw_1                             <- read.csv("ds_sg_sw.csv")
df_sub50AF_md_n_SNVs_by_bin_and_mut_1  <- read.csv("df_sub50AF_md_n_SNVs_by_bin_and_mut.csv")
df_sg_AA_HA_1                          <- read.csv("df_sg_AA_HA.csv")
df_AA_HA_1                             <- read.csv("df_AA_HA.csv")
df_AA_HA_paired_1                      <- read.csv("df_AA_HA_paired.csv")
df_mut_bins_prop_all_1                 <- read.csv("df_mut_bins_prop_all.csv")

## Merge_ID function
Add_Merge_ID_17_18_H1N1 <- function(data) {
  data$Merge_ID <- "17-18_H1N1"
  return(data)
}

## Execute
df_all_iSNVs_1                         <- Add_Merge_ID_17_18_H1N1(df_all_iSNVs_1)
df_sub50AF_1                           <- Add_Merge_ID_17_18_H1N1(df_sub50AF_1)
df_sub50AF_md_n_SNVs_by_subtype_1      <- Add_Merge_ID_17_18_H1N1(df_sub50AF_md_n_SNVs_by_subtype_1)
df_sub50AF_md_n_SNVs_by_vax_1          <- Add_Merge_ID_17_18_H1N1(df_sub50AF_md_n_SNVs_by_vax_1)
df_mut_bins_prop_noNeut_1              <- Add_Merge_ID_17_18_H1N1(df_mut_bins_prop_noNeut_1)
df_compiled_1                          <- Add_Merge_ID_17_18_H1N1(df_compiled_1)
df_divergence_1                        <- Add_Merge_ID_17_18_H1N1(df_divergence_1)
ds_divergence_1                        <- Add_Merge_ID_17_18_H1N1(ds_divergence_1)
df_prop_shared_1                       <- Add_Merge_ID_17_18_H1N1(df_prop_shared_1)
df_JTplot_1                            <- Add_Merge_ID_17_18_H1N1(df_JTplot_1)
df_household_1                         <- Add_Merge_ID_17_18_H1N1(df_household_1)
sg_1                                   <- Add_Merge_ID_17_18_H1N1(sg_1)
sg_sw_1                                <- Add_Merge_ID_17_18_H1N1(sg_sw_1)
sg_paired_1                            <- Add_Merge_ID_17_18_H1N1(sg_paired_1)
df_mut_bins_prop_1                     <- Add_Merge_ID_17_18_H1N1(df_mut_bins_prop_1)
ds_sg_sw_1                             <- Add_Merge_ID_17_18_H1N1(ds_sg_sw_1)
df_sub50AF_md_n_SNVs_by_bin_and_mut_1  <- Add_Merge_ID_17_18_H1N1(df_sub50AF_md_n_SNVs_by_bin_and_mut_1)
df_sg_AA_HA_1                          <- Add_Merge_ID_17_18_H1N1(df_sg_AA_HA_1)
df_AA_HA_1                             <- Add_Merge_ID_17_18_H1N1(df_AA_HA_1)
df_AA_HA_paired_1                      <- Add_Merge_ID_17_18_H1N1(df_AA_HA_paired_1)
df_mut_bins_prop_all_1                 <- Add_Merge_ID_17_18_H1N1(df_mut_bins_prop_all_1)

#### Import 18-19_H1N1 ####
## 18-19_H3N2
setwd(dir_1819_H1N1); getwd(); dir()

df_all_iSNVs_2                         <- read.csv("df_all_iSNVs.csv")
df_sub50AF_2                           <- read.csv("df_sub50AF.csv")
df_sub50AF_md_n_SNVs_by_subtype_2      <- read.csv("df_sub50AF_md_n_SNVs_by_subtype.csv")
df_sub50AF_md_n_SNVs_by_vax_2          <- read.csv("df_sub50AF_md_n_SNVs_by_vax.csv")
df_mut_bins_prop_noNeut_2              <- read.csv("df_mut_bins_prop_noNeut.csv")
df_compiled_2                          <- read.csv("df_compiled.csv")
df_divergence_2                        <- read.csv("df_divergence.csv")
ds_divergence_2                        <- read.csv("ds_divergence.csv")
df_prop_shared_2                       <- read.csv("df_prop_shared.csv")
df_JTplot_2                            <- read.csv("df_JTplot.csv")
df_household_2                         <- read.csv("df_household.csv")
sg_2                                   <- read.csv("sg.csv")
sg_sw_2                                <- read.csv("sg_sw.csv")
sg_paired_2                            <- read.csv("sg_paired.csv")
df_mut_bins_prop_2                     <- read.csv("df_mut_bins_prop.csv")
ds_sg_sw_2                             <- read.csv("ds_sg_sw.csv")
df_sub50AF_md_n_SNVs_by_bin_and_mut_2  <- read.csv("df_sub50AF_md_n_SNVs_by_bin_and_mut.csv")
df_sg_AA_HA_2                          <- read.csv("df_sg_AA_HA.csv")
df_AA_HA_2                             <- read.csv("df_AA_HA.csv")
df_AA_HA_paired_2                      <- read.csv("df_AA_HA_paired.csv")
df_mut_bins_prop_all_2                 <- read.csv("df_mut_bins_prop_all.csv")

## Merge_ID
Add_Merge_ID_18_19_H1N1 <- function(data) {
  data$Merge_ID <- "18-19_H1N1"
  return(data)
}

## Execute 
df_all_iSNVs_2                         <- Add_Merge_ID_18_19_H1N1(df_all_iSNVs_2)
df_sub50AF_2                           <- Add_Merge_ID_18_19_H1N1(df_sub50AF_2)
df_sub50AF_md_n_SNVs_by_subtype_2      <- Add_Merge_ID_18_19_H1N1(df_sub50AF_md_n_SNVs_by_subtype_2)
df_sub50AF_md_n_SNVs_by_vax_2          <- Add_Merge_ID_18_19_H1N1(df_sub50AF_md_n_SNVs_by_vax_2)
df_mut_bins_prop_noNeut_2              <- Add_Merge_ID_18_19_H1N1(df_mut_bins_prop_noNeut_2)
df_compiled_2                          <- Add_Merge_ID_18_19_H1N1(df_compiled_2)
df_divergence_2                        <- Add_Merge_ID_18_19_H1N1(df_divergence_2)
ds_divergence_2                        <- Add_Merge_ID_18_19_H1N1(ds_divergence_2)
df_prop_shared_2                       <- Add_Merge_ID_18_19_H1N1(df_prop_shared_2)
df_JTplot_2                            <- Add_Merge_ID_18_19_H1N1(df_JTplot_2)
df_household_2                         <- Add_Merge_ID_18_19_H1N1(df_household_2)
sg_2                                   <- Add_Merge_ID_18_19_H1N1(sg_2)
sg_sw_2                                <- Add_Merge_ID_18_19_H1N1(sg_sw_2)
sg_paired_2                            <- Add_Merge_ID_18_19_H1N1(sg_paired_2)
df_mut_bins_prop_2                     <- Add_Merge_ID_18_19_H1N1(df_mut_bins_prop_2)
ds_sg_sw_2                             <- Add_Merge_ID_18_19_H1N1(ds_sg_sw_2)
df_sub50AF_md_n_SNVs_by_bin_and_mut_2  <- Add_Merge_ID_18_19_H1N1(df_sub50AF_md_n_SNVs_by_bin_and_mut_2)
df_sg_AA_HA_2                          <- Add_Merge_ID_18_19_H1N1(df_sg_AA_HA_2)
df_AA_HA_2                             <- Add_Merge_ID_18_19_H1N1(df_AA_HA_2)
df_AA_HA_paired_2                      <- Add_Merge_ID_18_19_H1N1(df_AA_HA_paired_2)
df_mut_bins_prop_all_2                 <- Add_Merge_ID_18_19_H1N1(df_mut_bins_prop_all_2)

#### Merge H1N1 ####
df_all_iSNVs_H1N1                          <- rbind(df_all_iSNVs_1, df_all_iSNVs_2)
df_sub50AF_H1N1                            <- rbind(df_sub50AF_1, df_sub50AF_2)
df_sub50AF_md_n_SNVs_by_subtype_H1N1       <- rbind(df_sub50AF_md_n_SNVs_by_subtype_1, df_sub50AF_md_n_SNVs_by_subtype_2)
df_sub50AF_md_n_SNVs_by_vax_H1N1           <- rbind(df_sub50AF_md_n_SNVs_by_vax_1, df_sub50AF_md_n_SNVs_by_vax_2)
df_mut_bins_prop_noNeut_H1N1               <- rbind(df_mut_bins_prop_noNeut_1, df_mut_bins_prop_noNeut_2)
df_compiled_H1N1                           <- rbind(df_compiled_1, df_compiled_2)
df_divergence_H1N1                         <- rbind(df_divergence_1, df_divergence_2)
ds_divergence_H1N1                         <- rbind(ds_divergence_1, ds_divergence_2)
df_prop_shared_H1N1                        <- rbind(df_prop_shared_1, df_prop_shared_2)
df_JTplot_H1N1                             <- rbind(df_JTplot_1, df_JTplot_2)
df_household_H1N1                          <- rbind(df_household_1, df_household_2)
sg_H1N1                                    <- rbind(sg_1, sg_2)
sg_sw_H1N1                                 <- rbind(sg_sw_1, sg_sw_2)
sg_paired_H1N1                             <- rbind(sg_paired_1, sg_paired_2)
df_mut_bins_prop_H1N1                      <- rbind(df_mut_bins_prop_1, df_mut_bins_prop_2)
ds_sg_sw_H1N1                              <- rbind(ds_sg_sw_1, ds_sg_sw_2)
df_sub50AF_md_n_SNVs_by_bin_and_mut_H1N1   <- rbind(df_sub50AF_md_n_SNVs_by_bin_and_mut_1, df_sub50AF_md_n_SNVs_by_bin_and_mut_2)
df_sg_AA_HA_H1N1                           <- rbind(df_sg_AA_HA_1, df_sg_AA_HA_2)
df_AA_HA_H1N1                              <- rbind(df_AA_HA_1, df_AA_HA_2)
df_AA_HA_paired_H1N1                       <- rbind(df_AA_HA_paired_1, df_AA_HA_paired_2)
df_mut_bins_prop_all_H1N1                  <- rbind(df_mut_bins_prop_all_1, df_mut_bins_prop_all_2)

#### Merge H3N2 and H1N1 merged dfs ####
df_all_iSNVs                        <- rbind(df_all_iSNVs_H3N2, df_all_iSNVs_H1N1)
df_sub50AF                          <- rbind(df_sub50AF_H3N2, df_sub50AF_H1N1)
df_sub50AF_md_n_SNVs_by_subtype     <- rbind(df_sub50AF_md_n_SNVs_by_subtype_H3N2, df_sub50AF_md_n_SNVs_by_subtype_H1N1)
df_sub50AF_md_n_SNVs_by_vax         <- rbind(df_sub50AF_md_n_SNVs_by_vax_H3N2, df_sub50AF_md_n_SNVs_by_vax_H1N1)
df_mut_bins_prop_noNeut             <- rbind(df_mut_bins_prop_noNeut_H3N2, df_mut_bins_prop_noNeut_H1N1)
df_compiled                         <- rbind(df_compiled_H3N2, df_compiled_H1N1)
df_divergence                       <- rbind(df_divergence_H3N2, df_divergence_H1N1)
ds_divergence                       <- rbind(ds_divergence_H3N2, ds_divergence_H1N1)
df_prop_shared                      <- rbind(df_prop_shared_H3N2, df_prop_shared_H1N1)
df_JTplot                           <- rbind(df_JTplot_H3N2, df_JTplot_H1N1)
df_household                        <- rbind(df_household_H3N2, df_household_H1N1)
sg                                  <- rbind(sg_H3N2, sg_H1N1)
sg_sw                               <- rbind(sg_sw_H3N2, sg_sw_H1N1)
sg_paired                           <- rbind(sg_paired_H3N2, sg_paired_H1N1)
df_mut_bins_prop                    <- rbind(df_mut_bins_prop_H3N2, df_mut_bins_prop_H1N1)
ds_sg_sw                            <- rbind(ds_sg_sw_H3N2, ds_sg_sw_H1N1)
df_sub50AF_md_n_SNVs_by_bin_and_mut <- rbind(df_sub50AF_md_n_SNVs_by_bin_and_mut_H3N2, df_sub50AF_md_n_SNVs_by_bin_and_mut_H1N1)
df_sg_AA_HA                         <- rbind(df_sg_AA_HA_H3N2, df_sg_AA_HA_H1N1)
df_AA_HA                            <- rbind(df_AA_HA_H3N2, df_AA_HA_H1N1)
df_AA_HA_paired                     <- rbind(df_AA_HA_paired_H3N2, df_AA_HA_paired_H1N1)
df_mut_bins_prop_all                <- rbind(df_mut_bins_prop_all_H3N2, df_mut_bins_prop_all_H1N1)

#### SAVE ####
setwd(dir_s); getwd(); dir()

write.csv(df_all_iSNVs,                         "df_all_iSNVs.csv")
write.csv(df_sub50AF,                           "df_sub50AF.csv")
write.csv(df_sub50AF_md_n_SNVs_by_subtype,      "df_sub50AF_md_n_SNVs_by_subtype.csv")
write.csv(df_sub50AF_md_n_SNVs_by_vax,          "df_sub50AF_md_n_SNVs_by_vax.csv")
write.csv(df_mut_bins_prop_noNeut,              "df_mut_bins_prop_noNeut.csv")
write.csv(df_compiled,                          "df_compiled.csv")
write.csv(df_divergence,                        "df_divergence.csv")
write.csv(ds_divergence,                        "ds_divergence.csv")
write.csv(df_prop_shared,                       "df_prop_shared.csv")
write.csv(df_JTplot,                            "df_JTplot.csv")
write.csv(df_household,                         "df_household.csv")
write.csv(sg,                                   "sg.csv")
write.csv(sg_sw,                                "sg_sw.csv")
write.csv(sg_paired,                            "sg_paired.csv")
write.csv(df_mut_bins_prop,                     "df_mut_bins_prop.csv")
write.csv(ds_sg_sw,                             "ds_sg_sw.csv")
write.csv(df_sub50AF_md_n_SNVs_by_bin_and_mut,  "df_sub50AF_md_n_SNVs_by_bin_and_mut.csv")
write.csv(df_sg_AA_HA,                          "df_sg_AA_HA.csv")
write.csv(df_AA_HA,                             "df_AA_HA.csv")
write.csv(df_AA_HA_paired,                      "df_AA_HA_paired.csv")
write.csv(df_mut_bins_prop_all,                 "df_mut_bins_prop_all.csv")


