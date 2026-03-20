## Clear Global Environment
rm(list = ls())

## Directories
dir    <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/dfs",sep="")
dir_md <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/summary",sep="")
dir_s  <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/figures",sep="")

#### Session prep ####
## Install packages and load libraries as required
if(!require(ggbeeswarm)){
  install.packages("ggbeeswarm", dependencies = T)
  library(ggbeeswarm)
}
if(!require(scales)){
  install.packages("scales", dependencies = T)
  library(scales)
}
if(!require(ggpubr)){
  install.packages("ggpubr", dependencies = T)
  library(ggpubr)
}
if(!require(tidyverse)){
  install.packages("tidyverse", dependencies = T)
  library(tidyverse)
}
if(!require(lubridate)){
  install.packages("lubridate", dependencies = T)
  library(lubridate)
}
if(!require(ggplot2)){
  install.packages("ggplot2", dependencies = T)
  library(ggplot2)
}
if(!require(grid)){
  install.packages("grid", dependencies = T)
  library(grid)
}
if(!require(gridExtra)){
  install.packages("gridExtra", dependencies = T)
  library(gridExtra)
}
if(!require(cowplot)){
  install.packages("cowplot", dependencies = T)
  library(cowplot)
}
if(!require(reshape2)){
  install.packages("reshape2", dependencies = T)
  library(reshape2)
}
if(!require(vcfR)){
  install.packages("vcfR", dependencies = T)
  library(vcfR)
}

#### Plot standards ####
## boxplots and quasi-random dots
quasi_size    = 1.75
quasi_alpha   = 0.60
boxplot_alpha = 0.30

## values
y_nuc_div <- 0.00065
y_pinpis1  <- 0.001
y_pinpis2  <- 0.006
y_piX1 <- .0005
y_piX2 <- .006

## factors
factor_by_gene <- c("PB2", "PB1", "PA", "HA", "NP", "NA")
factor_by_segment <- c("PB2", "PB1", "PA", "HA", "NP", "Neuraminidase", "M1", "M2", "NS1", "NEP")
factor_by_gene_cleaned <- c("PB2", "PB1", "PA", "HA", "Anti. HA", "Nonanti. HA", "NP", "NA", "M1", "M2", "NS1", "NEP", "PA-X", "PB1-F2"   )

## theme
size = 7
axis_formatting <- theme(axis.text.x = element_text(size = size),
                         axis.text.y = element_text(size = size),
                         axis.title.x = element_text(size = size, margin = margin(t = 4)),
                         axis.title.y = element_text(size = size, margin = margin(r = 4)),
                         title = element_text(size = size))

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
# donor and recipient
palette_dr <- c("Donor" = "Black",
                "Recipient" = "Red")

# donor and recipient and within-host
palette_drw <- c("Donor" = "Black",
                 "Recipient" = "Red",
                 "Within-host" = "steelblue2")

# mutation types
palette_muts_NS_S_sw <- c("NS" = "#FF7F20",
                          "S" = "#4F7899")
palette_muts <- c("Nonsynonymous" = "#FF7F20",
                  "Synonymous" = "#4F7899",
                  "Nonsense" = "#7AD9C2")
palette_muts_NS_S_Stopgained <- c("Nonsynonymous" = "#FF7F20",
                                  "Synonymous" = "#4F7899",
                                  "Stop gained" = "#7AD9C2")

# mutation types without stop
palette_muts_NS_S <- c("Nonsynonymous" = "#FF7F20",
                       "Synonymous" = "#4F7899")
palette_muts_NS_S_2 <- c("Nonsynonymous, within-host" = "#FF7F20",
                         "Synonymous, within-host" = "#4F7899", 
                         "Nonsynonymous, between-host" = "#FDA564",
                         "Synonymous, between-host" = "#6C9EC6")

# antigenicity
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


# vax and unvax
palette_vax <- c("Vaccinated"   = "#8073ac",
                 "Unvaccinated" = "#d6604d",
                 "Either"       = "grey")
# pairs
palette_pairs <- c("Household" = "Red", 
                   "Random" = "Grey")

#### Functions ####
# Define the analyze_data function
analyze_data <- function(Name_of_data_frame, Group_1, Group_2, test) {
  # Extract the data
  data <- get(Name_of_data_frame)
  
  # Convert ID2 to character type
  data$ID2 <- as.character(data$ID2)
  
  # Subset data for each group
  group1 <- data[data$ID2 == Group_1, "V1"]
  group2 <- data[data$ID2 == Group_2, "V1"]
  
  # Check if any group has all zero values
  if (all(group1 == 0) || all(group2 == 0)) {
    skip_reason <- ifelse(all(group1 == 0), paste(Group_1, "has all zero values"), paste(Group_2, "has all zero values"))
    message("Skipping statistical analyses because ", skip_reason)
    
    # Fill out the result data frame with information on why it was skipped
    result_df <- data.frame(
      Name_of_data_frame = Name_of_data_frame,
      Group_1 = Group_1,
      Group_2 = Group_2,
      Skip_Reason = skip_reason,
      test_used = NA,
      p_value = NA,
      is_significant = NA
    )
    
    return(result_df)
  }
  
  # Perform the selected test
  if (test == "paired_t-test") {
    p_value <- t.test(group1, group2, paired = TRUE)$p.value
    test_used <- "Paired t-test"
  } else if (test == "one_way_anova") {
    p_value <- anova(lm(V1 ~ ID2, data = data))$`Pr(>F)`[1]
    test_used <- "One-way ANOVA"
  } else {
    stop("Invalid test specified. Choose 'paired_t-test' or 'one_way_anova'.")
  }
  
  # Apply multiple hypothesis correction
  num_tests <- length(group_combinations_H3)  # Adjust based on your actual situation
  p_value_corrected <- p.adjust(p_value, method = "bonferroni", n = num_tests)
  
  # Determine if the corrected p-value is significant
  is_significant <- p_value_corrected < 0.05
  
  # Create data frame
  result_df <- data.frame(
    Name_of_data_frame = Name_of_data_frame,
    Group_1 = Group_1,
    Group_2 = Group_2,
    Skip_Reason = NA,
    test_used = test_used,
    p_value = p_value_corrected,
    is_significant = is_significant
  )
  
  return(result_df)
}

## data summary
ds <- function(data, varname, groupnames) {
  require(plyr)
  summary_func <- function(x, col) {
    c(mean = mean(x[[col]], na.rm=TRUE),
      median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm = TRUE) / sqrt(sum(!is.na(x[[col]]))))
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

## bayesian bootstrap
b <- function(data){
  return(bayesboot(data, mean, R = 10))
  #return(bayesboot(data, mean, R = 10000))
}

#### Import data ####
## import metadata
setwd(dir_md); getwd(); dir()
md <- as.data.frame(read_csv("df_samplelist_metadata.csv"))

qc_pass <- c("REDACTED")

# pull all that pass QC
temp_md <- md[md$WSLH_ID %in% qc_pass,]; length(temp_md$WSLH_ID); round(length(temp_md$WSLH_ID)/384,2)
# sort by season and subtype
temp_md_1718H1N1 <- temp_md[temp_md$Season=="17-18" & temp_md$Subtype=="H1N1",]; length(temp_md_1718H1N1$WSLH_ID); round(length(temp_md_1718H1N1$WSLH_ID)/27,2)
temp_md_1718H3N2 <- temp_md[temp_md$Season=="17-18" & temp_md$Subtype=="H3N2",]; length(temp_md_1718H3N2$WSLH_ID); round(length(temp_md_1718H3N2$WSLH_ID)/163,2)
temp_md_1819H1N1 <- temp_md[temp_md$Season=="18-19" & temp_md$Subtype=="H1N1",]; length(temp_md_1819H1N1$WSLH_ID); round(length(temp_md_1819H1N1$WSLH_ID)/79,2)
temp_md_1819H3N2 <- temp_md[temp_md$Season=="18-19" & temp_md$Subtype=="H3N2",]; length(temp_md_1819H3N2$WSLH_ID); round(length(temp_md_1819H3N2$WSLH_ID)/115,2)
# calculate n_female
table(temp_md_1718H1N1$gender)[1]; round(unname(table(temp_md_1718H1N1$gender)[1]) / sum(unname(table(temp_md_1718H1N1$gender)[1]) + unname(table(temp_md_1718H1N1$gender)[2])),2)
table(temp_md_1718H3N2$gender)[1]; round(unname(table(temp_md_1718H3N2$gender)[1]) / sum(unname(table(temp_md_1718H3N2$gender)[1]) + unname(table(temp_md_1718H3N2$gender)[2])),2)
table(temp_md_1819H1N1$gender)[1]; round(unname(table(temp_md_1819H1N1$gender)[1]) / sum(unname(table(temp_md_1819H1N1$gender)[1]) + unname(table(temp_md_1819H1N1$gender)[2])),2)
table(temp_md_1819H3N2$gender)[1]; round(unname(table(temp_md_1819H3N2$gender)[1]) / sum(unname(table(temp_md_1819H3N2$gender)[1]) + unname(table(temp_md_1819H3N2$gender)[2])),2)
table(temp_md$gender)[1]; round(unname(table(temp_md$gender)[1]) / sum(unname(table(temp_md$gender)[1]) + unname(table(temp_md$gender)[2])),2)
# calculate n_under18
length(temp_md_1718H1N1$age[temp_md_1718H1N1$age<18]); round(length(temp_md_1718H1N1$age[temp_md_1718H1N1$age<18])/24,2)
length(temp_md_1718H3N2$age[temp_md_1718H3N2$age<18]); round(length(temp_md_1718H3N2$age[temp_md_1718H3N2$age<18])/119,2)
length(temp_md_1819H1N1$age[temp_md_1819H1N1$age<18]); round(length(temp_md_1819H1N1$age[temp_md_1819H1N1$age<18])/70,2)
length(temp_md_1819H3N2$age[temp_md_1819H3N2$age<18]); round(length(temp_md_1819H3N2$age[temp_md_1819H3N2$age<18])/70,2)
length(temp_md$age[temp_md$age<18]); round(length(temp_md$age[temp_md$age<18])/283,2)
# calculate n_vax
table(temp_md_1718H1N1$flu_vaccine)[2]; round(unname(table(temp_md_1718H1N1$flu_vaccine)[2]) / 24,2)
table(temp_md_1718H3N2$flu_vaccine)[2]; round(unname(table(temp_md_1718H3N2$flu_vaccine)[2]) / 119,2)
table(temp_md_1819H1N1$flu_vaccine)[2]; round(unname(table(temp_md_1819H1N1$flu_vaccine)[2]) / 70,2)
table(temp_md_1819H3N2$flu_vaccine)[2]; round(unname(table(temp_md_1819H3N2$flu_vaccine)[2]) / 70,2)
table(temp_md$flu_vaccine)[2]; round(unname(table(temp_md$flu_vaccine)[2]) / 283,2)
# calculate Ct mean and sd of each year
round(mean(temp_md_1718H1N1$Ct), 1); round(sd(temp_md_1718H1N1$Ct), 1)
round(mean(temp_md_1718H3N2$Ct), 1); round(sd(temp_md_1718H3N2$Ct), 1)
round(mean(temp_md_1819H1N1$Ct), 1); round(sd(temp_md_1819H1N1$Ct), 1)
round(mean(temp_md_1819H3N2$Ct), 1); round(sd(temp_md_1819H3N2$Ct), 1)
round(mean(temp_md$Ct), 1); round(sd(temp_md$Ct), 1)
min(temp_md$Ct); max(temp_md$Ct)

## import generated data
setwd(dir); getwd(); dir()
df_all_iSNVs <- read.csv("df_all_iSNVs.csv")                                               # all iSNVs from both 17-18 H3N2 and 18-19 H3N2
df_sub50AF <- read.csv("df_sub50AF.csv")                                                   # same as previous but only 50% frequency iSNVs and below
df_sub50AF_md_n_SNVs_by_subtype <- read.csv("df_sub50AF_md_n_SNVs_by_subtype.csv")         # number of iSNVs by subtype
df_sub50AF_md_n_SNVs_by_vax <- read.csv("df_sub50AF_md_n_SNVs_by_vax.csv")                 # number of iSNVs by subtype and vax status
df_sub50AF_md_n_SNVs_by_bin_and_mut <- read.csv("df_sub50AF_md_n_SNVs_by_bin_and_mut.csv") # number of iSNVs by mutation_type bin-proportion
df_mut_bins_prop_noNeut <- read.csv("df_mut_bins_prop_noNeut.csv")                         # count and mut_type-bin-proportion of iSNVs by subtype and mutation type w/o neutral expectation
df_mut_bins_prop <- read.csv("df_mut_bins_prop.csv")                                       # same as above, but with neutral expectation
df_mut_bins_prop_all <- read.csv("df_mut_bins_prop_all.csv")                               # same as df_sub50AF_md_n_SNVs_by_bin_and_mut but with proportions calculated per-sample
df_compiled <- read.csv("df_compiled.csv")                                                 # count and proportion of mutations shared by random pairings
df_household <- read.csv("df_household.csv")                                               # count and proportion of mutations shared by household pairings
df_divergence <- read.csv("df_divergence.csv")                                             # divergence and pi by sample-gene
ds_divergence <- read.csv("ds_divergence.csv")                                             # stats of df_divergence by gene
df_JTplot <- read.csv("df_JTplot.csv")                                                     # frequency of iSNVs shared between donor and recipient
sg <- read.csv("sg.csv")                                                                   # snpgenie output of gene-wise nucleotide diversity values
sg_sw <- read.csv("sg_sw.csv")                                                             # snpgenie output of sliding-window nucleotide diversity values
sg_paired <- read.csv("sg_paired.csv")                                                     # reconstructed sg to frame donor and recipient relationship
ds_sg_sw <- read.csv("ds_sg_sw.csv")                                                       # stats of sliding window pi stats for both 17-18 H3N2 and 18-19 H3N2
df_sg_AA_HA <- read.csv("df_sg_AA_HA.csv")                                                 # snpgenie output of sliding-window nucleotide diversity values, per AA
df_AA_HA <- read.csv("df_AA_HA.csv")                                                       # per-site piN-piS for each antigenic site and antigenic/non-antigenic for each sample
df_AA_HA_paired <- read.csv("df_AA_HA_paired.csv")                                         # reconstructed df_AA_HA to frame donor and recipient relationship
df_prop_shared <- read.csv("df_prop_shared.csv")                                           # count and proportion of mutations shared by random pairings for both 17-18 H3N2 and 18-19 H3N2
sg_sw_antigenic_H3 <- read_tsv("antigenic_regions_H3.tsv")                                 # antigenic regions as defined by the literature for H3
sg_sw_antigenic_H1 <- read_tsv("antigenic_regions_H1.tsv")                                 # antigenic regions as defined by the literature for H1

#### Figure 3ab - Plot: pairing, JT plot ####
## labels 
Fig5alabels <- c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
## create bin range
bin_range <- data.frame("Lower" = seq(0,0.95,0.05), "Upper" = seq(0.05,1,0.05))
bin_range <- rbind(bin_range, data.frame("Lower" = 1, "Upper" = 1))
bin_range$Bins <- paste(bin_range[,1], bin_range[,2], sep="-")

## recalculate df_prop_shared prop_var_shared for all pairs across seasons/subtypes
df_prop_shared_all_r <- df_prop_shared[df_prop_shared$Var2=="Random",]
df_prop_shared_all_h <- df_household
df_compiled_h <- df_prop_shared_all_h
df_compiled_r <- df_prop_shared_all_r
df_compiled_h$Bins <- ""
df_compiled_r$Bins <- ""

## add Prop_Var_Shared bins
# this isn't the fastest method, but it's adjustable (bin_range above)
# these nested loops parse every iSNV and bin it based on its AF
n <- length(df_compiled_h[,1])
m <- length(bin_range[,1])
start = Sys.time()
for (i in 1:n) {
  for (j in 1:m) {
    df_compiled_h[i,]$Bins[df_compiled_h[i,]$Prop_Var_Shared >= bin_range[j,1] && 
                             df_compiled_h[i,]$Prop_Var_Shared < bin_range[j,2]] <- bin_range[j,3]
    df_compiled_h[i,]$Bins[df_compiled_h[i,]$Prop_Var_Shared == 0.3] <- "0.3-0.35"
    df_compiled_h[i,]$Bins[df_compiled_h[i,]$Prop_Var_Shared == 1.0] <- "0.95-1"
    # not sure why 0.3 isn't mapping
  }
  out <- paste(i, "/", n, ": binning iSNVs..."); print(out)
}
finish <- Sys.time()
out1 <- paste("Start: ", start)
out2 <- paste("Finish:", finish)
print(out1); print(out2)
df_compiled_h$Pairing <- "Household"

## split df_compiled_h into bin and bin frequency
df_compiled_h <- as.data.frame(table(df_compiled_h$Bins, df_compiled_h$Pairing))
list_compiled <- split(df_compiled_h, with(df_compiled_h, Var2, drop=T))
n <- length(list_compiled)
m <- length(list_compiled[[1]][,1])
for (i in 1:n) {
  for (j in 1:m) {
    list_compiled[[i]]$Prop[j] <- sum(list_compiled[[i]][j,3]/sum(list_compiled[[i]][,3]))
  }
}
# back to df
df_compiled_h <- Reduce(full_join, list_compiled)
colnames(df_compiled_h)[1] <- "Bins"

## aggregate df_compiled_r into total bin and bin frequency
df_compiled_r <- aggregate(Freq~Var1, df_compiled_r, FUN=sum)
df_compiled_r$Var2 <- "Random"
colnames(df_compiled_r)[1] <- "Bins"

## split df_compiled_r into bin and bin frequency
list_compiled <- split(df_compiled_r, with(df_compiled_r, Var2, drop=T))
n <- length(list_compiled)
m <- length(list_compiled[[1]][,1])
for (i in 1:n) {
  for (j in 1:m) {
    list_compiled[[i]]$Prop[j] <- sum(list_compiled[[i]][j,2]/sum(list_compiled[[i]][,2]))
  }
}
# back to df
df_compiled_r <- Reduce(full_join, list_compiled)

## put df_prop_shared_all_r and df_prop_shared_all_h back together from df_compiled_h
df_prop_shared_all <- rbind(df_compiled_r, df_compiled_h)
df_prop_shared_all$Merge_ID <- "All"
df_prop_shared_all$Freq <- 0
df_compiled$X.1 <- NULL
df_compiled$X <- NULL
colnames(df_compiled)[1] <- "Bins"

df_compiled <- rbind(df_compiled, df_prop_shared_all[df_prop_shared_all$Var2=="Random",])

## Stats!
# calculate threshold for each merge_ID
n <- length(levels(as.factor(df_compiled$Merge_ID)))
p05_threshold <- list()
p05_pairs <- list()
for (i in 1:n) {
  j <- levels(as.factor(df_compiled$Merge_ID))[i]
  # calculate threshold for each merge_ID
  df_compiled_pvalue <- quantile(df_compiled$Prop[df_compiled$Merge_ID==j], 
                                 probs = c(0.50, 0.75, 0.9, 0.95, 0.99))
  p05_threshold[[i]] <- as.numeric(gsub("95%","",df_compiled_pvalue[4]))
  # calculate pairs that pass their own season/subtype threshold
  household_pairs_sigdif <- df_household$ID[df_household$Merge_ID==j &
                                              df_household$Prop_Var_Shared >= p05_threshold[[i]]]
  p05_pairs[[i]] <- household_pairs_sigdif
  # give names to each list entry
  names(p05_threshold)[[i]] <- j
  names(p05_pairs)[[i]] <- j
}

# pairs to keep for donor-recipient analyses
str_p05_pairs <- c(p05_pairs[[1]], p05_pairs[[2]], p05_pairs[[3]], p05_pairs[[4]])
# WSLH_IDs that pass genetic criterion
str_p05_pair_WSLHID <- separate(as.data.frame(str_p05_pairs), str_p05_pairs, c("Donor", "Recipient"), sep = "-")
str_p05_pair_WSLHID <- unique(c(str_p05_pair_WSLHID[,1], str_p05_pair_WSLHID[,2]))

p <- which(round(p05_threshold$All,2) < bin_range$Upper & round(p05_threshold$All,2) >= bin_range$Lower)
p <- p+.5



## Validated donor-recipient pairs over time between collection
temp <- md[md$WSLH_ID %in% str_p05_pair_WSLHID,]
df <- temp %>%
  dplyr::mutate(DoC = mdy(DoC))

# Calculate days from donor sample
df_diff <- df %>%
  dplyr::group_by(Num) %>%
  dplyr::mutate(Index_DoC = first(DoC[Case == "Index"])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Day_Diff = as.numeric(DoC - Index_DoC))

# Assign a unique Pair_ID per donor-recipient pair
df_diff <- df_diff %>%
  dplyr::group_by(Num) %>%
  dplyr::arrange(DoC) %>%
  dplyr::mutate(
    pair_number = cumsum(Case == "Contact"),
    Pair_ID = paste0(Num, "_", if_else(Case == "Index", lead(pair_number), pair_number))
  ) %>%
  dplyr::ungroup()

# Order by contact sampling delay
ordered_levels <- df_diff %>%
  dplyr::filter(Case == "Contact") %>%
  dplyr::group_by(Num) %>%
  dplyr::summarise(max_delay = max(Day_Diff, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(max_delay) %>%
  dplyr::pull(Num)

# Apply factor levels to Pair_ID and Num
df_diff <- df_diff %>%
  dplyr::mutate(Pair_ID = factor(Pair_ID, levels = ordered_levels))

df_diff <- df_diff %>%
  dplyr::mutate(Num = factor(Num, levels = ordered_levels))

# Adjust recipients with overlapping collection dates
df_diff <- df_diff %>%
  dplyr::group_by(Num) %>%
  dplyr::arrange(DoC) %>%
  dplyr::mutate(
    donor_days = list(Day_Diff[Case == "Index"]),
    # Offset if recipient shares day with a donor
    Day_Diff = if_else(
      Case == "Contact" & Day_Diff %in% unlist(donor_days),
      Day_Diff + 0.125,
      Day_Diff
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Num, Day_Diff) %>%
  dplyr::mutate(
    dup_count = sum(Case == "Contact"),
    contact_index = if_else(Case == "Contact", row_number(), NA_integer_),
    Day_Diff = if_else(
      Case == "Contact" & dup_count > 1,
      Day_Diff + 0.125 * (contact_index - 1),
      Day_Diff
    )
  ) %>%
  dplyr::ungroup()



# Plot
plot_dr_sampling <-  ggplot(df_diff, aes(x = Day_Diff, y = Num, color = Case)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray50") +
  geom_line(aes(group = Num), color = "gray60", linetype = "dotted") +
  geom_point(size = 2) +
  labs(
    x = "Days from donor sample collection date",
    y = NULL,
    title = NULL
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0,8,1), labels = seq(0,8,1)) + 
  scale_color_manual(values = c("Index" = "black", "Contact" = "red")) +  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )


## totals before the filter
length(unique(df_JTplot$Num)) # number of households (and donors)      -- should be 36
length(unique(df_JTplot$Donor_WSLH_ID)) # number of donors (to verify) -- should be 36
length(unique(df_JTplot$ID))  # number of donor-recipient pairs        -- should be 46

length(unique(sg_paired$Num)) # number of households (and donors)      -- should be 36
length(unique(sg_paired$Donor_WSLH_ID)) # number of donors (to verify) -- should be 36
length(unique(sg_paired$ID))  # number of donor-recipient pairs        -- should be 46

length(unique(df_AA_HA_paired$Num)) # number of households (and donors)      -- should be 36
length(unique(df_AA_HA_paired$Donor_WSLH_ID)) # number of donors (to verify) -- should be 36
length(unique(df_AA_HA_paired$ID))  # number of donor-recipient pairs        -- should be 46

# calculate donors, recipients, and donor-recipient pairs per Merge_ID
summary_df <- df_AA_HA_paired %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::summarise(
    n_donors   = n_distinct(Donor_WSLH_ID),
    n_recipients = n_distinct(Recipient_WSLH_ID),
    n_pairs    = n_distinct(ID)
  ); summary_df

# Count unique recipients per donor
donor_recipient_counts <- df_AA_HA_paired %>%
  dplyr::distinct(Donor_WSLH_ID, Recipient_WSLH_ID) %>%
  dplyr::group_by(Donor_WSLH_ID) %>%
  dplyr::summarise(n_recipients = dplyr::n_distinct(Recipient_WSLH_ID))

# Count how many donors fall into each "number of recipients" category
donor_transmission_summary <- donor_recipient_counts %>%
  dplyr::group_by(n_recipients) %>%
  dplyr::summarise(n_donors = dplyr::n())

# Donor-recipient count per season-subtype
donor_trans_by_season <- df_AA_HA_paired %>%
  dplyr::distinct(Merge_ID, Donor_WSLH_ID, Recipient_WSLH_ID) %>%
  dplyr::group_by(Merge_ID, Donor_WSLH_ID) %>%
  dplyr::summarise(n_recipients = dplyr::n_distinct(Recipient_WSLH_ID), .groups = "drop") %>%
  dplyr::group_by(Merge_ID, n_recipients) %>%
  dplyr::summarise(n_donors = dplyr::n(), .groups = "drop")



## remove invalid donor-recipient pairs from the following data frames
df_JTplot       <- df_JTplot[df_JTplot$ID %in% str_p05_pairs,]
sg_paired       <- sg_paired[sg_paired$ID %in% str_p05_pairs,]
df_AA_HA_paired <- df_AA_HA_paired[df_AA_HA_paired$ID %in% str_p05_pairs,]

## totals after the filter
length(unique(df_JTplot$Num)) # number of households (and donors)      -- should be 31
length(unique(df_JTplot$Donor_WSLH_ID)) # number of donors (to verify) -- should be 31
length(unique(df_JTplot$ID))  # number of donor-recipient pairs        -- should be 40

length(unique(sg_paired$Num)) # number of households (and donors)      -- should be 31
length(unique(sg_paired$Donor_WSLH_ID)) # number of donors (to verify) -- should be 31
length(unique(sg_paired$ID))  # number of donor-recipient pairs        -- should be 40

length(unique(df_AA_HA_paired$Num)) # number of households (and donors)      -- should be 31
length(unique(df_AA_HA_paired$Donor_WSLH_ID)) # number of donors (to verify) -- should be 31
length(unique(df_AA_HA_paired$ID))  # number of donor-recipient pairs        -- should be 40

## totals per season
# calculate donors, recipients, and donor-recipient pairs per Merge_ID
summary_df <- df_AA_HA_paired %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::summarise(
    n_donors   = n_distinct(Donor_WSLH_ID),
    n_recipients = n_distinct(Recipient_WSLH_ID),
    n_pairs    = n_distinct(ID)
  ); summary_df

# Count unique recipients per donor
donor_recipient_counts <- df_AA_HA_paired %>%
  dplyr::distinct(Donor_WSLH_ID, Recipient_WSLH_ID) %>%
  dplyr::group_by(Donor_WSLH_ID) %>%
  dplyr::summarise(n_recipients = dplyr::n_distinct(Recipient_WSLH_ID))

# Count how many donors fall into each "number of recipients" category
donor_transmission_summary <- donor_recipient_counts %>%
  dplyr::group_by(n_recipients) %>%
  dplyr::summarise(n_donors = dplyr::n())

# Donor-recipient count per season-subtype
donor_trans_by_season <- df_AA_HA_paired %>%
  dplyr::distinct(Merge_ID, Donor_WSLH_ID, Recipient_WSLH_ID) %>%
  dplyr::group_by(Merge_ID, Donor_WSLH_ID) %>%
  dplyr::summarise(n_recipients = dplyr::n_distinct(Recipient_WSLH_ID), .groups = "drop") %>%
  dplyr::group_by(Merge_ID, n_recipients) %>%
  dplyr::summarise(n_donors = dplyr::n(), .groups = "drop")




## if no value for an x-axis bin, give it the value 0
df_prop_shared <- df_prop_shared_all
df_zero_r <- data.frame("Bins" = bin_range$Bins, "Var2" = "Random", "Freq" = 0, "Prop" = 0)
df_zero_h <- data.frame("Bins" = bin_range$Bins, "Var2" = "Household", "Freq" = 0, "Prop" = 0)
df_zero <- rbind(df_zero_r, df_zero_h)
df_zero$ID <- paste(df_zero$Bins, df_zero$Var2, sep="-")
df_prop_shared$ID <- paste(df_prop_shared$Bins, df_prop_shared$Var2, sep="-")
df_prop_shared <- merge(df_prop_shared, df_zero, by = "ID", all = T)
df_prop_shared$Bins.x <- df_prop_shared$Bins.y
df_prop_shared$Var2.x <- df_prop_shared$Var2.y
df_prop_shared$Prop.x[is.na(df_prop_shared$Bins.x)] <- 0
df_prop_shared$Freq.x[is.na(df_prop_shared$Bins.x)] <- 0
df_prop_shared$Prop.x[is.na(df_prop_shared$Prop.x)] <- 0
df_prop_shared$Freq.x[is.na(df_prop_shared$Freq.x)] <- 0
df_prop_shared$X <- NULL
df_prop_shared$Prop.y <- NULL
df_prop_shared$Freq.y <- NULL
df_prop_shared$Var2.y <- NULL
df_prop_shared$Bins.y <- NULL
#df_prop_shared <- merge(df_prop_shared, df_zero, by = "ID", all = T)
colnames(df_prop_shared) <- c("ID", "Bins", "Freq", "Var2", "Prop", "Merge_ID")

## Prop pairs by prop var shared, factor by random and pairs
plot_pairing <- ggplot(data = df_prop_shared, 
                aes(x = Bins, y = Prop, fill = Var2)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_fill_manual(values = palette_pairs) + 
  scale_x_discrete(breaks = function(x){x[c(T,F)]}, labels = Fig5alabels, expand = c(0.01, 0.01)) + 
  scale_y_continuous(breaks = seq(0, 1, .10), labels = seq(0, 1, .1), expand = c(0.01, 0.01)) + 
  labs(x = "Proportion of variants shared", y = "Proportion of pairs") + 
  geom_segment(aes(x=p, xend=p,
                   y=0.0, yend=1.0), 
               linetype=2) + 
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.90, 0.85)) + 
  axis_formatting + 
  legend_formatting + 
  background_formatting


## JT plot in groups
temp <- data.frame("iSNV" = paste(df_all_iSNVs$REF, df_all_iSNVs$POS, df_all_iSNVs$ALT, 
                                  "-", df_all_iSNVs$Feature_ID, sep = ""),
                   "Ann" = df_all_iSNVs$Ann)
df_JTplot_ann <- merge(df_JTplot, temp, by = "iSNV")
df_JTplot_ann$Ann <- as.factor(df_JTplot_ann$Ann); levels(df_JTplot_ann$Ann)
df_JTplot_ann$Double_AF <- df_JTplot_ann$Donor_AF + df_JTplot_ann$Recipient_AF
df_JTplot_ann <- df_JTplot_ann[df_JTplot_ann$Double_AF!=0,]

## JT plot
# remove duplicates
df_JTplot_ann$ID2 <- paste(df_JTplot_ann$ID, df_JTplot_ann$Merge_ID, df_JTplot_ann$iSNV, df_JTplot_ann$Num, df_JTplot_ann$Donor_AF, df_JTplot_ann$Recipient_AF)
df_JTplot_ann <- df_JTplot_ann[!duplicated(df_JTplot_ann$ID2),]
# add grouping
df_JTplot_ann$Group <- ""
df_JTplot_ann$Group[df_JTplot_ann$Donor_AF==0 & df_JTplot_ann$Recipient_AF>0] <- "In recipient"
df_JTplot_ann$Group[df_JTplot_ann$Donor_AF>0  & df_JTplot_ann$Recipient_AF==0] <- "In donor"
df_JTplot_ann$Group[df_JTplot_ann$Donor_AF>0  & df_JTplot_ann$Recipient_AF>0] <- "Shared"

Fig5b <- ggplot(data = df_JTplot_ann, color = "black",
                aes(x = Donor_AF, y = Recipient_AF)) + 
  geom_point(alpha = 0.5, stroke = NA) + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) + 
  scale_x_continuous(
    minor_breaks = seq(0, 1, 0.1),
    breaks = seq(0, 1, 0.1),
    labels = ifelse(seq(0, 10, 1) %% 2 == 0, seq(0, 1, 0.1), ""),
    expand = c(0.01, 0.01)
  ) + 
  scale_y_continuous(
    minor_breaks = seq(0, 1, 0.1),
    breaks = seq(0, 1, 0.1),
    labels = ifelse(seq(0, 10, 1) %% 2 == 0, seq(0, 1, 0.1), ""),
    expand = c(0.01, 0.01)
  ) + 
  labs(x="iSNV frequency (Donor)", y="iSNV frequency (Recipient)") + 
  theme(legend.title = element_blank(), legend.key = element_blank(), legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) + 
  axis_formatting + legend_formatting + background_formatting

#### Figure 4   - Plot: (antigenicity) piN_minus_piS between hosts ####
### plotting πN-πS per site for donors, recipients, and all samples (antigenic vs. non-antigenic)
# function for no-zero analysis
plot_antigenic_piN_minus_piS_nozeros <- function(season_subtype, 
                                                 df_wh = df_AA_HA, 
                                                 df_bw = df_AA_HA_paired,
                                                 y_limits = c(NA, NA), 
                                                 title) {
  ## filter
  # between host
  df_bw <- df_bw %>%
    # 1) keep only the season and only the “all” rows,
    #    and at least one of Donor/Recipient is non‐zero
    filter(
      Merge_ID == season_subtype,
      Gene      %in% c("Antigenic_all", "Non-antigenic_all"),
      Donor_perSite_piN_minus_piS  != 0 |
        Recipient_perSite_piN_minus_piS != 0
    ) %>%
    # 2) recode Region & pull out just the two columns we care about
    transmute(
      Region    = recode(
        Gene,
        "Antigenic_all"     = "Antigenic",
        "Non-antigenic_all" = "Non-antigenic"
      ),
      Donor     = Donor_perSite_piN_minus_piS,
      Recipient = Recipient_perSite_piN_minus_piS
    ) %>%
    # 3) pivot to long form
    pivot_longer(
      cols      = c(Donor, Recipient),
      names_to  = "Group",
      values_to = "Value"
    )
  # within host
  df_wh <- df_wh[df_wh$perSite_piN_minus_piS!=0,]
  
  # 1) build the “All” (within-host) data
  df_comm <- df_wh %>%
    filter(Merge_ID == season_subtype) %>%
    transmute(
      Region = if_else(Antigenicity == "Antigenic", "Antigenic", "Non-antigenic"),
      Group  = "All",
      Value  = perSite_piN_minus_piS
    )
  
  ## 2) build the donor/recipient data
  #df_bw_long <- df_bw %>%
  #  filter(
  #    Merge_ID == season_subtype,
  #    Gene %in% c("Antigenic_all", "Non-antigenic_all")
  #  ) %>%
  #  transmute(
  #    Region = recode(
  #      Gene,
  #      "Antigenic_all"     = "Antigenic",
  #      "Non-antigenic_all" = "Non-antigenic"
  #    ),
  #    Donor     = Donor_perSite_piN_minus_piS,
  #    Recipient = Recipient_perSite_piN_minus_piS
  #  ) %>%
  #  pivot_longer(
  #    cols      = c(Donor, Recipient),
  #    names_to  = "Group",
  #    values_to = "Value"
  #  )
  
  # 3) combine & factor
  df_plot <- bind_rows(df_comm, df_bw) %>%
    mutate(
      Region = factor(Region, levels = c("Antigenic", "Non-antigenic")),
      Group  = factor(Group,  levels = c("Donor", "Recipient", "All"))
    )
  
  # 4) normality tests by Region × Group
  safe_shapiro <- function(x) {
    if (length(unique(x)) == 1) {
      return(list(W = NA_real_, p_value = NA_real_, is_normal = FALSE))
    }
    test <- shapiro.test(x)
    list(W = unname(test$statistic), p_value = test$p.value, is_normal = (test$p.value > 0.05))
  }
  
  normality_results <- df_plot %>%
    dplyr::group_by(Region, Group) %>%
    dplyr::reframe(
      n = n(),
      stats = list(safe_shapiro(Value))
    ) %>%
    unnest_wider(stats)
  
  #normality_results <- df_plot %>%
  #  group_by(Region, Group) %>%
  #  dplyr::summarise(
  #    n       = n(),
  #    W       = shapiro.test(Value)$statistic,
  #    p_value = shapiro.test(Value)$p.value,
  #    is_normal = (p_value > 0.05),
  #    .groups = "drop"
  #  )
  print(normality_results)
  
  # 5) paired Wilcoxon tests (Donor vs Recipient) *within each Region*
  wilcox_results <- df_bw %>%
    group_by(Region) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(
        Value[Group=="Donor"],
        Value[Group=="Recipient"],
        paired = TRUE
      )$statistic,
      p_value   = wilcox.test(
        Value[Group=="Donor"],
        Value[Group=="Recipient"],
        paired = TRUE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  print(wilcox_results)
  # one‐sample tests vs. zero for every Region × Group 
  one_sample_results <- df_plot %>%
    dplyr::group_by(Region, Group) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(Value, mu = 0)$statistic,
      p_value   = wilcox.test(Value, mu = 0)$p.value,
      label     = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      .groups = "drop"
    )
  print(one_sample_results)
  
  # 6) Build positions for one‐sample labels
  y_top  <- if (!is.na(y_limits[2])) y_limits[2] else max(df_plot$Value * 1e4, na.rm = TRUE)
  text_y <- y_top * 0.98  # 98% of the top
  
  one_sample_df <- one_sample_results %>%
    filter(label != "ns") %>%
    mutate(
      x0 = as.numeric(Region),          # 1 = Antigenic, 2 = Non-antigenic
      x  = x0 + case_when(              # offsets that match the dodge
        Group == "Donor"     ~ -0.27,   # left
        Group == "Recipient" ~  0.00,   # centre
        Group == "All" ~  0.27    # right
      ),
      y  = text_y
    )
  
  # 7) draw the plot
  exponent <- -4
  p <- ggplot(df_plot, aes(x = Region, y = Value * 1e4, fill = Group)) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_quasirandom(
      aes(color = Group, group = Group),
      width       = 0.1,
      dodge.width = 0.8,
      size        = quasi_size,
      alpha       = quasi_alpha,
      shape       = 16,
      stroke      = 0
    ) +
    geom_boxplot(
      position      = position_dodge(width = 0.8),
      width         = 0.6,
      outlier.shape = NA,
      alpha         = boxplot_alpha
    ) +
    geom_text(
      data        = one_sample_df,
      mapping     = aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size        = 4,
      vjust       = 1) +
    scale_fill_manual(values = c("Donor"="black", "Recipient"="red", "All"="steelblue2")) +
    scale_color_manual(values = c("Donor"="black", "Recipient"="red", "All"="steelblue2")) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = title,
      x     = "",
      y     = bquote("\u03c0N-\u03c0S per site (×10"^.(exponent)*")")
    ) +
    theme(
      legend.title    = element_blank(),
      legend.key      = element_blank(),
      legend.position = "none"
    ) +
    axis_formatting + legend_formatting + background_formatting
  
  return(p)
}
# plot the data
plot_piN_minus_piS_bw_H3N2_1718_noZeros <- plot_antigenic_piN_minus_piS_nozeros("17-18_H3N2", y_limits = c(-30, 20), title = "2017–18 H3N2")
plot_piN_minus_piS_bw_H3N2_1819_noZeros <- plot_antigenic_piN_minus_piS_nozeros("18-19_H3N2", y_limits = c(-30, 20), title = "2018–19 H3N2")
plot_piN_minus_piS_bw_H1N1_1718_noZeros <- plot_antigenic_piN_minus_piS_nozeros("17-18_H1N1", y_limits = c(-30, 20), title = "2017–18 H1N1")
plot_piN_minus_piS_bw_H1N1_1819_noZeros <- plot_antigenic_piN_minus_piS_nozeros("18-19_H1N1", y_limits = c(-30, 20), title = "2018–19 H1N1")

# clean the plots
p1 <- plot_piN_minus_piS_bw_H3N2_1718_noZeros + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
p2 <- plot_piN_minus_piS_bw_H3N2_1819_noZeros + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
p3 <- plot_piN_minus_piS_bw_H1N1_1718_noZeros + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
p4 <- plot_piN_minus_piS_bw_H1N1_1819_noZeros + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
# make it a figure
plot2_piN_minus_piS_bw_all_noZeros <- plot_grid(p1, p2, p3, p4,
                                                labels = c("A", "B", "C", "D"), 
                                                label_size = Size_adjust, 
                                                hjust = LR_adjust, vjust = UD_adjust, 
                                                ncol = 2, nrow = 2, rel_heights = c(1, 1, 1, 1))

## create the function
plot_antigenic_piN_minus_piS <- function(season_subtype = "",
                                         plot_title     = season_subtype,
                                         df_wh          = df_AA_HA, 
                                         df_bw          = df_AA_HA_paired,
                                         y_limits       = c(NA, NA)) {
  # --- subset by season or all
  if (is.null(season_subtype) || season_subtype == "") {
    df_wh_sub <- df_wh
    df_bw_sub <- df_bw
    plot_title <- plot_title %||% "All seasons"
  } else {
    df_wh_sub <- df_wh %>% filter(Merge_ID == season_subtype)
    df_bw_sub <- df_bw %>% filter(Merge_ID == season_subtype)
    plot_title <- plot_title %||% season_subtype
  }
  
  # --- All (within‐host)
  df_comm <- df_wh_sub %>%
    transmute(
      Region = if_else(Antigenicity == "Antigenic", "Antigenic", "Non-antigenic"),
      Group  = "All",
      Value  = perSite_piN_minus_piS
    )
  
  # --- donor/recipient long
  df_bw_long <- df_bw_sub %>%
    filter(Gene %in% c("Antigenic_all", "Non-antigenic_all")) %>%
    transmute(
      Region    = recode(Gene,
                         "Antigenic_all"     = "Antigenic",
                         "Non-antigenic_all" = "Non-antigenic"),
      Donor     = Donor_perSite_piN_minus_piS,
      Recipient = Recipient_perSite_piN_minus_piS
    ) %>%
    pivot_longer(c(Donor, Recipient),
                 names_to  = "Group",
                 values_to = "Value")
  
  # --- combine & factor
  df_plot <- bind_rows(df_comm, df_bw_long) %>%
    mutate(
      Region = factor(Region, levels = c("Antigenic", "Non-antigenic")),
      Group  = factor(Group,  levels = c("Donor", "Recipient", "All"))
    )
  
  # --- 1) Medians
  medians <- df_plot %>%
    dplyr::group_by(Region, Group) %>%
    dplyr::summarise(
      median_per10k = median(Value * 1e4, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- 2) Region‐wise comparisons with mixed pairing/unpaired logic
  region_comparisons <- df_plot %>%
    dplyr::group_by(Region) %>%
    dplyr::group_modify(~ {
      df   <- .x
      donor    <- df$Value[df$Group == "Donor"]
      recipient<- df$Value[df$Group == "Recipient"]
      All<- df$Value[df$Group == "All"]
      
      tibble::tibble(
        comparison = c("Donor vs Recipient",
                       "Donor vs All",
                       "Recipient vs All"),
        statistic = c(
          wilcox.test(donor,     recipient, paired   = TRUE)$statistic,
          wilcox.test(donor,     All, paired   = FALSE)$statistic,
          wilcox.test(recipient, All, paired   = FALSE)$statistic
        ),
        p_value   = c(
          wilcox.test(donor,     recipient, paired   = TRUE)$p.value,
          wilcox.test(donor,     All, paired   = FALSE)$p.value,
          wilcox.test(recipient, All, paired   = FALSE)$p.value
        )
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      label = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  # --- X) Region comparison within each Group (antigenic vs non-antigenic)
  region_by_group <- df_plot %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      # pull out the two sets as lists
      a         = list(Value[Region == "Antigenic"]),
      n         = list(Value[Region == "Non-antigenic"]),
      # decide pairing once per group
      paired    = first(Group) %in% c("Donor", "Recipient"),
      # run the test with that scalar paired flag
      statistic = wilcox.test(
        a[[1]],
        n[[1]],
        paired = first(Group) %in% c("Donor", "Recipient")
      )$statistic,
      p_value   = wilcox.test(
        a[[1]],
        n[[1]],
        paired = first(Group) %in% c("Donor", "Recipient")
      )$p.value,
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  print(region_by_group)
  
  
  
  # --- 3) Normality
  normality_results <- df_plot %>%
    dplyr::group_by(Region, Group) %>%
    dplyr::summarise(
      n         = n(),
      W         = if (length(unique(Value)) == 1) NA_real_ else shapiro.test(Value)$statistic,
      p_value   = if (length(unique(Value)) == 1) NA_real_ else shapiro.test(Value)$p.value,
      is_normal = if (length(unique(Value)) == 1) FALSE else shapiro.test(Value)$p.value > 0.05,
      .groups   = "drop"
    )
  
  #normality_results <- df_plot %>%
  #  dplyr::group_by(Region, Group) %>%
  #  dplyr::summarise(
  #    n         = n(),
  #    W         = shapiro.test(Value)$statistic,
  #    p_value   = shapiro.test(Value)$p.value,
  #    is_normal = p_value > 0.05,
  #    .groups   = "drop"
  #  )
  print(normality_results)
  
  # --- 4) Donor vs Recipient (paired) per Region
  wilcox_results <- df_bw_long %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(
        Value[Group == "Donor"],
        Value[Group == "Recipient"],
        paired = TRUE
      )$statistic,
      p_value   = wilcox.test(
        Value[Group == "Donor"],
        Value[Group == "Recipient"],
        paired = TRUE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))

  # --- 5) One-sample vs zero for all Region × Group
  one_sample_results <- df_plot %>%
    dplyr::group_by(Region, Group) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(Value, mu = 0)$statistic,
      p_value   = wilcox.test(Value, mu = 0)$p.value,
      label     = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      .groups = "drop"
    )

  # --- 6) Build label positions as before ---
  y_top  <- if (!is.na(y_limits[2])) y_limits[2] else max(df_plot$Value * 1e4, na.rm = TRUE)
  text_y <- y_top * 0.98
  one_sample_df <- one_sample_results %>%
    filter(label != "ns") %>%
    mutate(
      x0 = as.numeric(Region),
      x  = x0 + case_when(
        Group == "Donor"     ~ -0.27,
        Group == "Recipient" ~  0.00,
        Group == "All" ~  0.27
      ),
      y  = text_y
    )
  
  # --- 7) Plot
  exponent <- -4
  p <- ggplot(df_plot, aes(x = Region, y = Value * 1e4, fill = Group)) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_quasirandom(aes(color = Group, group = Group),
                     width       = 0.1,
                     dodge.width = 0.8,
                     size        = quasi_size+.5,
                     alpha       = quasi_alpha,
                     shape       = 16,
                     stroke      = 0) +
    geom_boxplot(position      = position_dodge(0.8),
                 width         = 0.6,
                 outlier.shape = NA,
                 alpha         = boxplot_alpha) +
    geom_text(data        = one_sample_df,
              mapping     = aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size        = 4,
              vjust       = 1) +
    scale_fill_manual(values = c("Donor"="black",
                                 "Recipient"="red",
                                 "All"="steelblue2")) +
    scale_color_manual(values = c("Donor"="black",
                                  "Recipient"="red",
                                  "All"="steelblue2")) +
    coord_cartesian(ylim = y_limits) +
    labs(title = plot_title,
         x     = NULL,
         y     = bquote("\u03c0N-\u03c0S per site (×10"^.(exponent)*")")) +
    theme(legend.position = "none") +
    axis_formatting + legend_formatting + background_formatting
  
  # --- Distribution breakdown for "All" group: negative, zero, positive
  value_distribution <- df_plot %>%
    filter(Group == "All") %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise(
      n_total  = n(),
      n_neg    = sum(Value < 0, na.rm = TRUE),
      n_zero   = sum(Value == 0, na.rm = TRUE),
      n_pos    = sum(Value > 0, na.rm = TRUE),
      pct_neg  = 100 * n_neg / n_total,
      pct_zero = 100 * n_zero / n_total,
      pct_pos  = 100 * n_pos / n_total,
      .groups  = "drop"
    )
  
  print(value_distribution)
  
  # --- return everything ---
  list(
    plot                  = p,
    medians               = medians,
    region_comparisons    = region_comparisons,
    donor_recipient_tests = wilcox_results,
    one_sample_tests      = one_sample_results
  )
}


## plot the data –
plot_piN_minus_piS_bw           <- plot_antigenic_piN_minus_piS("", plot_title = "", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_H3N2_1718 <- plot_antigenic_piN_minus_piS("17-18_H3N2", plot_title = "2017–18 H3N2", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_H3N2_1819 <- plot_antigenic_piN_minus_piS("18-19_H3N2", plot_title = "2018–19 H3N2", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_H1N1_1718 <- plot_antigenic_piN_minus_piS("17-18_H1N1", plot_title = "2017–18 H1N1", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_H1N1_1819 <- plot_antigenic_piN_minus_piS("18-19_H1N1", plot_title = "2018–19 H1N1", y_limits = c(-30, 20))
# clean the plots
p1 <- plot_piN_minus_piS_bw_H3N2_1718$plot + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
p2 <- plot_piN_minus_piS_bw_H3N2_1819$plot + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
p3 <- plot_piN_minus_piS_bw_H1N1_1718$plot + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
p4 <- plot_piN_minus_piS_bw_H1N1_1819$plot + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"))
# make it a figure
plot2_piN_minus_piS_bw_all <- plot_grid(p1, p2, p3, p4,
                                        labels = c("A", "B", "C", "D"), 
                                        label_size = Size_adjust, 
                                        hjust = LR_adjust, vjust = UD_adjust, 
                                        ncol = 2, nrow = 2, rel_heights = c(1, 1, 1, 1))

### plotting πN-πS per site for donors, recipients, and All samples (each antigenic region)
## create the function
plot_piN_minus_piS_by_site <- function(season_subtype,
                                       plot_title = season_subtype,
                                       df_wh   = df_AA_HA,
                                       df_bw   = df_AA_HA_paired,
                                       y_limits = c(NA, NA)) {
  # 1) All ("within-host") by site, dropping the “all” row
  df_comm <- df_wh %>%
    filter(Merge_ID == season_subtype,
           Antigenic_Site != "all") %>%
    transmute(
      Site  = Antigenic_Site,
      Group = "All",
      Value = perSite_piN_minus_piS
    )
  
  # 2) donor vs recipient by site, dropping the “all” row
  df_bw_long <- df_bw %>%
    filter(Merge_ID == season_subtype,
           grepl("^Antigenic_", Gene),
           Gene != "Antigenic_all") %>%
    transmute(
      Site      = sub("^Antigenic_", "", Gene),
      Donor     = Donor_perSite_piN_minus_piS,
      Recipient = Recipient_perSite_piN_minus_piS
    ) %>%
    pivot_longer(
      cols      = c(Donor, Recipient),
      names_to  = "Group",
      values_to = "Value"
    )
  
  # 3) combine & set factor levels so sites appear in order
  site_levels <- unique(df_comm$Site)
  df_plot <- bind_rows(df_comm, df_bw_long) %>%
    mutate(
      Site  = factor(Site,  levels = site_levels),
      Group = factor(Group, levels = c("Donor", "Recipient", "All"))
    )
  
  # 4) Normality tests by Site × Group
  normality_results <- df_plot %>%
    group_by(Site, Group) %>%
    dplyr::summarise(
      n        = n(),
      W        = if (length(Value) >= 3 && length(unique(Value)) > 1) {
        shapiro.test(Value)$statistic
      } else {
        NA_real_
      },
      p_value  = if (length(Value) >= 3 && length(unique(Value)) > 1) {
        shapiro.test(Value)$p.value
      } else {
        NA_real_
      },
      is_normal = ifelse(is.na(p_value), FALSE, p_value > 0.05),
      .groups  = "drop"
    )
  print(normality_results)
  
  # 5) Paired Wilcoxon (Donor vs Recipient) *within each Site*
  wilcox_results <- df_bw_long %>%
    group_by(Site) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(
        Value[Group == "Donor"],
        Value[Group == "Recipient"],
        paired = TRUE
      )$statistic,
      p_value   = wilcox.test(
        Value[Group == "Donor"],
        Value[Group == "Recipient"],
        paired = TRUE
      )$p.value,
      .groups   = "drop"
    ) %>%
    mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  print(wilcox_results)
  
  # 6) One‐sample Wilcoxon vs zero for every Site × Group
  one_sample_results <- df_plot %>%
    group_by(Site, Group) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(Value, mu = 0)$statistic,
      p_value   = wilcox.test(Value, mu = 0)$p.value,
      label     = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      .groups   = "drop"
    )
  print(one_sample_results)
  
  # 7) Build positions for one‐sample labels
  y_top  <- if (!is.na(y_limits[2])) y_limits[2] else max(df_plot$Value * 1e4, na.rm = TRUE)
  text_y <- y_top * 0.98  # 98% of the top
  
  one_sample_df <- one_sample_results %>%
    filter(label != "ns") %>%
    mutate(
      x0 = as.numeric(Site),            # 1 = Antigenic, 2 = Non-antigenic
      x  = x0 + case_when(              # offsets that match the dodge
        Group == "Donor"     ~ -0.27,   # left
        Group == "Recipient" ~  0.00,   # centre
        Group == "All" ~  0.27    # right
      ),
      y  = text_y
    )
  
  # 8) Finally, draw it
  exponent <- -4
  p <- ggplot(df_plot, aes(x = Site, y = Value * 1e4, fill = Group)) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_quasirandom(aes(color = Group, group = Group),
                     width       = 0.1,
                     dodge.width = 0.8,
                     size        = quasi_size,
                     alpha       = quasi_alpha,
                     shape       = 16,
                     stroke      = 0) +
    geom_boxplot(position      = position_dodge(width = 0.8),
                 width         = 0.6,
                 outlier.shape = NA,
                 alpha         = boxplot_alpha) +
    geom_text(data        = one_sample_df,
              mapping     = aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              size        = 4,
              vjust       = 1) + 
    scale_fill_manual(values = c("Donor"="black", "Recipient"="red", "All"="steelblue2")) +
    scale_color_manual(values = c("Donor"="black", "Recipient"="red", "All"="steelblue2")) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = plot_title,
      x     = "",
      y     = bquote("\u03c0N-\u03c0S per site (×10"^.(exponent)*")")
    ) +
    theme(
      legend.title    = element_blank(),
      legend.key      = element_blank(),
      legend.position = "none",
      axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))
      ) +
    axis_formatting + legend_formatting + background_formatting
  return(p)
}

## plot the data
plot_piN_minus_piS_bw_by_site_H3N2_1718 <- plot_piN_minus_piS_by_site("17-18_H3N2", "2017–18 H3N2", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_by_site_H3N2_1819 <- plot_piN_minus_piS_by_site("18-19_H3N2", "2018–19 H3N2", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_by_site_H1N1_1718 <- plot_piN_minus_piS_by_site("17-18_H1N1", "2017–18 H1N1", y_limits = c(-30, 20))
plot_piN_minus_piS_bw_by_site_H1N1_1819 <- plot_piN_minus_piS_by_site("18-19_H1N1", "2018–19 H1N1", y_limits = c(-30, 20))
# clean the plots
p1 <- plot_piN_minus_piS_bw_by_site_H3N2_1718 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p2 <- plot_piN_minus_piS_bw_by_site_H3N2_1819 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p3 <- plot_piN_minus_piS_bw_by_site_H1N1_1718 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p4 <- plot_piN_minus_piS_bw_by_site_H1N1_1819 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
# make it a figure
plot_piN_minus_piS_bw_by_site <- plot_grid(p1, p2, p3, p4,
                                           labels = c("A", "B", "C", "D"), 
                                           label_size = Size_adjust, 
                                           hjust = LR_adjust, vjust = UD_adjust, 
                                           ncol = 1, nrow = 4, rel_heights = c(1, 1, 1, 1))

#### Figure 1a  - Plot: (whole-genome) iSNVs across the genome ####
#df_all_iSNVs$Merge_ID <- as.factor(df_all_iSNVs$Merge_ID); summary(df_all_iSNVs); plot(table(df_all_iSNVs$WSLH_ID)); table <- as.data.frame(table(df_all_iSNVs$WSLH_ID)); max(table$Freq)
#df_sub50AF$Merge_ID <- as.factor(df_sub50AF$Merge_ID); summary(df_sub50AF); plot(table(df_sub50AF$WSLH_ID)); table <- as.data.frame(table(df_sub50AF$WSLH_ID)); max(table$Freq)

## Make NA not the problematic kind of NA
df_sub50AF$Gene_ID[is.na(df_sub50AF$Gene_ID)] <- "NA"
df_sub50AF$Gene_Name[is.na(df_sub50AF$Gene_Name)] <- "NA"

## Only include some genes
df_sub50AF <- df_sub50AF[df_sub50AF$Gene_Name %in% factor_by_gene,]
# change factor level order
df_sub50AF$Gene_Name <- factor(df_sub50AF$Gene_Name, levels = factor_by_gene)

## change factor level order for Ann
df_sub50AF$Ann <- factor(df_sub50AF$Ann, levels = c("Nonsynonymous", "Synonymous", "Stop gained"))
table <- as.data.frame(table(df_all_iSNVs$Ann[df_all_iSNVs$Merge_ID=="18-19_H3N2"]))

## Fig1a: Manhattan plot
function_manhattan_plot <- function(year_subtype, title, legend_position_PA) {
  # define plot groups and ranges
  genes <- factor_by_gene
  x_axis_start <- c(0,0,0,0,0,0)
  x_axis_end <- c(2316,2316,2209,1737,1541,1441)
  by <- c(500,500,500,500,500,500)
  spec <- data.frame(genes, x_axis_start, x_axis_end, by)
  # create a list to put the gene plots in 
  plots <- list()
  for (gene in genes) {
    # filter to selected year_subtype
    if (year_subtype=="17-19_H3N2") { # for combined H3 group
      df_plot     <- df_sub50AF[df_sub50AF$Gene_Name==gene & df_sub50AF$Merge_ID=="17-18_H3N2" | df_sub50AF$Merge_ID=="18-19_H3N2",]
      }
    if (year_subtype!="17-19_H3N2") {
      df_plot     <- df_sub50AF[df_sub50AF$Gene_Name==gene & df_sub50AF$Merge_ID==year_subtype,]
      }
    x_axis_start<- spec$x_axis_start[spec$genes==gene]
    x_axis_plot <- spec$x_axis_end[spec$genes==gene]
    by_plot     <- spec$by[spec$genes==gene]
    name_plot   <- paste(gene,"_","plot")
    name_plot   <- gsub(" ","",name_plot)
  
    ## y-axis formatting
    # y-axis line and title
    if (gene == "PB2"| gene == "HA"){
      y_aesthetics = theme(axis.line.y = element_line(colour="grey"),
                           axis.text.y = element_text(hjust=0.5, size = 6),
                           axis.title.y = element_text(size = 6, margin = margin(r = 4)))
      } else {
        y_aesthetics = theme(axis.line.y=element_blank(),
                             axis.ticks.y= element_blank(),
                             axis.text.y=element_blank(),
                             axis.title.y = element_blank())
      }
    # x-axis title
    if (gene == "PB1"){
      y_aesthetics_x = theme(axis.title.x = element_text(size = 6))
      x_label <- "Nucleotide position"
    } else {
      y_aesthetics_x = theme(axis.title.x = element_text(size = 6))
      x_label <- ""
    }
    
    ## legend formatting
    if (gene == "PA"){
      legend_aesthetics = theme(legend.background = element_blank(),
                                legend.title = element_blank(),
                                legend.key = element_blank(),
                                legend.position = legend_position_PA)
      } else {
        legend_aesthetics = theme(legend.background = element_blank(),
                                  legend.title = element_blank(),
                                  legend.key = element_blank(),
                                  legend.position = "none")
      }
    # plot
    plot <- ggplot(df_plot, aes(x = STA, y = AF, color = Ann)) + 
      geom_point(size = 1.5, alpha = 0.75, stroke = NA) + 
      geom_hline(yintercept = 3, linetype = "solid", color = "grey") + 
      coord_cartesian(ylim = c(0, 50), xlim = c(x_axis_start, x_axis_plot)) + 
      labs(y = "iSNV frequency (%)", x = x_label, title=gene) + 
      scale_x_continuous(breaks = seq(signif(x_axis_start,3),
                                      signif(x_axis_plot,3),by_plot),
                         labels = seq(signif(x_axis_start,3),
                                      signif(x_axis_plot,3),by_plot)) + 
      scale_color_manual(values = palette_muts_NS_S_Stopgained) + 
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 6, vjust = -2),
            axis.line.x = element_line(color="grey"),
            axis.text.x = element_text(size = 6)) + 
      guides(color = guide_legend(override.aes=list(shape = 15, size = 5))) + 
      y_aesthetics + y_aesthetics_x + legend_aesthetics + legend_formatting
    plots[[name_plot]] <- plot
    }
    if (title == "") {
      Fig2a_top <- plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow = 1, rel_widths = c(1.1, 1, 1.6))
      Fig2a_bot <- plot_grid(plots[[4]], plots[[5]], plots[[6]], nrow = 1, rel_widths = c(1.5, 1.5, 1))
      Fig1a <- plot_grid(Fig2a_top, Fig2a_bot, nrow = 2, ncol = 1)
      rm(Fig2a_top)
      rm(Fig2a_bot)
      return(Fig1a)
    }
    if (title != "") {
      Fig2a_top <- plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow = 1, rel_widths = c(1.2, 1, 1))
      Fig2a_bot <- plot_grid(plots[[4]], plots[[5]], plots[[6]], nrow = 1, rel_widths = c(1.75, 1.5, 1))
      Fig1a <- plot_grid(Fig2a_top, Fig2a_bot, nrow = 2, ncol = 1)
      rm(Fig2a_top)
      rm(Fig2a_bot)
      title <- ggdraw() + draw_label(title, fontface='bold', size = 10, vjust = 1.25)
      Fig1a <- plot_grid(title, Fig1a, ncol=1, rel_heights=c(0.1, 1))
      return(Fig1a)
      }
    }
plot_manhattan <- function_manhattan_plot("17-19_H3N2", "", "right")
plot_manhattan_1718_H3N2 <- function_manhattan_plot("17-18_H3N2", "2017–18 H3N2", "none")
plot_manhattan_1819_H3N2 <- function_manhattan_plot("18-19_H3N2", "2018–19 H3N2", "none")
plot_manhattan_1718_H1N1 <- function_manhattan_plot("17-18_H1N1", "2017–18 H1N1", "none")
plot_manhattan_1819_H1N1 <- function_manhattan_plot("18-19_H1N1", "2018–19 H1N1", "none")

## FigS1a
# Prepare the subset of interest: just vaccinated and unvaccinated
df_plot <- df_sub50AF_md_n_SNVs_by_vax %>%
  filter(flu_vaccine %in% c("Vaccinated", "Unvaccinated")) %>%
  mutate(
    x_axis = flu_vaccine,
    Merge_ID = recode(Merge_ID,
                      "17-18_H3N2" = "2017–18 H3N2",
                      "18-19_H3N2" = "2018–19 H3N2",
                      "17-18_H1N1" = "2017–18 H1N1",
                      "18-19_H1N1" = "2018–19 H1N1"
    )
  )

df_plot$Group <- df_plot$flu_vaccine


# plot
FigS1a_dist <- ggplot(df_plot, aes(x = Count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "white") +
  scale_x_continuous(
    breaks = seq(0, max(df_plot$Count, na.rm = TRUE) + 1, 1),
    labels = function(x) ifelse(x %% 5 == 0, as.character(x), ""),
    expand = c(0.01, 0.01)
  ) +
  labs(x = "# of iSNVs per sample", y = "# of samples") +
  coord_cartesian(xlim = c(0,max(df_plot$Count))) + 
  axis_formatting +
  legend_formatting +
  background_formatting +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25)
  )


# Create a new 'Either' group by combining vax + unvax per subtype
df_all_group <- df_plot %>%
  dplyr::filter(flu_vaccine %in% c("Vaccinated", "Unvaccinated")) %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::mutate(Group = "Either") %>%
  dplyr::ungroup()

# Combine with the original data
df_plot_combined <- bind_rows(df_plot, df_all_group)

# Set factor levels so 'Either' appears last
df_plot_combined$Group <- factor(df_plot_combined$Group, levels = c("Unvaccinated", "Vaccinated", "Either"))

# Assess normality
normality_results <- df_plot_combined %>%
  dplyr::group_by(Merge_ID, x_axis) %>%
  dplyr::summarise(
    n = n(),
    W = shapiro.test(Count)$statistic,
    p_value = shapiro.test(Count)$p.value,
    is_normal = p_value > 0.05,
    .groups = "drop"
  )

# Perform stats
wilcox_results <- df_plot_combined %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::filter(all(c("Vaccinated", "Unvaccinated") %in% unique(x_axis))) %>%
  dplyr::summarise(
    n = n(),
    statistic = wilcox.test(Count ~ x_axis)$statistic,
    p_value = wilcox.test(Count ~ x_axis)$p.value,
    .groups = "drop"
  ) %>%
  mutate(significant = p_value < 0.05,
         label = case_when(
           p_value < 0.001 ~ "***",
           p_value < 0.01  ~ "**",
           p_value < 0.05  ~ "*",
           TRUE            ~ "ns"
         ))

# Create stats for plot
box_heights <- df_plot_combined %>%
  dplyr::group_by(Merge_ID, x_axis) %>%
  dplyr::summarise(box_top = quantile(Count, 0.75) + 1.5 * IQR(Count), .groups = "drop")

max_box_heights <- box_heights %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::summarise(y_pos = max(box_top, na.rm = TRUE) * 1.05, .groups = "drop")

sig_labels <- wilcox_results %>%
  dplyr::filter(label != "ns") %>%
  dplyr::left_join(max_box_heights, by = "Merge_ID") %>%
  dplyr::mutate(x = as.numeric(factor(Merge_ID, levels = unique(df_plot$Merge_ID))))

sig_labels <- sig_labels %>%
  mutate(
    x_center = as.numeric(factor(Merge_ID)),
    x_start = x_center - 0.2,
    x_end   = x_center + 0.2
  )

## Plot standards
ymax <- 200
tick_vals <- seq(0, ymax, 1)

df_plot_combined$Merge_ID <- factor(df_plot_combined$Merge_ID, 
                                    levels = c("2017–18 H3N2",
                                               "2018–19 H3N2",
                                               "2017–18 H1N1",
                                               "2018–19 H1N1"))

## brief calculations
med_by_group <- df_plot_combined %>%
  # 1. season × existing groups
  dplyr::group_by(Merge_ID, Group) %>%
  dplyr::summarize(
    median_Count = median(Count, na.rm = TRUE),
    .groups = "drop"
  )

med_either     <- df_plot_combined %>%
  # 2. season × “Either”
  dplyr::group_by(Merge_ID) %>%
  dplyr::summarize(
    median_Count = median(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = "Either")

med_overall_vg <- df_plot_combined %>%
  # 3a. overall × Vaccinated & Unvaccinated
  dplyr::group_by(Group) %>%
  dplyr::summarize(
    median_Count = median(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Merge_ID = "All seasons")

med_overall_e  <- df_plot_combined %>%
  # 3b. overall × Either
  dplyr::summarize(
    median_Count = median(Count, na.rm = TRUE)
  ) %>%
  mutate(
    Group    = "Either",
    Merge_ID = "All seasons"
  )

# 4. bind and order
medians <- bind_rows(
  med_by_group,
  med_either,
  med_overall_vg,
  med_overall_e
) %>%
  arrange(
    factor(Merge_ID,
           levels = c("2017–18 H3N2",
                      "2018–19 H3N2",
                      "2017–18 H1N1",
                      "2018–19 H1N1",
                      "All seasons")),
    factor(Group, levels = c("Vaccinated", "Unvaccinated", "Either"))
  )

medians


## Plot
FigS1b_vax <- ggplot(df_plot_combined, aes(x = Merge_ID, y = Count, fill = Group)) +
  geom_quasirandom(
    aes(color = Group),
    width       = 0.1,
    dodge.width = 0.8,
    size        = quasi_size,
    alpha       = quasi_alpha,
    shape       = 16,
    stroke      = 0
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), 
               width = 0.6, 
               alpha = 0.6) +
  scale_fill_manual(values = palette_vax) +
  scale_color_manual(values = palette_vax) + 
  #scale_y_continuous(
  #  breaks = tick_vals,
  #  labels = ifelse(tick_vals %% 5 == 0, as.character(tick_vals), ""),
  #  expand = c(0.01, 0.01)) + 
  scale_x_discrete(labels = c("2017–18 H3N2" = "2017–18\nH3N2", 
                              "2018–19 H3N2" = "2018–19\nH3N2", 
                              "2017–18 H1N1" = "2017–18\nH1N1", 
                              "2018–19 H1N1" = "2018–19\nH1N1")) + 
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1, base = 10),
    limits = c(0, ymax),
    breaks = c(0, 1, 10, 100),
    minor_breaks = c(2, 3, 4, 5, 6, 7, 8, 9, 20, 30, 40, 50, 60, 70, 80, 90)
  ) +
  geom_segment(data = sig_labels,
               aes(x = x_start, xend = x_end,
                   y = y_pos, yend = y_pos),
               inherit.aes = FALSE,
               linewidth = 0.5, color = "black") +
  geom_text(data = sig_labels,
            aes(x = x_center, y = y_pos + 0.01 * max(df_plot$Count, na.rm = TRUE), label = label),
            inherit.aes = FALSE,
            size = 4, vjust = 0) + 
  coord_cartesian(ylim = c(0, ymax)) + 
  labs(x = "", y = "# iSNVs per sample") +
  axis_formatting +
  legend_formatting +
  background_formatting +
  theme(
    legend.title = element_blank(), 
    legend.key = element_blank(),
    legend.position = c(0.825, 0.8),  # Centered horizontally, just above the plot
    legend.direction = "vertical",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    panel.grid.major.y = element_line(color = "grey70", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "grey90", linewidth = 0.25)
    )

## FigS1b_pi
# Step 1: Total π per sample
df_total_pi <- aggregate(df_divergence$pi, list(df_divergence$WSLH_ID), FUN = sum)

# Step 2: Annotate with vax/subtype
df_total_pi_annotated <- df_total_pi %>%
  dplyr::rename(Sample_ID = Group.1, Count = x) %>%
  dplyr::inner_join(df_sub50AF_md_n_SNVs_by_vax, by = c("Sample_ID" = "WSLH_ID")) %>%
  dplyr::filter(flu_vaccine %in% c("Vaccinated", "Unvaccinated")) %>%
  dplyr::mutate(
    Group = flu_vaccine,
    Merge_ID = recode(Merge_ID,
                      "17-18_H3N2" = "2017–18 H3N2",
                      "18-19_H3N2" = "2018–19 H3N2",
                      "17-18_H1N1" = "2017–18 H1N1",
                      "18-19_H1N1" = "2018–19 H1N1")
  )

# Step 3: Add "Either" group
df_pi_either <- df_total_pi_annotated %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::mutate(Group = "Either") %>%
  dplyr::ungroup()

df_pi_plot <- bind_rows(df_total_pi_annotated, df_pi_either)
df_pi_plot$Group <- factor(df_pi_plot$Group, levels = c("Unvaccinated", "Vaccinated", "Either"))

# Step 4: Order factors and scale π
scaling_factor <- 1e4
exponent <- -log10(scaling_factor)
df_pi_plot$Count_scaled <- df_pi_plot$Count.x * scaling_factor
y_height <- 200

# Step 5: Stats (exclude "Either")
normality_results <- df_pi_plot %>%
  dplyr::filter(Group != "Either") %>%
  dplyr::group_by(Merge_ID, Group) %>%
  dplyr::summarise(
    n = n(),
    W = shapiro.test(Count_scaled)$statistic,
    p_value = shapiro.test(Count_scaled)$p.value,
    is_normal = p_value > 0.05,
    .groups = "drop"
  )

wilcox_results <- df_pi_plot %>%
  dplyr::filter(Group %in% c("Vaccinated", "Unvaccinated")) %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::filter(all(c("Vaccinated", "Unvaccinated") %in% unique(Group))) %>%
  dplyr::summarise(
    n = n(),
    statistic = wilcox.test(Count_scaled ~ Group)$statistic,
    p_value = wilcox.test(Count_scaled ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# Step 6: Significance label positions
box_heights <- df_pi_plot %>%
  dplyr::filter(Group %in% c("Vaccinated", "Unvaccinated")) %>%
  dplyr::group_by(Merge_ID, Group) %>%
  dplyr::summarise(box_top = quantile(Count_scaled, 0.75) + 1.5 * IQR(Count_scaled), .groups = "drop")

max_box_heights <- box_heights %>%
  dplyr::group_by(Merge_ID) %>%
  dplyr::summarise(y_pos = max(box_top, na.rm = TRUE) * 1.05, .groups = "drop")

sig_labels <- wilcox_results %>%
  dplyr::filter(label != "ns") %>%
  dplyr::left_join(max_box_heights, by = "Merge_ID") %>%
  dplyr::mutate(x = as.numeric(factor(Merge_ID, levels = levels(df_pi_plot$Merge_ID))))

# Step 7: Final plot
df_pi_plot$Merge_ID <- factor(
  df_pi_plot$Merge_ID,
  levels = c("2017–18 H3N2", "2018–19 H3N2", "2017–18 H1N1", "2018–19 H1N1")
)

FigS1b_pi <- ggplot(df_pi_plot, aes(x = Merge_ID, y = Count_scaled)) +
  geom_boxplot(aes(fill = Group), outlier.shape = NA, position = position_dodge(width = 0.8), width = 0.6) +
  geom_jitter(aes(color = Group),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1) +
  scale_fill_manual(values = palette_vax) +
  scale_color_manual(values = palette_vax) +
  scale_y_continuous(limits = c(0, y_height)) +
  labs(x = "", y = bquote("\u03c0 (×10"^.(exponent) * ")")) +
  geom_segment(data = sig_labels,
               aes(x = x - 0.2, xend = x + 0.2, y = y_pos, yend = y_pos),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_text(data = sig_labels,
            aes(x = x, y = y_pos + 0.01 * y_height, label = label),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "right") +
  axis_formatting +
  legend_formatting +
  background_formatting

## Print stats
#print(normality_results)
#print(wilcox_results)
#
## Return plot
#FigS1b_pi


#### Figure 1bc - Plot: (whole-genome) iSNV frequency spectra ####
## df for neutral expectation
df_mut_bins_neutral <- df_mut_bins_prop[df_mut_bins_prop$Mutation_type=="Neutral expectation" & 
                                          df_mut_bins_prop$Merge_ID=="18-19_H1N1",] #there's nothing special about 18-19_H1N1 here
df_mut_bins_neutral$X.1 <- NULL
df_mut_bins_neutral$X <- NULL
df_mut_bins_neutral$Total <- ""
df_mut_bins_neutral$Merge_ID <- ""

## factor order
df_mut_bins_prop_noNeut$Mutation_type <- factor(df_mut_bins_prop_noNeut$Mutation_type,
                                                levels = c("Nonsynonymous", "Synonymous", "Stop gained"))

## palette
palette_muts <- c("Nonsynonymous" = "#FF7F20",
                  "Synonymous" = "#4F7899",
                  "Stop gained" = "#7AD9C2")

# 1) compute your summaries once
df_summary <- df_mut_bins_prop_all %>% 
  dplyr::group_by(Merge_ID, Bins, Mutation_type) %>% 
  dplyr::summarise(
    n         = n(),
    mean_prop = mean(Prop, na.rm=TRUE),
    sd_prop   = sd(  Prop, na.rm=TRUE),
    .groups   = "drop"
  )

# 2) pick the season you want, e.g. 2017-18 H3N2
season   <- "17-18_H3N2"
raw      <- df_mut_bins_prop_all %>% filter(Merge_ID == season)
summary4 <- df_summary            %>% filter(Merge_ID == season)

# 3) make sure Bins is a factor in the order you like
bin_levels <- c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%")
raw$Bins               <- factor(raw$Bins,      levels=bin_levels)
summary4$Bins          <- factor(summary4$Bins, levels=bin_levels)
raw$Mutation_type      <- factor(raw$Mutation_type, levels=names(palette_muts))
summary4$Mutation_type <- factor(summary4$Mutation_type, levels=names(palette_muts))

# 4) the plot
plot_mut_bins_prop_1718_H3N2 <- ggplot(raw,
       aes(x = Bins, y = Prop, colour = Mutation_type, fill = Mutation_type,
           group = interaction(Bins, Mutation_type))) +
  # add in the neutral expectation dots and lines
  geom_point(data = df_mut_bins_neutral, 
             aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")),
                 y = Mean_prop), size = 2.5, stroke = NA) + 
  geom_line(data = df_mut_bins_neutral, 
            aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")), 
                y = Mean_prop, group = 1), linewidth = 1.0, alpha = quasi_alpha) +
  # raw proportions
  geom_quasirandom(
    width       = .2, 
    dodge.width = .8,
    size        = quasi_size, 
    alpha       = quasi_alpha,
    shape       = 21,
    stroke      = 0
  ) +
  # overlay ±1 SD bars
  geom_errorbar(
    data = summary4,
    inherit.aes = FALSE,
    aes(x    = Bins,
        ymin = mean_prop - sd_prop,
        ymax = mean_prop + sd_prop,
        group = Mutation_type),
    position = position_dodge(width = 0.8),
    width    = 0, size = 0.5
  ) +
  # overlay the means
  geom_point(
    data = summary4,
    inherit.aes = FALSE,              # drop x=Bins, y=Prop, etc.
    aes(x = Bins, y = mean_prop,
        fill = Mutation_type),
    position = position_dodge(width = 0.8),
    shape    = 23, size = 2, stroke = 0.5
  ) +
  # formatting
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = palette_muts) +
  scale_fill_manual(values  = palette_muts) +
  coord_cartesian(ylim = c(0,1)) + 
  labs(
    title = "2017–18 H3N2",
    x     = NULL,
    y     = "Mean proportion of variants"
  ) +
  axis_formatting +
  background_formatting +
  theme(legend.position = "none")



# 2) pick the season you want, e.g. 2017-18 H3N2
season   <- "18-19_H3N2"
raw      <- df_mut_bins_prop_all %>% filter(Merge_ID == season)
summary4 <- df_summary            %>% filter(Merge_ID == season)

# 3) make sure Bins is a factor in the order you like
bin_levels <- c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%")
raw$Bins               <- factor(raw$Bins,      levels=bin_levels)
summary4$Bins          <- factor(summary4$Bins, levels=bin_levels)
raw$Mutation_type      <- factor(raw$Mutation_type, levels=names(palette_muts))
summary4$Mutation_type <- factor(summary4$Mutation_type, levels=names(palette_muts))

# 4) the plot
plot_mut_bins_prop_1819_H3N2 <- ggplot(raw,
                                       aes(x = Bins, y = Prop, colour = Mutation_type, fill = Mutation_type,
                                           group = interaction(Bins, Mutation_type))) +
  # add in the neutral expectation dots and lines
  geom_point(data = df_mut_bins_neutral, 
             aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")),
                 y = Mean_prop), size = 2.5, stroke = NA) + 
  geom_line(data = df_mut_bins_neutral, 
            aes(x = factor(Bins, levels = c("3-10%", "10-20%", "20-30%", "30-40%", "40-50%")), 
                y = Mean_prop, group = 1), linewidth = 1.0, alpha = quasi_alpha) +
  # raw proportions
  geom_quasirandom(
    width       = .2, 
    dodge.width = .8,
    size        = quasi_size, 
    alpha       = quasi_alpha,
    shape       = 21,
    stroke      = 0
  ) +
  # overlay ±1 SD bars
  geom_errorbar(
    data = summary4,
    inherit.aes = FALSE,
    aes(x    = Bins,
        ymin = mean_prop - sd_prop,
        ymax = mean_prop + sd_prop,
        group = Mutation_type),
    position = position_dodge(width = 0.8),
    width    = 0, size = 0.5
  ) +
  # overlay the means
  geom_point(
    data = summary4,
    inherit.aes = FALSE,              # drop x=Bins, y=Prop, etc.
    aes(x = Bins, y = mean_prop,
        fill = Mutation_type),
    position = position_dodge(width = 0.8),
    shape    = 23, size = 2, stroke = 0.5
  ) +
  # formatting
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = palette_muts) +
  scale_fill_manual(values  = palette_muts) +
  coord_cartesian(ylim = c(0,1)) + 
  labs(
    title = "2018–19 H3N2",
    x     = NULL,
    y     = "Mean proportion of variants"
  ) +
  axis_formatting +
  background_formatting +
  theme(legend.position = "none")



## Fig2supp
function_Fig1suppY <- function(i, title) {
  df_plot <- df_all_iSNVs %>% 
    filter(Merge_ID == i)
  
  Fig1suppY <- ggplot(df_plot, aes(x = AF)) + 
    stat_bin(
      color = "black", fill = "white", binwidth = 0.1,
      aes(y = log10(cumsum(after_stat(count))))
    ) + 
    scale_y_continuous(
      breaks = 0:4,
      labels = as.character(0:4),
      limits = c(0, 4)
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
    labs(
      x = "Within-host iSNV frequency", 
      y = expression("Cumulative iSNV count ("*log[10]*")"), 
      title = title
    ) + 
    axis_formatting + legend_formatting + background_formatting + 
    theme(
      legend.title = element_blank(), 
      legend.key = element_blank(), 
      legend.position = "none",
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7)
    )
  
  return(Fig1suppY)
}
Fig1suppC <- function_Fig1suppY("17-18_H3N2", "2017–18 H3N2")
Fig1suppD <- function_Fig1suppY("18-19_H3N2", "2018–19 H3N2")
Fig1suppG <- function_Fig1suppY("17-18_H1N1", "2017–18 H1N1")
Fig1suppH <- function_Fig1suppY("18-19_H1N1", "2018–19 H1N1")

#### Figure 2ab - Plot: (gene-wise) within-host PiN and PiS ####
## Make NA not the problematic kind of NA
sg$product[is.na(sg$product)] <- "Neuraminidase"
sg$product[sg$product=="HA_antigenic"] <- "HA\nAntigenic"
sg$product[sg$product=="HA_nonantigenic"] <- "HA\nNonantigenic"
sg$product <- as.character(sg$product)
sg$product[sg$product=="Neuraminidase"] <- "NA"

## Only include some genes
sg <- sg[sg$product %in% factor_by_gene,]

## raw for piN and piS within hosts
yheight <- 15
plot_raw_piN_piS <- function(year_subtype = "", title = "", y_height = yheight) {
  # Filter and clean data
  df <- sg
  if (year_subtype != "") {
    df <- df[df$Merge_ID == year_subtype, ]
  }
  df$product[df$product == "HA_antigenic"] <- "HA\nAntigenic"
  df$product[df$product == "HA_nonantigenic"] <- "HA\nNonantigenic"
  df$product[is.na(df$product)] <- "NA"
  df$product[df$product == "Neuraminidase"] <- "NA"
  df$product <- as.character(df$product)
  df <- df[df$product %in% factor_by_gene, ]
  
  # Reshape
  molten <- df %>%
    dplyr::select(file, product, piN, piS) %>%
    pivot_longer(cols = c(piN, piS), names_to = "Mutation", values_to = "Pi") %>%
    dplyr::mutate(
      Mutation = recode(Mutation, piN = "Nonsynonymous", piS = "Synonymous"),
      product = factor(product, levels = factor_by_gene)
    )
  
  # Normality tests
  normality_results <- molten %>%
    dplyr::group_by(product, Mutation) %>%
    dplyr::summarise(
      n = n(),
      W = shapiro.test(Pi)$statistic,
      p_value = shapiro.test(Pi)$p.value,
      is_normal = p_value > 0.05,
      .groups = "drop"
    )
  
  # Wilcoxon tests
  wide_pi <- molten %>%
    pivot_wider(names_from = Mutation, values_from = Pi) %>%
    dplyr::filter(!is.na(Nonsynonymous) & !is.na(Synonymous))
  
  wilcox_results <- wide_pi %>%
    dplyr::group_by(product) %>%
    dplyr::summarise(
      n = n(),
      statistic = wilcox.test(Nonsynonymous, Synonymous, paired = TRUE)$statistic,
      p_value = wilcox.test(Nonsynonymous, Synonymous, paired = TRUE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  
  # Calculate y-positions for significance labels
  box_heights <- molten %>%
    dplyr::group_by(product, Mutation) %>%
    dplyr::summarise(box_top = quantile(Pi, 0.75, na.rm = TRUE) + 1.5 * IQR(Pi, na.rm = TRUE), .groups = "drop")
  
  max_box_heights <- box_heights %>%
    dplyr::group_by(product) %>%
    dplyr::summarise(y_pos = max(box_top) * 1.05 * 1e4, .groups = "drop")
  
  sig_labels <- wilcox_results %>%
    dplyr::filter(label != "ns") %>%
    dplyr::left_join(max_box_heights, by = "product") %>%
    dplyr::mutate(x = as.numeric(factor(product, levels = factor_by_gene)))
  
  # Plot
  scaling_factor <- 1e4
  exponent <- -log10(scaling_factor)
  
  p <- ggplot(molten, aes(x = product, y = Pi * scaling_factor, fill = Mutation)) +
    geom_quasirandom(
      aes(color = Mutation),
      width       = 0.1,
      dodge.width = 0.8,
      size        = quasi_size-1,
      alpha       = quasi_alpha,
      shape       = 16,
      stroke      = 0
    ) +
    geom_boxplot(outlier.shape = NA, 
                 position = position_dodge(width = 0.8), 
                 width = 0.6,
                 alpha = boxplot_alpha) +
    scale_fill_manual(values = palette_muts_NS_S) +
    scale_color_manual(values = palette_muts_NS_S) +
    scale_y_continuous(limits = c(0, y_height)) +
    geom_segment(data = sig_labels,
                 aes(x = x - 0.2, xend = x + 0.2, y = y_pos, yend = y_pos),
                 inherit.aes = FALSE, linewidth = 0.5, color = "black") +
    geom_text(data = sig_labels,
              aes(x = x, y = y_pos + 0.01 * y_height, label = label),
              inherit.aes = FALSE, size = 4, vjust = 0) +
    labs(x = "", y = bquote("\u03c0 (×10"^.(exponent) * ")"), title = title) +
    theme(
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    axis_formatting +
    legend_formatting +
    background_formatting
  
  print(normality_results)
  print(wilcox_results)
  return(p)
}
plot_piN_piS_wi <- plot_raw_piN_piS("", "", yheight)
plot_piN_piS_wi <- plot_piN_piS_wi + theme(plot.margin = margin(t = 0, r = 4, b = 0, l = 4, "pt"))
plot_piN_piS_wi_1718_H3N2 <- plot_raw_piN_piS("17-18_H3N2", "2017–18 H3N2", yheight)
plot_piN_piS_wi_1718_H1N1 <- plot_raw_piN_piS("17-18_H1N1", "2017–18 H1N1", yheight)
plot_piN_piS_wi_1819_H3N2 <- plot_raw_piN_piS("18-19_H3N2", "2018–19 H3N2", yheight)
plot_piN_piS_wi_1819_H1N1 <- plot_raw_piN_piS("18-19_H1N1", "2018–19 H1N1", yheight)

## raw for piN_minus_piS within hosts
plot_raw_wh_piN_minus_piS <- function(year_subtype = "", title = "", y_min = -6, y_max = 6, seq_by = 2) {
  # Filter data
  df <- sg
  if (year_subtype != "") {
    df <- df[df$Merge_ID == year_subtype, ]
  }
  
  df$product[df$product == "HA_antigenic"] <- "HA\nAntigenic"
  df$product[df$product == "HA_nonantigenic"] <- "HA\nNonantigenic"
  df$product[is.na(df$product)] <- "NA"
  df$product[df$product == "Neuraminidase"] <- "NA"
  df$product <- as.character(df$product)
  df <- df[df$product %in% factor_by_gene, ]
  
  # Rename and clean
  df <- df %>%
    dplyr::select(file, product, piN_minus_piS) %>%
    dplyr::mutate(product = factor(product, levels = factor_by_gene))
  
  # Normality tests
  normality_results <- df %>%
    dplyr::group_by(product) %>%
    dplyr::summarise(
      n = n(),
      W = shapiro.test(piN_minus_piS)$statistic,
      p_value = shapiro.test(piN_minus_piS)$p.value,
      is_normal = p_value > 0.05,
      .groups = "drop"
    )
  
  # Wilcoxon test vs zero
  wilcox_results <- df %>%
    dplyr::group_by(product) %>%
    dplyr::summarise(
      n = n(),
      statistic = wilcox.test(piN_minus_piS, mu = 0)$statistic,
      p_value = wilcox.test(piN_minus_piS, mu = 0)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  # y-position of significance label
  box_heights <- df %>%
    dplyr::group_by(product) %>%
    dplyr::summarise(box_top = quantile(piN_minus_piS, 0.75) + 1.5 * IQR(piN_minus_piS), .groups = "drop")
  sig_labels <- wilcox_results %>%
    dplyr::filter(label != "ns") %>%
    left_join(box_heights, by = "product") %>%
    dplyr::mutate(
      x = as.numeric(factor(product, levels = factor_by_gene)),
      #y_pos = box_top * 1.05 * 1e4
      y_pos = .0005
    )
  
  # Plot
  scaling_factor <- 1e4
  exponent <- -log10(scaling_factor)
  
  p <- ggplot(df, aes(x = product, y = piN_minus_piS * scaling_factor)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
    geom_quasirandom(
      width       = 0.25,
      size        = quasi_size-1,
      alpha       = quasi_alpha,
      shape       = 16,
      stroke      = 0,
      color       = "black"
    ) +
    geom_boxplot(outlier.shape = NA, 
                 width = 0.6, 
                 fill = "grey70", 
                 alpha = boxplot_alpha) +
    geom_text(data = sig_labels,
              aes(x = x, y = y_pos * scaling_factor, label = label),
              inherit.aes = FALSE, size = 4, vjust = 1) +
    scale_y_continuous(limits = c(y_min, y_max),
                       breaks = seq(y_min, y_max, seq_by)) +
    labs(x = "", y = bquote("\u03c0N - \u03c0S (×10"^.(exponent) * ")"), title = title) +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      legend.key = element_blank()
    ) +
    axis_formatting +
    legend_formatting +
    background_formatting
  
  print(normality_results)
  print(wilcox_results)
  return(p)
}
plot_piN_minus_piS_wi <- plot_raw_wh_piN_minus_piS("", "")
plot_piN_minus_piS_wi <- plot_piN_minus_piS_wi + theme(plot.margin = margin(t = 0, r = 4, b = 0, l = 4, "pt"))
plot_piN_minus_piS_wi_1718_H3N2 <- plot_raw_wh_piN_minus_piS("17-18_H3N2", "2017–18 H3N2")
plot_piN_minus_piS_wi_1718_H1N1 <- plot_raw_wh_piN_minus_piS("17-18_H1N1", "2017–18 H1N1")
plot_piN_minus_piS_wi_1819_H3N2 <- plot_raw_wh_piN_minus_piS("18-19_H3N2", "2018–19 H3N2")
plot_piN_minus_piS_wi_1819_H1N1 <- plot_raw_wh_piN_minus_piS("18-19_H1N1", "2018–19 H1N1")

#### Figure 2cd - Plot: (hemagglutinin) sliding-window nucleotide diversity ####
# generate molten data
function_molten_ds_sg_sw <- function(season_subtype) {
  # choose subtype
  sg_sw1 <- sg_sw[sg_sw$Merge_ID==season_subtype,]
  
  # calculate piN-piS
  sg_sw1 <- transform(sg_sw1, piN_minus_piS=piN-piS)
  ds_sg_sw_piN <- ds(sg_sw1[sg_sw1$product=="HA",], varname="piN", groupnames=c("first_site"))
  ds_sg_sw_piS <- ds(sg_sw1[sg_sw1$product=="HA",], varname="piS", groupnames=c("first_site"))
  ds_sg_sw_piNpiS <- ds(sg_sw1[sg_sw1$product=="HA",], varname="piN_minus_piS", groupnames=c("first_site"))
  ds_sg_sw <- merge(ds_sg_sw_piN, ds_sg_sw_piS, by = "first_site")
  ds_sg_sw <- merge(ds_sg_sw,     ds_sg_sw_piNpiS, by = "first_site")
  
  
  # 1) Pivot the means
  molten_mean <- ds_sg_sw %>%
    select(first_site, piN, piS, piN_minus_piS) %>%
    pivot_longer(
      cols      = -first_site,
      names_to  = "Pi",
      values_to = "Mean"
    )
  
  # 2) Pivot the standard errors
  molten_se <- ds_sg_sw %>%
    select(first_site, se.x, se.y, se) %>%
    pivot_longer(
      cols      = -first_site,
      names_to  = "Pi",
      values_to = "SE"
    )
  
  # 3) Pivot the medians
  molten_med <- ds_sg_sw %>%
    select(first_site, median.x, median.y, median) %>%
    pivot_longer(
      cols      = -first_site,
      names_to  = "Pi",
      values_to = "Median"
    )
  
  # 4) Normalize the Pi codes and join all three
  #    se.x → piN, se.y → piS, se → piN_minus_piS
  recode_pi <- c(
    piN                = "piN",
    piS                = "piS",
    piN_minus_piS      = "piN_minus_piS",
    `se.x`             = "piN",
    `se.y`             = "piS",
    `se`               = "piN_minus_piS",
    `median.x`         = "piN",
    `median.y`         = "piS",
    `median`           = "piN_minus_piS"
  )
  
  molten <- molten_mean %>%
    mutate(Pi = recode(Pi, !!!recode_pi)) %>%
    left_join(
      molten_se  %>% mutate(Pi = recode(Pi, !!!recode_pi)),
      by = c("first_site", "Pi")
    ) %>%
    left_join(
      molten_med %>% mutate(Pi = recode(Pi, !!!recode_pi)),
      by = c("first_site", "Pi")
    )
  
  # 5) Turn Pi into nice labels
  molten <- molten %>%
    mutate(
      Pi = factor(
        Pi,
        levels = c("piN", "piS", "piN_minus_piS"),
        labels = c("NS",  "S",   "Diff")
      )
    )
  
  
  return(molten)
  
  
  
  
  ## melt it by pi values
  #molten_ds_sg_sw <- melt(ds_sg_sw,
  #                        id.vars = c("first_site"),
  #                        measure.vars = c("piN", "piS", "piN_minus_piS"),
  #                        variable.name = "Pi")
  #colnames(molten_ds_sg_sw)[2] <- "Pi1"
  #colnames(molten_ds_sg_sw)[3] <- "Mean"
  #molten_ds_sg_sw$Pi1 <- gsub("piN", "Nonsynonymous", molten_ds_sg_sw$Pi1)
  #molten_ds_sg_sw$Pi1 <- gsub("piS", "Synonymous", molten_ds_sg_sw$Pi1)
  #
  ## melt it by pi value SEs
  #molten_ds_sg_sw_se <- melt(ds_sg_sw,
  #                           id.vars = c("first_site"),
  #                           measure.vars = c("se.x", "se.y", "se"),
  #                           variable.name = "SE")
  #colnames(molten_ds_sg_sw_se)[1] <- "first_site2"
  #colnames(molten_ds_sg_sw_se)[2] <- "Pi"
  #colnames(molten_ds_sg_sw_se)[3] <- "SE"
  #molten_ds_sg_sw_se$Pi <- gsub("se.x", "Nonsynonymous", molten_ds_sg_sw_se$Pi)
  #molten_ds_sg_sw_se$Pi <- gsub("se.y", "Synonymous",    molten_ds_sg_sw_se$Pi)
  #molten_ds_sg_sw_se$Pi <- gsub("se",   "piN_minus_piS", molten_ds_sg_sw_se$Pi)
  #
  ### cbind the two molten dfs
  #molten_ds_sg_sw <- cbind(molten_ds_sg_sw, molten_ds_sg_sw_se)
  #molten_ds_sg_sw$Pi1 <- NULL
  #molten_ds_sg_sw$first_site2 <- NULL
  #
  ## gsub NS and S to simplify legend
  #molten_ds_sg_sw$Pi <- gsub("Nonsynonymous", "NS", molten_ds_sg_sw$Pi)
  #molten_ds_sg_sw$Pi <- gsub("Synonymous",     "S", molten_ds_sg_sw$Pi)
  #molten_ds_sg_sw$Pi <- as.factor(molten_ds_sg_sw$Pi)
  #return(molten_ds_sg_sw)
}
molten_ds_sg_sw_1718H3N2 <- function_molten_ds_sg_sw("17-18_H3N2")
molten_ds_sg_sw_1718H1N1 <- function_molten_ds_sg_sw("17-18_H1N1")
molten_ds_sg_sw_1819H3N2 <- function_molten_ds_sg_sw("18-19_H3N2")
molten_ds_sg_sw_1819H1N1 <- function_molten_ds_sg_sw("18-19_H1N1")

# generate piX sliding windows without bootstrapping
function_plot1_piX <- function(molten_ds_sg_sw,
                               season_subtype   = "17-18_H3N2",
                               title            = "2017–18 H3N2",
                               y_height         = NA) {
  
  ## --- constants -----------------------------------------------------------
  scaling_factor <- 1e4                       # multiply π by 10-4
  exponent       <- -log10(scaling_factor)    # used in y-axis label
  
  # colour palette already defined in your workspace
  palette_sw <- palette_muts_NS_S_sw          # NS vs S
  
  ## --- keep only πN & πS (drop πN–πS curve) -------------------------------
  df_plot <- molten_ds_sg_sw %>%
    dplyr::filter(Pi %in% c("NS", "S"))
  
  ## --- ggplot --------------------------------------------------------------
  p <- ggplot(df_plot, aes(x = first_site, colour = Pi, fill = Pi)) +
    geom_ribbon(aes(ymin = (Mean  - SE) * scaling_factor,
                    ymax = (Mean  + SE) * scaling_factor),
                alpha = 0.25, linewidth = 0) +
    geom_line(aes(y = Mean   * scaling_factor), linewidth = 0.8) +
    #geom_line(aes(y = Median * scaling_factor),
    #          linewidth = 0.6, linetype = "dashed") +
    scale_colour_manual(values = palette_sw, name = "Mutation\ntype") +
    scale_fill_manual(values   = palette_sw, name = "Mutation\ntype") +
    scale_y_continuous(limits = c(0, y_height)) +
    labs(x   = "Hemagglutinin nucleotide position",
         y   = bquote("\u03c0 (×10"^.(exponent)*")"),
         title = title) +
    axis_formatting +
    legend_formatting +
    background_formatting +
    theme(legend.key     = element_blank(),
          legend.position = "none")
  
  return(p)
}
yheight <- 30
plot2_sliding_window1_1718H3N2_1 <- function_plot1_piX(molten_ds_sg_sw_1718H3N2, "17-18_H3N2", "2017–18 H3N2", NA)
plot2_sliding_window1_1718H3N2   <- function_plot1_piX(molten_ds_sg_sw_1718H3N2, "17-18_H3N2", "2017–18 H3N2", yheight)
plot2_sliding_window1_1718H1N1   <- function_plot1_piX(molten_ds_sg_sw_1718H1N1, "17-18_H1N1", "2017–18 H1N1", yheight)
plot2_sliding_window1_1819H3N2   <- function_plot1_piX(molten_ds_sg_sw_1819H3N2, "18-19_H3N2", "2018–19 H3N2", yheight)
plot2_sliding_window1_1819H1N1   <- function_plot1_piX(molten_ds_sg_sw_1819H1N1, "18-19_H1N1", "2018–19 H1N1", yheight)

function_plot2_piN_minus_piS <- function(molten_ds_sg_sw,
                                         season_subtype,
                                         title,
                                         y_min = NA) {
  # choose palette & antigenic map based on subtype
  palette_HA_anti_s <- if (grepl("H3N2", season_subtype)) palette_H3_anti_s else palette_H1_anti_s
  sg_sw_antigenic  <- if (grepl("H3N2", season_subtype)) sg_sw_antigenic_H3 else sg_sw_antigenic_H1
  
  # scaling
  scaling_factor <- 1e4
  exponent       <- -log10(scaling_factor)
  
  # get only the Diff (πN–πS) rows
  df <- molten_ds_sg_sw %>%
    dplyr::filter(Pi == "Diff")
  
  # build plot
  p <- ggplot() +
    # 1) antigenic blocks
    geom_rect(
      data = sg_sw_antigenic %>% filter(!is.na(Antigenic_Site)),
      inherit.aes = FALSE,
      aes(
        xmin   = Start,
        xmax   = Stop,
        ymin   = -Inf,
        ymax   = Inf,
        fill   = factor(Antigenic_Site)
      ),
      alpha = 0.25
    ) +
    # 2) horizontal zero line
    geom_hline(yintercept = 0, color = "grey") +
    # 3) SE ribbon
    geom_ribbon(
      data = df,
      aes(
        x    = first_site,
        ymin = (Mean - SE) * scaling_factor,
        ymax = (Mean + SE) * scaling_factor
      ),
      fill     = "black",
      alpha    = 0.2,
      inherit.aes = FALSE
    ) +
    # 4) mean trend line
    geom_line(
      data = df,
      aes(x = first_site, y = Mean * scaling_factor),
      color  = "black",
      size   = 0.8
    ) +
    # 5) median dashed line
    #geom_line(
    #  data = df,
    #  aes(x = first_site, y = Median * scaling_factor),
    #  color     = "black",
    #  linetype  = "dashed",
    #  size      = 0.6
    #) +
    # scales & labels
    scale_fill_manual(values = palette_HA_anti_s, name = "Antigenic\nRegion") +
    scale_y_continuous(limits = c(y_min, 5)) +
    labs(
      x     = "Hemagglutinin nucleotide position",
      y     = bquote("\u03c0N - \u03c0S (×10"^.(exponent)*")"),
      title = title
    ) +
    # styling
    axis_formatting +
    legend_formatting +
    background_formatting +
    theme(
      legend.key      = element_blank(),
      legend.position = "none"
    )
  
  return(p)
}
yheight <- -25
plot2_sliding_window2_1718H3N2_1 <- function_plot2_piN_minus_piS(molten_ds_sg_sw_1718H3N2, "17-18_H3N2", "2017–18 H3N2", NA)
plot2_sliding_window2_1718H3N2   <- function_plot2_piN_minus_piS(molten_ds_sg_sw_1718H3N2, "17-18_H3N2", "2017–18 H3N2", yheight)
plot2_sliding_window2_1718H1N1   <- function_plot2_piN_minus_piS(molten_ds_sg_sw_1718H1N1, "17-18_H1N1", "2017–18 H1N1", yheight)
plot2_sliding_window2_1819H3N2   <- function_plot2_piN_minus_piS(molten_ds_sg_sw_1819H3N2, "18-19_H3N2", "2018–19 H3N2", yheight)
plot2_sliding_window2_1819H1N1   <- function_plot2_piN_minus_piS(molten_ds_sg_sw_1819H1N1, "18-19_H1N1", "2018–19 H1N1", yheight)

#### Figure 3cd - Plot: (gene-wise) piX between hosts ####
# NAs aren't always NAs, sometimes they're a gene!
sg_paired$Gene[is.na(sg_paired$Gene)] <- "NA"

# function to plot the raw data
function_bw_piN_raw <- function(year_subtype, title, y_height, seq_by) {
  # filter to selected year_subtype
  if (year_subtype=="") { # for combined analysis group
    # melt the data by PiN
    molten_sg_paired_piN <- melt(sg_paired,
                                 id.vars = c("ID", "Gene"),
                                 measure.vars = c("Donor_piN", "Recipient_piN"),
                                 variable.name = "PiN")
  }
  if (year_subtype!="") {
    # melt the data by PiN
    molten_sg_paired_piN <- melt(sg_paired[sg_paired$Merge_ID==year_subtype,],
                                 id.vars = c("ID", "Gene"),
                                 measure.vars = c("Donor_piN", "Recipient_piN"),
                                 variable.name = "PiN")
  }
  # clarify the grouping
  molten_sg_paired_piN$Group[molten_sg_paired_piN$PiN=="Donor_piN"] <- "Donor"
  molten_sg_paired_piN$Group[molten_sg_paired_piN$PiN=="Recipient_piN"] <- "Recipient"
  
  # factor by Gene
  molten_sg_paired_piN$Gene <- as.factor(molten_sg_paired_piN$Gene)
  levels <- levels(molten_sg_paired_piN$Gene)
  
  # change factor level order
  molten_sg_paired_piN$Gene <- factor(molten_sg_paired_piN$Gene, levels = factor_by_gene)
  
  # define the scaling factor and negative label
  scaling_factor <- 10^4
  exponent <- -log10(scaling_factor)
  
  # run normality tests
  normality_results <- molten_sg_paired_piN %>%
    dplyr::mutate(
      Gene = as.character(Gene),
      Group = as.character(Group)
    ) %>%
    dplyr::group_by(Gene, Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      W = shapiro.test(value)$statistic,
      p_value = shapiro.test(value)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(is_normal = p_value > 0.05)
  
  # prepare for stats
  piN_wide <- molten_sg_paired_piN %>%
    filter(Group %in% c("Donor", "Recipient")) %>%
    select(ID, Gene, Group, value) %>%
    pivot_wider(names_from = Group, values_from = value) %>%
    drop_na(Donor, Recipient)
  # run the stats
  wilcox_results <- piN_wide %>%
    group_by(Gene) %>%
    dplyr::summarise(
      n = dplyr::n(),
      statistic = wilcox.test(Donor, Recipient, paired = TRUE)$statistic,
      p_value = wilcox.test(Donor, Recipient, paired = TRUE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(significant = p_value < 0.05)
  # simplify the stats
  wilcox_results_simplified <- wilcox_results %>%
    mutate(label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))
  
  # create significance labels
  box_heights <- molten_sg_paired_piN %>%
    dplyr::group_by(Gene, Group) %>%
    dplyr::summarise(box_top = quantile(value, 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE),
                     .groups = "drop")
  max_box_heights <- box_heights %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(y_pos = max(box_top, na.rm = TRUE) * 1.05 * 1e4, .groups = "drop")
  sig_labels <- wilcox_results_simplified %>%
    filter(label != "ns") %>%
    left_join(max_box_heights, by = "Gene") %>%
    mutate(x = match(Gene, levels(molten_sg_paired_piN$Gene)))

  # plot it
  plot_piN_bw_box <- ggplot(molten_sg_paired_piN, 
                            aes(x = Gene, y = value * scaling_factor)) +
    #geom_jitter(aes(color = Group),
    #            position = position_dodge(width = 0.8), size = 0.8, alpha = 0.5) +
    geom_quasirandom(
      aes(color = Group),
      width       = 0.1,
      dodge.width = 0.8,
      size        = quasi_size-0.5,
      alpha       = quasi_alpha,
      shape       = 16,
      stroke      = 0
    ) +
    geom_boxplot(aes(fill = Group),
                 outlier.shape = NA, 
                 position = position_dodge(width = 0.8), 
                 width = 0.6,
                 alpha = boxplot_alpha) +
    scale_fill_manual(values = palette_dr) + 
    scale_color_manual(values = palette_dr) + 
    scale_y_continuous(limits = c(0, y_height),
                       breaks = seq(0, y_height, by = seq_by)) +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.position = "none") + 
    labs(x = "", y = bquote("\u03c0N (×10"^.(exponent) * ")"), 
         title = if (year_subtype == "") NULL else title) + 
    axis_formatting + 
    legend_formatting + 
    background_formatting
  plot_piN_bw_box <- plot_piN_bw_box +
    # Bracket lines
    geom_segment(data = sig_labels,
                 aes(x = x - 0.2,
                     xend = x + 0.2,
                     y = y_pos,
                     yend = y_pos),
                 inherit.aes = FALSE,
                 linewidth = 0.5,
                 color = "black") +
    # Asterisk labels
    geom_text(data = sig_labels,
              aes(x = x, y = y_pos + 0.01 * y_height, label = label),
              inherit.aes = FALSE,
              size = 4,
              vjust = 0)
  # print normality and stats tests and export plot
  print(normality_results)
  print(""); print(wilcox_results)
  return(plot_piN_bw_box)
}
function_bw_piS_raw <- function(year_subtype, title, y_height, seq_by) {
  # filter to selected year_subtype
  if (year_subtype=="") { # for combined analysis group
    # melt the data by PiS
    molten_sg_paired_piS <- melt(sg_paired,
                                 id.vars = c("ID", "Gene"),
                                 measure.vars = c("Donor_piS", "Recipient_piS"),
                                 variable.name = "PiS")
  }
  if (year_subtype!="") {
    # melt the data by PiS
    molten_sg_paired_piS <- melt(sg_paired[sg_paired$Merge_ID==year_subtype,],
                                 id.vars = c("ID", "Gene"),
                                 measure.vars = c("Donor_piS", "Recipient_piS"),
                                 variable.name = "PiS")
  }
  # clarify the grouping
  molten_sg_paired_piS$Group[molten_sg_paired_piS$PiS=="Donor_piS"] <- "Donor"
  molten_sg_paired_piS$Group[molten_sg_paired_piS$PiS=="Recipient_piS"] <- "Recipient"
  
  # factor by Gene
  molten_sg_paired_piS$Gene <- as.factor(molten_sg_paired_piS$Gene)
  levels <- levels(molten_sg_paired_piS$Gene)  
  
  # change factor level order
  molten_sg_paired_piS$Gene <- factor(molten_sg_paired_piS$Gene, levels = factor_by_gene)
  
  # define the scaling factor and negative label
  scaling_factor <- 10^4
  exponent <- -log10(scaling_factor)
  
  # run normality tests
  normality_results <- molten_sg_paired_piS %>%
    dplyr::mutate(
      Gene = as.character(Gene),
      Group = as.character(Group)
    ) %>%
    dplyr::group_by(Gene, Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      W = shapiro.test(value)$statistic,
      p_value = shapiro.test(value)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(is_normal = p_value > 0.05)
  
  # run stats
  piS_wide <- molten_sg_paired_piS %>%
    filter(Group %in% c("Donor", "Recipient")) %>%
    select(ID, Gene, Group, value) %>%
    pivot_wider(names_from = Group, values_from = value) %>%
    drop_na(Donor, Recipient)
  wilcox_results <- piS_wide %>%
    group_by(Gene) %>%
    dplyr::summarise(
      n = dplyr::n(),
      statistic = wilcox.test(Donor, Recipient, paired = TRUE)$statistic,
      p_value = wilcox.test(Donor, Recipient, paired = TRUE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(significant = p_value < 0.05)
  # simplify the stats
  wilcox_results_simplified <- wilcox_results %>%
    mutate(label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))
  
  # create significance labels
  box_heights <- molten_sg_paired_piS %>%
    dplyr::group_by(Gene, Group) %>%
    dplyr::summarise(box_top = quantile(value, 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE),
                     .groups = "drop")
  max_box_heights <- box_heights %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(y_pos = max(box_top, na.rm = TRUE) * 1.05 * 1e4, .groups = "drop")
  sig_labels <- wilcox_results_simplified %>%
    filter(label != "ns") %>%
    left_join(max_box_heights, by = "Gene") %>%
    mutate(x = match(Gene, levels(molten_sg_paired_piS$Gene)))
  
  ## plot it
  plot_piS_bw_box <- ggplot(molten_sg_paired_piS, 
                            aes(x = Gene, y = value * scaling_factor)) +
    geom_quasirandom(
      aes(color = Group),
      width       = 0.1,
      dodge.width = 0.8,
      size        = quasi_size-0.5,
      alpha       = quasi_alpha,
      shape       = 16,
      stroke      = 0
    ) +
    geom_boxplot(aes(fill = Group),
                 outlier.shape = NA, 
                 position = position_dodge(width = 0.8), 
                 width = 0.6,
                 alpha = boxplot_alpha) + 
    scale_fill_manual(values = palette_dr) + 
    scale_color_manual(values = palette_dr) + 
    scale_y_continuous(limits = c(0, y_height),
                       breaks = seq(0, y_height, by = seq_by)) +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.position = "none") + 
    labs(x = "", y = bquote("\u03c0S (×10"^.(exponent) * ")"), 
         title = if (year_subtype == "") NULL else title) + 
    axis_formatting + 
    legend_formatting + 
    background_formatting
  plot_piS_bw_box <- plot_piS_bw_box +
    # Bracket lines
    geom_segment(data = sig_labels,
                 aes(x = x - 0.2,
                     xend = x + 0.2,
                     y = y_pos,
                     yend = y_pos),
                 inherit.aes = FALSE,
                 linewidth = 0.5,
                 color = "black") +
    # Asterisk labels
    geom_text(data = sig_labels,
              aes(x = x, y = y_pos + 0.01 * y_height, label = label),
              inherit.aes = FALSE,
              size = 4,
              vjust = 0)
  
  plot_piS_bw_jitter <- ggplot(molten_sg_paired_piS, aes(x = Gene,
                                                  y = value * scaling_factor,
                                                  fill = Group)) +
    geom_jitter(aes(color = Group),
                position = position_dodge(width = 0.8), size = 0.8, alpha = 0.5) +
    #geom_boxplot(outlier.shape = NA, 
    #             position = position_dodge(width = 0.8), width = 0.6) +
    #scale_fill_manual(values = palette_dr) + 
    scale_color_manual(values = palette_dr) + 
    scale_y_continuous(limits = c(0, y_height)) +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.position = "none") + 
    labs(x = "", y = bquote("\u03c0S (×10"^.(exponent) * ")"), title = title) + 
    axis_formatting + 
    legend_formatting + 
    background_formatting
  #if (year_subtype=="") {plot_piS_bw <- plot_piS_bw + labs(x = "", y = bquote("\u03c0S (×10"^.(exponent) * ")"))}
  #if (year_subtype!="") {plot_piS_bw <- plot_piS_bw + labs(x = "", y = bquote("\u03c0S (×10"^.(exponent) * ")"), title = title)}
  #if (year_subtype=="") {plot_piN_bw <- plot_piN_bw + labs(x = "", y = bquote("\u03c0N (×10"^.(exponent) * ")"))}
  #if (year_subtype!="") {plot_piN_bw <- plot_piN_bw + labs(x = "", y = bquote("\u03c0N (×10"^.(exponent) * ")"), title = title)}
  # make a figure out of the two plots
  plot_piS_bw_raw <- plot_grid(plot_piS_bw_box, plot_piS_bw_jitter, nrow = 1,
                               labels = c("A", "B"), label_size = Size_adjust, 
                               hjust = LR_adjust, vjust = UD_adjust)
  # print normality and stats tests and export plot
  print(normality_results)
  print(""); print(wilcox_results)
  return(plot_piS_bw_box)
}
# piN
plot_piN_bw_raw <- function_bw_piN_raw("", "", 10, 2)
plot_piN_bw_raw_1718_H3N2 <- function_bw_piN_raw("17-18_H3N2", "2017–18 H3N2", 15, 2.5)
plot_piN_bw_raw_1819_H3N2 <- function_bw_piN_raw("18-19_H3N2", "2018–19 H3N2", 15, 2.5)
plot_piN_bw_raw_1718_H1N1 <- function_bw_piN_raw("17-18_H1N1", "2017–18 H1N1", 15, 2.5)
plot_piN_bw_raw_1819_H1N1 <- function_bw_piN_raw("18-19_H1N1", "2018–19 H1N1", 15, 2.5)

# piS
plot_piS_bw_raw <- function_bw_piS_raw("", "", 10, 2)
plot_piS_bw_raw_1718_H3N2 <- function_bw_piS_raw("17-18_H3N2", "2017–18 H3N2", 15, 2.5)
plot_piS_bw_raw_1819_H3N2 <- function_bw_piS_raw("18-19_H3N2", "2018–19 H3N2", 15, 2.5)
plot_piS_bw_raw_1718_H1N1 <- function_bw_piS_raw("17-18_H1N1", "2017–18 H1N1", 15, 2.5)
plot_piS_bw_raw_1819_H1N1 <- function_bw_piS_raw("18-19_H1N1", "2018–19 H1N1", 15, 2.5)

#### Figure 3e  - Plot: (gene-wise) piN-piS between hosts ####
# function to derive piN-piS for a group
function_bw_piN_minus_piS_raw <- function(year_subtype, title, legend, y_height, seq_by, sig_height) {
  # Filter and melt
  if (year_subtype == "") {
    molten <- melt(sg_paired,
                   id.vars = c("ID", "Gene"),
                   measure.vars = c("Donor_piN_minus_piS", "Recipient_piN_minus_piS"),
                   variable.name = "piN_minus_piS")
  } else {
    molten <- melt(sg_paired[sg_paired$Merge_ID == year_subtype, ],
                   id.vars = c("ID", "Gene"),
                   measure.vars = c("Donor_piN_minus_piS", "Recipient_piN_minus_piS"),
                   variable.name = "piN_minus_piS")
  }
  
  molten$Group <- ifelse(molten$piN_minus_piS == "Donor_piN_minus_piS", "Donor", "Recipient")
  molten$Gene <- factor(molten$Gene)
  molten$Gene <- factor(molten$Gene, levels = factor_by_gene)
  
  # Normality check
  normality_results <- molten %>%
    dplyr::group_by(Gene, Group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      W = shapiro.test(value)$statistic,
      p_value = shapiro.test(value)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(is_normal = p_value > 0.05)
  
  # Prepare for Wilcoxon paired test
  paired <- molten %>%
    dplyr::select(ID, Gene, Group, value) %>%
    tidyr::pivot_wider(names_from = Group, values_from = value) %>%
    tidyr::drop_na()
  
  wilcox_results <- paired %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      n = dplyr::n(),
      statistic = wilcox.test(Donor, Recipient, paired = TRUE)$statistic,
      p_value = wilcox.test(Donor, Recipient, paired = TRUE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::mutate(significant = p_value < 0.05)
  
  wilcox_results_simplified <- wilcox_results %>%
    dplyr::mutate(label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ))
  
  # Prepare for Wilcoxon one-sample test
  one_sample_results <- molten %>%
    dplyr::group_by(Gene, Group) %>%
    dplyr::summarise(
      n         = n(),
      statistic = wilcox.test(value, mu = 0)$statistic,
      p_value   = wilcox.test(value, mu = 0)$p.value,
      .groups   = "drop"
    ) %>%
    mutate(
      label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
  print(one_sample_results)
  
  # compute a y‐position just above the top of each box
  box_heights <- molten %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      box_top = quantile(value, 0.75, na.rm = TRUE) + 1.5 * IQR(value, na.rm = TRUE),
      .groups = "drop"
    )
  scaling_factor <- 1e4
  exponent <- -log10(scaling_factor)
  box_heights$box_top <- sig_height/scaling_factor
  
  sig_labels <- one_sample_results %>%
    left_join(box_heights, by = "Group") %>%
    mutate(
      # map Donor→1, Recipient→2
      x     = ifelse(Group=="Donor", 1, 2),
      y_pos = box_top
    )
  
  # Plot
  plot_piN_piS_box <- ggplot(molten, 
                             aes(x = Gene, y = value * scaling_factor)) +
    geom_hline(yintercept = 0, color = "grey") + 
      geom_quasirandom(
        aes(color = Group),
        width       = 0.1,
        dodge.width = 0.8,
        size        = quasi_size-0.5,
        alpha       = quasi_alpha,
        shape       = 16,
        stroke      = 0
      ) + 
    geom_boxplot(aes(fill = Group),
                 outlier.shape = NA, 
                 position = position_dodge(width = 0.8), 
                 width = 0.6,
                 alpha = boxplot_alpha) + 
    #geom_text(data = sig_labels,
    #          aes(x = Gene, y = y_pos * scaling_factor, label = label),
    #          inherit.aes = FALSE,
    #          position = position_dodge(width = 0.8), 
    #          size = 4,
    #          vjust = 0) + 
    geom_text(
      data        = sig_labels,
      aes(
        x     = Gene,
        y     = y_pos * scaling_factor,
        label = label,
        group = Group
      ),
      position = position_dodge(width = 0.8),
      size     = 2,
      vjust    = 1
    ) + 
    scale_fill_manual(values = palette_dr) +
    scale_color_manual(values = palette_dr) +
    scale_y_continuous(limits = c(y_height[1], y_height[2]),
                       breaks = seq(y_height[1], y_height[2], seq_by)) +
    labs(x = "", y = bquote("\u03c0N - \u03c0S (×10"^.(exponent) * ")"),
         title = if (year_subtype == "") NULL else title) +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.position = legend,
          legend.justification = c(0.5, 1),          # anchor at top‐center of legend box
          legend.direction     = "horizontal") +        # entries laid out left→right)
    axis_formatting +
    legend_formatting +
    background_formatting 
  
  print(normality_results)
  print(""); print(wilcox_results)
  return(plot_piN_piS_box)
}
plot_piN_piS_bw <- function_bw_piN_minus_piS_raw("", "", c(0.5, 1), c(-6, 6), 2, 3.5)
y_height <- c(-20, 5)
plot_piN_piS_bw_1718_H3N2 <- function_bw_piN_minus_piS_raw("17-18_H3N2", "2017–18 H3N2", "none", y_height, 5, 5)
plot_piN_piS_bw_1819_H3N2 <- function_bw_piN_minus_piS_raw("18-19_H3N2", "2018–19 H3N2", "none", y_height, 5, 5)
plot_piN_piS_bw_1718_H1N1 <- function_bw_piN_minus_piS_raw("17-18_H1N1", "2017–18 H1N1", "none", y_height, 5, 5)
plot_piN_piS_bw_1819_H1N1 <- function_bw_piN_minus_piS_raw("18-19_H1N1", "2018–19 H1N1", "none", y_height, 5, 5)

#### All plots ####
setwd(dir_s); dir()

### Calculations in manuscript
# % of iSNVs detected below 10% AF
length(df_sub50AF$WSLH_ID[df_sub50AF$AF<10])/length(df_sub50AF$WSLH_ID) # 0.7743506
# median numbers of iSNVs
medians
# prop_var_shared vertical line
p05_threshold

### Fig1
## main
# iSNV frequency spectrum plot
Fig1_iSNV_freq_spec <- plot_grid(plot_mut_bins_prop_1718_H3N2, plot_mut_bins_prop_1819_H3N2, nrow = 1,
                                 labels = c("B", "C"), label_size = Size_adjust, 
                                 hjust = LR_adjust, vjust = UD_adjust)
# manhattan plot with freq spec plot
Fig1_manhattan_spec <- plot_grid(plot_manhattan, Fig1_iSNV_freq_spec,
                                 labels = c("A", ""), 
                                 label_size = Size_adjust, 
                                 rel_heights = c(1.5, 1),
                                 nrow = 2, ncol = 1,
                                 hjust = LR_adjust, 
                                 vjust = UD_adjust)
# save plot
ggsave("F1-manhattan_spec.png", Fig1_manhattan_spec,
       width = 6.5, height = 6, 
       units = "in", device='png', dpi=600) # plot_manhattan is only 17-19 H3N2
## supp
# distribution of iSNV counts, by vax/unvax/either, cumulative iSNV counts
FigS1 <- plot_grid(FigS1a_dist, FigS1b_vax,
                   Fig1suppC, Fig1suppD, Fig1suppG, Fig1suppH,
                   labels = c("A", "B", "C", "D", "E", "F"), ncol = 2,
                   label_size = Size_adjust, hjust = LR_adjust, vjust = UD_adjust, align = "v")
ggsave("S1-vax_counts_pi.png", FigS1,
       width = 6.5, height = 6, 
       units = "in", device='png', dpi=600)

### Fig2
## main
# gene-wise piX and piN-piS 
plot_genewise_wh <- plot_grid(plot_piN_piS_wi, plot_piN_minus_piS_wi,
                              labels = c("A","B"), label_size = Size_adjust, 
                              align = "v",
                              hjust = LR_adjust, vjust = UD_adjust, ncol = 2)
# get legend for sliding window plots
plot_sw_legend1 <- plot2_sliding_window1_1718H3N2_1
plot_sw_legend1 <- plot_sw_legend1 + theme(legend.position="right")
plot_sw_legend1 <- ggpubr::get_legend(plot_sw_legend1)
plot_sw_legend1 <- as_ggplot(plot_sw_legend1)

plot_sw_legend2 <- plot2_sliding_window2_1718H3N2_1
plot_sw_legend2 <- plot_sw_legend2 + theme(legend.position="right")
plot_sw_legend2 <- ggpubr::get_legend(plot_sw_legend2)
plot_sw_legend2 <- as_ggplot(plot_sw_legend2)

plot2_sliding_window12 <- plot_grid(plot2_sliding_window1_1718H3N2_1, 
                                    plot2_sliding_window2_1718H3N2_1, 
                                    align="v", axis ="lr", labels = c("C", "D"),
                                    label_size = Size_adjust, hjust = LR_adjust, vjust = UD_adjust,
                                    nrow = 2)
plot2_sliding_window12_legends <- plot_grid(plot_sw_legend1, plot_sw_legend2, 
                                            align="v", axis ="l", nrow = 2,
                                            vjust = UD_adjust)
plot2_sliding_window12 <- plot_grid(plot2_sliding_window12, 
                                    plot2_sliding_window12_legends, 
                                    align="v", axis ="lr",
                                    hjust = LR_adjust, vjust = UD_adjust,
                                    rel_widths = c(1, .1), nrow = 1)
# plot gene-wise plots with sw plots
plot2 <- plot_grid(plot_genewise_wh, plot2_sliding_window12,
                   axis = "lr",
                   rel_heights = c(1, 2), nrow = 2)
ggsave("F2-genewise_wh_sw_HA.png", plot2,
       width = 6.5, height = 5, 
       units = "in", device='png', dpi=600)

## supp
# piN, piS, piN-piS by season
# get right legend for NS and S
plot_legend_NS_Sr <- plot_piN_piS_wi
plot_legend_NS_Sr <- plot_legend_NS_Sr + theme(legend.position="right", legend.title = element_blank())
plot_legend_NS_Sr <- ggpubr::get_legend(plot_legend_NS_Sr)
plot_legend_NS_Sr <- as_ggplot(plot_legend_NS_Sr)
# supp piX and piN-piS
p1 <- plot_piN_piS_wi_1718_H3N2       + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p2 <- plot_piN_piS_wi_1819_H3N2       + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p3 <- plot_piN_minus_piS_wi_1718_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p4 <- plot_piN_minus_piS_wi_1819_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p5 <- plot_piN_piS_wi_1718_H1N1       + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p6 <- plot_piN_piS_wi_1819_H1N1       + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p7 <- plot_piN_minus_piS_wi_1718_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p8 <- plot_piN_minus_piS_wi_1819_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
plot_genewise_wh_all <- plot_grid(p1, p2, 
                                  p5, p6, 
                                  p3, p4, 
                                  p7, p8,
                                  align="h", axis ="lr",
                                  labels = c("A","B", 
                                             "C", "D", 
                                             "E", "F", 
                                             "G", "H"), label_size = Size_adjust, 
                                  hjust = LR_adjust, vjust = UD_adjust, ncol = 2)
# add legends 
plot_genewise_wh_all_legends <- as_ggplot(ggpubr::get_legend(plot_piN_piS_wi_1718_H3N2 + 
                                                                    theme(legend.position = "bottom",
                                                                          plot.margin = margin(t = 0, r = 4, b = 2, l = 4, "pt"))))
plot_genewise_wh_all_wLegends <- plot_grid(plot_genewise_wh_all, plot_genewise_wh_all_legends,
                                                nrow = 2, rel_heights = c(1, .05))
ggsave("S2-genewise_wh.png", plot_genewise_wh_all_wLegends,
       width = 6.5, height = 6, 
       units = "in", device='png', dpi=600)
# get H1 antigenicity legend for sliding window plots
plot_sw_legend_H1 <- plot2_sliding_window2_1819H1N1
plot_sw_legend_H1 <- plot_sw_legend_H1 + theme(legend.position="right", legend.title = element_blank())
plot_sw_legend_H1 <- ggpubr::get_legend(plot_sw_legend_H1)
plot_sw_legend_H1 <- as_ggplot(plot_sw_legend_H1)
# get H3 antigenicity legend for sliding window plots
plot_sw_legend_H3 <- plot2_sliding_window2_1819H3N2
plot_sw_legend_H3 <- plot_sw_legend_H3 + theme(legend.position="right", legend.title = element_blank())
plot_sw_legend_H3 <- ggpubr::get_legend(plot_sw_legend_H3)
plot_sw_legend_H3 <- as_ggplot(plot_sw_legend_H3)
# supp sw HA piX and piN-piS
plot_sw_wh_all <- plot_grid(plot2_sliding_window1_1718H3N2, plot2_sliding_window1_1819H3N2,
                            plot2_sliding_window2_1718H3N2, plot2_sliding_window2_1819H3N2,
                            plot2_sliding_window1_1718H1N1, plot2_sliding_window1_1819H1N1, 
                            plot2_sliding_window2_1718H1N1, plot2_sliding_window2_1819H1N1,
                            labels = c("A","B", "C", "D", "E", "F", "G", "H"), 
                            align="v", axis ="lr", label_size = Size_adjust, 
                            hjust = LR_adjust, vjust = UD_adjust, ncol = 2)
# add legends
plot_sw_wh_all_legends <- plot_grid(plot_legend_NS_Sr, plot_sw_legend_H3, 
                                    plot_legend_NS_Sr, plot_sw_legend_H1, ncol = 1)
plot_sw_wh_all <- plot_grid(plot_sw_wh_all, plot_sw_wh_all_legends, ncol = 2, rel_widths = c(1, .175))
# sliding window of HA for all season-subtypes
ggsave("S3-sw_wh.png", plot_sw_wh_all,
       width = 6.5, height = 6, 
       units = "in", device='png', dpi=600)

### Fig3
## main
# D-R pairing, JT plot, and π 
# left
plot3_l <- plot_grid(plot_pairing, 
                     plot_grid(plot_piN_bw_raw, plot_piS_bw_raw, nrow = 1, ncol = 2, 
                               labels = c("C", "D"), label_size = Size_adjust, 
                               hjust = LR_adjust, vjust = UD_adjust),
                     labels = c("A", ""), rel_widths = c(1, 1, 1), nrow = 2,
                     label_size = Size_adjust, hjust = LR_adjust, vjust = UD_adjust)
# right
plot3_r <- plot_grid(Fig5b, plot_piN_piS_bw, align = "hv",
                     labels = c("B", "E"), rel_widths = c(1, 1), ncol = 1,
                     label_size = Size_adjust, hjust = LR_adjust, vjust = UD_adjust)
# both
plot3 <- plot_grid(plot3_l, plot3_r, align = "l",
                   rel_widths = c(2, 1), ncol = 2,
                   hjust = LR_adjust, vjust = UD_adjust)
plot3_legend <- as_ggplot(
  ggpubr::get_legend(plot_piN_piS_bw + 
                       theme(legend.position = "bottom",
                             plot.margin = margin(t = 0, r = 4, b = 2, l = 4, "pt"))))
# save
ggsave("F3-prop_shared_JT_pi_hh.png", plot3,
       width = 6.5, height = 3.75, 
       units = "in", device='png', dpi=600)
## supp 5
# H3N2 piN, piS, diff b/w
p1 <- plot_piN_bw_raw_1718_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p2 <- plot_piN_bw_raw_1819_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p3 <- plot_piS_bw_raw_1718_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p4 <- plot_piS_bw_raw_1819_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p5 <- plot_piN_piS_bw_1718_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p6 <- plot_piN_piS_bw_1819_H3N2 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
plotS5 <- plot_grid(p1, p2, p3, p4, p5, p6,
                    labels = c("A","B", "C", "D", "E", "F"), 
                    align="v", axis ="lr", label_size = Size_adjust,
                    hjust = LR_adjust, vjust = UD_adjust, ncol = 2)
plotS5_wLegend <- plot_grid(plotS5, plot3_legend,
                            nrow = 2, rel_heights = c(1, .075))
ggsave("S5-plot_bw_all_H3.png", plotS5_wLegend,
       width = 6.5, height = 5, 
       units = "in", device='png', dpi=600)

## supp 6
# H1N1 piN, piS, diff b/w
p1 <- plot_piN_bw_raw_1718_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p2 <- plot_piN_bw_raw_1819_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p3 <- plot_piS_bw_raw_1718_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p4 <- plot_piS_bw_raw_1819_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p5 <- plot_piN_piS_bw_1718_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
p6 <- plot_piN_piS_bw_1819_H1N1 + theme(plot.margin = margin(t = 4, r = 4, b = 4, l = 4, "pt"), axis.title.x = element_blank())
plotS6 <- plot_grid(p1, p2, p3, p4, p5, p6,
                    labels = c("A","B", "C", "D", "E", "F"), 
                    align="v", axis ="lr", label_size = Size_adjust,
                    hjust = LR_adjust, vjust = UD_adjust, ncol = 2)
plotS6_wLegend <- plot_grid(plotS6, plot3_legend,
                            nrow = 2, rel_heights = c(1, .075))
ggsave("S6-plot_bw_all_H1.png", plotS6_wLegend,
       width = 6.5, height = 5, 
       units = "in", device='png', dpi=600)

### Fig4
## main
# piN-piS for antigenic and non-antigenic regions of HA
plot2_piN_minus_piS_bw_all_legend <- as_ggplot(ggpubr::get_legend(plot_piN_minus_piS_bw_H3N2_1718$plot + 
                                                                    theme(legend.position = "bottom",
                                                                          legend.title = element_blank(),
                                                                          plot.margin = margin(t = 0, r = 4, b = 2, l = 4, "pt"))))
plot2_piN_minus_piS_bw_all_wLegend <- plot_grid(plot2_piN_minus_piS_bw_all, plot2_piN_minus_piS_bw_all_legend,
                                                nrow = 2, rel_heights = c(1, .05))
# piN-piS for all
plot_piN_minus_piS_bw_all_wLegend <- plot_grid(plot_piN_minus_piS_bw$plot, plot2_piN_minus_piS_bw_all_legend,
                                               nrow = 2, rel_heights = c(1, .05))

ggsave("S7-antigenic_diff_all.png", plot2_piN_minus_piS_bw_all_wLegend,
       width = 6.5, height = 5, 
       units = "in", device='png', dpi=600)
# without zeros
ggsave("F4-antigenic_diff_all_viruses.png", plot_piN_minus_piS_bw_all_wLegend,
       width = 6.5, height = 4, 
       units = "in", device='png', dpi=600)
## supp 8
# piN-piS for all regions of HA
plot_piN_minus_piS_bw_by_site_wLegend <- plot_grid(plot_piN_minus_piS_bw_by_site, plot2_piN_minus_piS_bw_all_legend,
                                                   nrow = 2, rel_heights = c(1, .05))
#ggsave("S8-antigenic_diff_all_regions.png", plot_piN_minus_piS_bw_by_site_wLegend,
#       width = 6.5, height = 6.5, 
#       units = "in", device='png', dpi=600)

## supp 9
# timing from donor sampling to recipient sampling
ggsave("S9-donor_recip_timeline.png", plot_dr_sampling,
       width = 6.5, height = 4, 
       units = "in", device='png', dpi=600)


#### ####

