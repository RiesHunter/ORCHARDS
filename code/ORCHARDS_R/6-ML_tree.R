## Clear Global Environment
rm(list = ls())

#### Session prep ####
## Install packages and load libraries as required
# Check if the 'ape' package is installed; if not, install and load it
if(!require(ape)){
  install.packages("ape", dependencies = T)
  library(ape)
}
# Repeat the above for other required packages: 'cowplot', 'phangorn', 'ggplot2', and 'ggtree'
if(!require(cowplot)){
  install.packages("cowplot", dependencies = T)
  library(cowplot)
}
if(!require(phangorn)){
  install.packages("phangorn", dependencies = T)
  library(phangorn)
}
if(!require(ggplot2)){
  install.packages("ggplot2", dependencies = T)
  library(ggplot2)
}
if(!require(ggtree)){
  install.packages("ggtree", dependencies = T)
  library(ggtree)
}
if(!require(RColorBrewer)){
  install.packages("RColorBrewer", dependencies = T)
  library(RColorBrewer)
}
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")  # Install Bioconductor manager if not available
  BiocManager::install("ggtree")  # Install the 'ggtree' package for visualizing phylogenetic trees
  library(ggtree)  # Load the 'ggtree' package
}

#### Import data ####
# Set the directory path for output figures
dir_s <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/figures")
# Set the directory path for input data (specifically for df files)
dir <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/figs/dfs")
# Set the directory for metadata files
dir_md <- paste("/Volumes/PrimateFS/Lab-Friedrich/TCF lab/Current Lab Members/Hunter_Ries/ORCHARDS/data/summary", sep="")
setwd(dir_md)  # Change working directory to the metadata directory
getwd()  # Get the current working directory
dir()  # List files in the current directory

# Import metadata from a CSV file into a data frame
md <- as.data.frame(read_csv("df_samplelist_metadata.csv"))
# Define households with specific sequence identifiers
households <- data.frame("seq" = c("REDACTED"))

# Merge the households data frame with the metadata based on the "seq" and "WSLH_ID" columns
households_md <- merge(households, md, by.x = "seq", by.y = "WSLH_ID")
# Create a new data frame with relevant columns
df_md_hh <- households_md[,c(1,15)]
colnames(df_md_hh) <- c("seq", "cat")


#### Functions ####
# Define a function to root a phylogenetic tree at a specified outgroup
root_tree <- function (x, y) {
  x <- root(x, outgroup = y, 
            resolve.root = TRUE)
  # x <- midpoint.root(x)  # Uncomment to use midpoint rooting
  print(is.rooted(x))
  return(x)
}
# Define a function to calculate likelihood for a given phylogenetic tree and grouping variable
likelihood <- function(x, y) {
  temp <- pml(x, phydat, site.rate = "gamma", model = "GTR")
  temp2 <- data.frame("logLik" = temp$logLik, "group" = y)
  return(temp2)
}
# Define a function to calculate the total length of the tree
treelength <- function(x, y) {
  temp <- sum(node.depth.edgelength(x))
  temp2 <- data.frame("treelength" = temp, "group" = y)
  return(temp2)
}

#### Plot the maximum likelihood trees ####
setwd(dir)  # Change working directory to the data frame directory
getwd()  # Get the current working directory
dir()  # List files in the current directory
# Prepare to read and process multiple trees in a loop
p <- list()  # Initialize an empty list to store plots
for (i in 1:4) {
  # Create a data frame to store information about trees, fasta files, outgroups, and seasons
  table_season_subtype <- data.frame(
    "tree" = c("17-18_H1N1-clustalw.raxml.bestTree", "18-19_H1N1-clustalw.raxml.bestTree", 
               "17-18_H3N2-clustalw.raxml.bestTree", "18-19_H3N2-clustalw.raxml.bestTree"),
    "fasta" = c("17-18_H1N1-clustalw_HA.fasta", "18-19_H1N1-clustalw_HA.fasta", 
                "17-18_H3N2-clustalw_HA.fasta", "18-19_H3N2-clustalw_HA.fasta"),
    "outgroup" = c("H1N1-A_Michigan_45_2015|Reference_sequence", "H1N1-A_Michigan_45_2015|Reference_sequence", 
                   "H3N2-A_Hong_Kong_4801_2014|Reference_sequence", "H3N2-A_Singapore_INFIMH-16-0019_2016|Reference_sequence"),
    "season" = c("2017–18 H1N1", "2018–19 H1N1", "2017–18 H3N2", "2018–19 H3N2"),
    "scale_height" = c(10, 50, 55, 30))
  # Read the best tree from file
  tr_bestTree <- read.tree(paste(table_season_subtype$tree[i]))
  # Read the associated sequence data in FASTA format
  phydat <- read.phyDat(file=table_season_subtype$fasta[i], format="fasta")
  # Root the tree using the specified outgroup
  ClustalW_best <- root_tree(tr_bestTree, table_season_subtype$outgroup[i])
  # Calculate the likelihood for the tree
  df_likelihood <- likelihood(ClustalW_best, "")
  # Calculate the total tree length
  df_treelength <- treelength(ClustalW_best, "")


  
  # Convert the 'cat' column to a factor for proper handling in plotting
  df_md_hh$cat <- as.factor(df_md_hh$cat)
  # Create a named vector to map specific colors to each level for use in the plot
  levels <- levels(df_md_hh$cat)
  
  ## Colors
  # Pull from multiple qualitative palettes
  set1  <- brewer.pal(9, "Set1")
  set2  <- brewer.pal(8, "Set2")
  set3  <- brewer.pal(12, "Set3")
  dark2 <- brewer.pal(8, "Dark2")
  paired <- brewer.pal(12, "Paired")
  
  # Combine and deduplicate (Set3 and Paired may overlap)
  all_colors <- unique(c(set1, set2, set3, dark2, paired))
  
  # Randomly sample X colors
  set.seed(1224)  # for reproducibility
  group_colors <- sample(all_colors, length(levels))
  
  # Add names
  names(group_colors) <- levels(factor(df_md_hh$cat))
  
  # Create a new data frame with the outgroup information from the current iteration
  df_md <- data.frame("seq" = c(table_season_subtype$outgroup[i]),
                      "cat" = c("Outgroup"))
  # Create a complete color palette with a designated color for the outgroup
  group_colors <- c(group_colors, "white")
  names(group_colors) <- c(levels(factor(df_md_hh$cat)), df_md$cat)
  
  # Combine the outgroup data with household metadata, effectively adding back the household info
  df_md <- rbind(df_md, df_md_hh)
  # Convert the 'cat' column to a factor for proper handling in plotting
  df_md$cat <- as.factor(df_md$cat)
  # Create a named vector to map specific colors to each level for use in the plot
  levels <- levels(df_md$cat)
  
  ## add color for uncategorized tips
  # Extract all tip labels from the tree
  all_tips <- ClustalW_best$tip.label
  # Check which are missing from df_md
  tips_missing <- setdiff(all_tips, df_md$seq)
  # Add them as "Uncategorized"
  df_uncat <- data.frame(seq = tips_missing, cat = "Uncategorized")
  # Combine with your original metadata
  df_md <- rbind(df_md, df_uncat)
  # Ensure cat is a factor with correct levels
  df_md$cat <- as.character(df_md$cat)
  df_md$cat[is.na(df_md$cat)] <- "Uncategorized"
  df_md$cat <- factor(df_md$cat)
  # Add color for "Uncategorized" in the palette
  group_colors_full <- c(group_colors, "black")
  names(group_colors_full)[length(group_colors_full)] <- "Uncategorized"
  
  ## Create the phylogenetic tree plot using ggtree, with the current tree and metadata
  p[[i]] <- ggtree(ClustalW_best) %<+% df_md + 
    geom_tippoint(aes(fill = cat), 
                  size = 1,
                  shape = 21,
                  stroke = 0.2,
                  color = "black",
                  position = "identity") +
    theme_tree2() + 
    labs(title = paste(table_season_subtype$season[i])) + 
    theme(legend.position = "none", 
          plot.title = element_text(size = 8, hjust = 0.5)) + 
    scale_fill_manual(values = group_colors_full) 
}

## Combine multiple plots into a grid
Size_adjust = 12         # Adjust size of labels in the plot grid
LR_adjust = -0.5         # Adjust horizontal position of labels
UD_adjust = 1.1          # Adjust vertical position of labels 

## highlight the clade in 18-19 H1N1
tree_1819 <- read.tree("18-19_H1N1-clustalw.raxml.bestTree")
tree_1819 <- root(tree_1819,
                  outgroup = "H1N1-A_Michigan_45_2015|Reference_sequence",
                  resolve.root = TRUE)
#define your clade members again
clade_tips <- c(
  "REDACTED"
)
# mrca node
mrca_node <- getMRCA(tree_1819, clade_tips)
# plot
p2_highlighted <- p[[2]] +
  annotate("rect",
           xmin  = .0105,
           xmax  = .0120,
           ymin  = 05,
           ymax  = 23,
           fill   = NA,
           colour = "red",
           size   = 0.8)

# Combine individual plots into a 2x2 grid with specified labels and positions
plots <- plot_grid(p[[3]], p[[4]], 
                   p[[1]], p2_highlighted, 
                   nrow = 2, ncol = 2,
                   labels = c("A", "B", "C", "D"),  
                   label_size = Size_adjust, hjust = LR_adjust,  vjust = UD_adjust)

# Reset working directory to the original figure directory
setwd(dir_s); dir()

## Save the combined plot as a PDF file
ggsave("S4-ML_trees.png", plots,
       width = 6.5, height = 7, 
       units = "in", device="png", dpi = 600)

#### ####






library(ape)       # trees & distances
library(phangorn)  # reading the alignment (phydat)
library(ggplot2)   # plotting
library(dplyr)     # summary stats


tree_file  <- "18-19_H1N1-clustalw.raxml.bestTree"
fasta_file <- "18-19_H1N1-clustalw_HA.fasta"
ref_tip    <- "H1N1-A_Michigan_45_2015|Reference_sequence"


tr  <- read.tree(tree_file)
tr  <- root(tr, outgroup = ref_tip, resolve.root = TRUE)

ha <- read.phyDat(fasta_file, format = "fasta")


mat        <- cophenetic(tr)          # full pairwise patristic matrix
dist_ref   <- mat[, ref_tip]          # column for the reference
dist_ref   <- dist_ref[names(dist_ref) != ref_tip]  # drop self‑distance

df_dist <- data.frame(seq = names(dist_ref),
                      dist = dist_ref)  # substitutions per site

p <- ggplot(df_dist, aes(x = dist)) +
  geom_histogram(binwidth = 0.0005, colour = "black", fill = "#3182bd") +
  theme_minimal(base_size = 11) +
  labs(title = "2018–19 H1N1: HA divergence from Michigan/45/2015",
       x = "Patristic distance (substitutions per site)",
       y = "Number of genomes")

print(p)



stats <- df_dist %>%
  dplyr::summarise(n      = n(),
            mean   = mean(dist),
            median = median(dist),
            sd     = sd(dist),
            min    = min(dist),
            max    = max(dist),
            p95    = quantile(dist, 0.95))

print(stats)

ha_len <- sum(attr(ha, "weight"))
stats_nt <- mutate(stats,
                   mean_nt   = mean   * ha_len,
                   median_nt = median * ha_len,
                   max_nt    = max    * ha_len)
print(stats_nt)





tr <- read.tree("18-19_H1N1-clustalw.raxml.bestTree")
tr <- root(tr, outgroup="H1N1-A_Michigan_45_2015|Reference_sequence",
           resolve.root=TRUE)
node  <- getMRCA(tr, c("REDACTED"))
clade <- extract.clade(tr, node)
tips  <- clade$tip.label

df_dist <- data.frame(seq = names(dist_ref), dist = dist_ref)

clade_tips <- c(
  "REDACTED"
)

df_dist <- df_dist %>%
  mutate(
    in_clade = ifelse(seq %in% clade_tips, "Clade member", "Other sample")
  )




p_dist <- ggplot(df_dist, aes(x = seq, y = dist, color = in_clade)) +
  geom_point(size = 2) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  ) +
  scale_color_manual(
    name   = NULL,
    values = c("Clade member" = "#de2d26", "Other sample" = "grey70")
  ) +
  labs(
    x     = "Sample",
    y     = "Patristic distance (subs/site)",
    title = "2018–19 H1N1 HA: divergence from Michigan/45"
  )








# 2. files and outgroup
tree_file <- "18-19_H1N1-clustalw.raxml.bestTree"
ref_tip    <- "H1N1-A_Michigan_45_2015|Reference_sequence"

# 3. read & root
tr <- read.tree(tree_file)
tr <- root(tr, outgroup = ref_tip, resolve.root = TRUE)

# 4. define your focal clade
clade_tips <- c(
  "REDACTED"
)

# 5. make a tip‑metadata table
df_meta <- data.frame(seq = tr$tip.label) %>%
  mutate(
    in_clade = ifelse(seq %in% clade_tips, "Clade member", "Other sample")
  )

# 6. choose colors
clade_colors <- c(
  "Clade member" = "#de2d26",
  "Other sample" = "grey70"
)

# 7. plot
p_tree <- ggtree(tr) %<+% df_meta +
  geom_tippoint(aes(color = in_clade), size = 2) +
  scale_color_manual(values = clade_colors) +
  theme_tree2() +
  theme(
    legend.position = c(0.825, 0.90),
    legend.text = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5)
  ) +
  labs(
    title = "2018–19 H1N1 HA phylogeny",
    color = NULL
  )



plot_grid(p_tree, p_dist)



