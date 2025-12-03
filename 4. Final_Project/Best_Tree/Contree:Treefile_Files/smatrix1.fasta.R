############################################################
# 0. Load packages
############################################################
library(ape)
library(ggtree)
library(dplyr)
library(ggplot2)

############################################################
# 1. Read ML tree
############################################################
ml_tree <- read.tree("smatrix1.fasta.treefile")
ml_tree$node.label <- as.character(ml_tree$node.label)  

############################################################
# 2. Reroot tree
############################################################
outgroup_name <- "Cynocephalus_volans"
if(!(outgroup_name %in% ml_tree$tip.label)) stop("Outgroup not found!")
ml_tree <- root(ml_tree, outgroup = outgroup_name, resolve.root = TRUE)

############################################################
# 3. Read traits
############################################################
traits <- read.table("50_taxa_list.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(traits) <- c("species", "lucidum")
traits$lucidum[grep("^Homo_sapiens", traits$species)] <- "N"

df <- data.frame(species = ml_tree$tip.label) %>%
  left_join(traits, by="species")
df$lucidum[is.na(df$lucidum)] <- "?"

############################################################
# 4. Build ggtree object
############################################################
p <- ggtree(ml_tree, layout = "rectangular") %<+% df +
  geom_tiplab(aes(color = lucidum), size = 7) +
  scale_color_manual(
    values = c("L" = "green", "N" = "red", "?" = "gray"),
    labels = c("L" = "Lucidum", "N" = "No Lucidum", "?" = "Unknown"),
    name = "Lucidum status"
  ) +
  guides(color = guide_legend(
    override.aes = list(shape = 15, size = 10),  # bigger legend boxes
    keywidth = 2,      # width of legend key
    keyheight = 2,     # height of legend key
    title.theme = element_text(size = 22),  # larger legend title
    label.theme = element_text(size = 20)   # larger legend text
  )) +
  theme_tree2() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    axis.text.x = element_text(size = 18),     
    axis.title.x = element_text(size = 20)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.6)))

############################################################
# 5. Add bootstrap support values
############################################################
# Extract ggtree internal node coordinates
tree_data <- p$data
internal_nodes <- tree_data[!tree_data$isTip, ]

# Add bootstrap labels (from ml_tree$node.label)
internal_nodes$label <- ml_tree$node.label

# Add labels to plot
p <- p + geom_text2(
  data = internal_nodes,
  aes(x = x, y = y, label = label),
  hjust = -0.2,
  vjust = 0.5,
  size = 6,
  color = "black"
)

############################################################
# 6. Display plot
############################################################
p

