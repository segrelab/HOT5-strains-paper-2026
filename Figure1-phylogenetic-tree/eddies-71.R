library(ape)
setwd("~/npsegre/melisa2026")

#write the tree
# Read the newick string from txt
newick_string <- readLines("Bac_71_phylogenomic-tree.txt")

# Parse and plot
tree <- read.tree(text = newick_string)
plot(tree)
write.tree(tree, file = "tip.label2.nwk")

# Read the tree
# tree <- read.tree("new-tree.newick")
tip.label2 <- read.tree("tip.label2.nwk")

# Inspect it
print(tip.label2)
summary(tip.label2)
tree$tip.label2    # see your taxon names

# Basic plot
par(family = "sans") #set font
plot(tip.label2, main = "Phylogenetic Tree", cex=0.1)
add.scale.bar()
