library(ape)
setwd("~/npsegre/melisa2026")

# Read the tree
# tree <- read.tree("new-tree.newick")
tip.label3 <- read.tree("Rhodobacters-SCGs-clean.fa.contree")

# Inspect it
print(tip.label3)
summary(tip.label3)
tree$tip.label3    # see your taxon names

# Basic plot
par(family = "sans") #set font
plot(tip.label3, main = "Phylogenetic Tree", cex=0.1)
add.scale.bar()
