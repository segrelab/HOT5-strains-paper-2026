library(ape)
setwd("~/npsegre/melisa2026")

# Read the tree
# tree <- read.tree("new-tree.newick")
tip.label <- read.tree("new-tree-labels.newick")

# Inspect it
print(tip.label)
summary(tip.label)
tip.label$tip.label    # see your taxon names

# Basic plot
par(family = "sans") #set font
plot(tip.label, main = "Phylogenetic Tree")
tip.label$tip.label[42] <- "HOT5 Strain E06"
tip.label$tip.label[41] <- "HOT5 Strain B10"
tip.label$tip.label[40] <- "HOT5 Strain F03"
tip.label$tip.label[39] <- "Marinovum algicola DG898"
tip.label$tip.label[38] <- "Roseobacter litoralis Och149"
tip.label$tip.label[36] <- "Rosoebacter denitrificans OCh114"
tip.label$tip.label[35] <- "Leisingera aquaemixtae"
tip.label$tip.label[34] <- "Leisingera adaeponensis DSM23529"
tip.label$tip.label[33] <- "Leisingera methylhalidivorans DSM23566"
tip.label$tip.label[32] <- "Phaoebacter sp. 11ANDIMAR09"
tip.label$tip.label[31] <- "PseudophaeobacterarcticusDSM23566"
tip.label$tip.label[30] <- "Roseobactersp. SK209-2-6"   
tip.label$tip.label[29] <- "Epibacterium multivorans"    
tip.label$tip.label[28] <- "Ruegeria mobilis F1926"   
tip.label$tip.label[27] <- "Epibacterium ulvae"     
tip.label$tip.label[26] <- "Phaeobacter inhibens DSM16374"
tip.label$tip.label[25] <- "Phaeobacter inhibens DSM17395"     
tip.label$tip.label[24] <- "Phaeobacter inhibens 2.10"
tip.label$tip.label[23] <- "Phaeobacter gallaeciensis DSM26640" 
tip.label$tip.label[22] <- "Phaeobacteritalicus" 
tip.label$tip.label[21] <- "Ruegeria marina"   
tip.label$tip.label[20] <- "Ruegeriapomeroyi DSS-3"
tip.label$tip.label[18] <- "Pelagicola litoralis"    
tip.label$tip.label[17] <- "Phaeobacter sp. CECT7735"   
tip.label$tip.label[16] <- "Roseovarius mucosus DSM17069" 
tip.label$tip.label[15] <- "Roseovarius sp.217"   
tip.label$tip.label[14] <- "HOT5 Strain B08"
tip.label$tip.label[13] <- "HOT5 Strain C03"
tip.label$tip.label[12] <- "Roseovarius atlanticus"
tip.label$tip.label[10] <- "Oceanicola sp. S124"
tip.label$tip.label[8] <- "Rhodobacter sp. LPB0142"
tip.label$tip.label[7] <- "Thalassobaculums alexigens DSM19539"
tip.label$tip.label[5] <- "Alteromonas macleodii"   
tip.label$tip.label[4] <- "Alteromonas macleodii ATCC27126"  
tip.label$tip.label[3] <- "Alteromonas australica" 
tip.label$tip.label[2] <- "Alteromonas sp. LOR" 
tip.label$tip.label[1] <- "Alteromonas stellipolaris"

# Remove duplicate tips (identical genomes, ANI distance = 0.0)
tip.label <- drop.tip(tip.label, c(6, 9, 11, 19, 37))

tip_fonts <- rep(1, length(tip.label$tip.label))
tip_fonts[tip.label$tip.label %in% c("HOT5 Strain B08", "HOT5 Strain C03",
                                      "HOT5 Strain B10", "HOT5 Strain E06",
                                      "HOT5 Strain F03")] <- 2

plot(tip.label, font = tip_fonts, cex = 0.6)
add.scale.bar()

