library(phytools)
library(readxl)
library(tidyverse)
library(RColorBrewer)

fld <- "/Users/raphael/Library/CloudStorage/Box-Box/KenyaAtlasComparison/"
tree <- read.tree(paste0(fld,"data/phylo/addTaxaComplete2Sep2022.tre"))
sp_lost_kept_gain <- read_csv(paste0(fld,"export/sp_lost_kept_gain.csv"))
sp_ebird = read_xlsx(paste0(fld,"data/eBird/sp_ebird.xlsx")) %>% 
  transmute(
    tre_label = str_replace_all(scientific_name," ","_"), 
    SEQ
    ) %>% 
  filter(SEQ %in% sp_lost_kept_gain$SEQ) # remove species which have been observed since the old atlas and SEQ=0
 
# Prune the tree to keep only species from sp_old/sp_lost_kept
pruned <- drop.tip(tree, setdiff(tree$tip.label, sp_ebird$tre_label))

# Find the equivalent seq of the tip of the tree
pruned_label_seq = sp_ebird$SEQ[match(pruned$tip.label, sp_ebird$tre_label)]

# Merge tree tip which have the same SEQ
while (any(duplicated(pruned_label_seq))){
  
  # Select the first duplicate
  seq_dup <- head(pruned_label_seq[duplicated(pruned_label_seq)],1)
  # Find the species name which need to be merged
  sp_dup <- pruned$tip.label[pruned_label_seq==seq_dup]
  print(sp_dup)
  
  # Find the node of the closer parent
  c <- getMRCA(pruned, sp_dup)
  stopifnot(!is.null(c))
  
  # Add a tip at the parent node
  pruned = bind.tip(pruned, "temp", edge.length=NULL, where=c, position=0)
  # delete the two children nodes
  pruned = drop.tip(pruned, sp_dup)
  pruned$tip.label[pruned$tip.label=="temp"] <- sp_dup[1]
  
  # Update the seq with the new tree
  pruned_label_seq = sp_ebird$SEQ[match(pruned$tip.label, sp_ebird$tre_label)]
}

# Update the tips label
pruned$tip.label <- pruned_label_seq

# Add random name to node label to be able to read in matlab
n <- length(pruned$node.label)
a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
pruned$node.label <- paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))

# Export tree
write.tree(pruned, paste0(fld,"data/phylo/kenyan_bird.tre"))




## Pagel's lambda
pruned2 <- drop.tip(pruned, setdiff(pruned$tip.label, sp_lost_kept_gain %>% .$SEQ))
id = match(pruned2$tip.label, sp_lost_kept_gain$SEQ)
x <- sp_lost_kept_gain$gain[id] - sp_lost_kept_gain$lost[id]
lambda<-phylosig(pruned2, x, method="lambda", test=T)



# Plot
# 
pruned2 <- drop.tip(pruned, setdiff(pruned$tip.label, sp_lost_kept_gain %>% filter((old+new)>10) %>% .$SEQ))

# Use sp_old scientific name (rather than ebird)
id = match(pruned2$tip.label, sp_lost_kept_gain$SEQ)
x <- sp_lost_kept_gain$gain[id]-sp_lost_kept_gain$lost[id]
pruned2$tip.label <- paste0(sp_lost_kept_gain$CommonName[id], ": ", ifelse(sign(x)==1,"+",""), x)
pruned2$tip.label <- paste0(sp_lost_kept_gain$Family[id],' - ', sp_lost_kept_gain$CommonName[id])
names(x) <- pruned2$tip.label
m <- contMap(pruned2, x, plot=F, fsize=0.1, lims=c(-20, 20))
m <-setMap(m, brewer.pal(11,'RdYlGn'))


setEPS()
postscript("whatever.eps")
plot(m, type="fan", fsize=0.1, lwd=1, outline=F)
dev.off()

plotTree(pruned, x)
