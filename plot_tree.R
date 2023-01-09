


library(phytools)
library(readxl)
library(tidyverse)
tree <- read.tree("/Users/raphael/Library/CloudStorage/Box-Box/KenyaAtlasComparison/data/addTaxaComplete2Sep2022.tre")
lost_kept_gain_sp <- read_xlsx("/Users/raphael/Library/CloudStorage/Box-Box/KenyaAtlasComparison/export/lost_kept_gain_sp.xlsx")

# Map old atlas to eBird taxa
sp_ebird = readxl::read_xlsx('/Users/raphael/Library/CloudStorage/Box-Box/KenyaAtlasComparison/data/eBird/sp_ebird.xlsx')%>%
  filter(SEQ!=0) %>%
  mutate(
    scientifique_name = str_replace_all(scientifique_name," ","_")
  ) %>%
  filter(scientifique_name %in% tree$tip.label) %>%
  left_join(lost_kept_gain_sp) %>%
  group_by(scientifique_name) %>%
  summarise(diff_score = sum(diff_score))

pruned <- drop.tip(tree, setdiff(tree$tip.label, sp_ebird$scientifique_name))
x <- sp_ebird$diff_score
names(x) <- sp_ebird$scientifique_name

contMap(pruned, x, type="fan",fsize=0.9)

plotTree(pruned, x)

write.tree(pruned, "/Users/raphael/Library/CloudStorage/Box-Box/KenyaAtlasComparison/data/kenyan_bird.tre")
