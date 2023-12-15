############################# 
#### Trees for Figure 1B ####
#############################

library(phytools)
library(ggtree)
parathaaDir <- ("C:/Users/mshort/Documents/proj/parathaa/")

#V1V2
load(file.path(parathaaDir, "output/20231204_Kingdom_V1V2/resultTree_bestThresholds.RData"))
in.jplace <- read.jplace(file.path(parathaaDir, "output/20231204_Kingdom_V1V2/merged_sub.jplace"))

in.tree.data <- resultData$tax_bestcuts
##Add unique labels for interior nodes:
in.tree.data <- in.tree.data %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
  select(!label) %>%
  rename(label=label_new)

load("output/Figures/synth_nomult_arc/V1V2_full_comparisons.RData")
nm <- "JOUE01000003.131499.133012"
pind <- 2
plotTree <- as.phylo(in.tree.data)
plc <- in.jplace@placements[which(in.jplace@placements$name==nm),]
plc <- plc %>% filter(like_weight_ratio > 0.5*max(like_weight_ratio))
plc

print(pind)
plotTree <- bind.tip(plotTree, tip.label=paste0("Placement_", pind), edge.length=plc$pendant_length[pind], where=plc$node[pind], position=plc$distal_length[pind])
plc1 <- which(plotTree$tip.label==paste0("Placement_", pind))
print(plc1)
test.tree <- tree_subset(plotTree, node=plc1, levels_back = 4)
if(nrow(as_tibble(test.tree))>50)
  test.tree <- tree_subset(plotTree, node=plc1, levels_back = 2)

labelInfo <- in.tree.data %>% select(label, Kingdom, Phylum, Class, Order, Family, Genus, Species)
plotTree2 <- left_join(as_tibble(test.tree), labelInfo)
truName <- compare.synth %>% filter(AccID==nm) %>% select(Species.silva) %>% as.character()
parName <- compare.synth %>% filter(AccID==nm) %>% select(Species.parathaa) %>% as.character()

plotData3 <- as.treedata(plotTree2) %>%
  mutate(Label1 = if_else(label==paste0("Placement_", pind), " Query", NA)) %>%
  mutate(Species = gsub(";", ";\n", Species)) 
plotTree <- as.phylo(in.tree.data)


t1 <- ggtree(plotData3, aes(color=Species)) + geom_tippoint(size=5) + geom_nodepoint(size=5) + 
  geom_tiplab(aes(label =Label1), size=6, show.legend=F) + 
  labs(title=paste0("Source: ", truName),  subtitle = paste0("Parathaa: ", parName)) +
  geom_treescale(fontsize = 6, linesize = 1, width = 0.005, x = .1) + theme(text=element_text(size=16)) +
  theme(legend.position = "bottom") + guides(color=guide_legend(nrow=3,byrow=TRUE)) + xlim(c(0,0.15)) + scale_color_manual(values=c("#F8766D", "#C49A00" ,"#53B400" ))
t1
ggsave(t1, filename = "output/Figures/Example2_v1v2.pdf", width = 8, height = 6, units="in")


#V4V5
load(file.path(parathaaDir, "output/20231203_kindom_fix_synth/resultTree_bestThresholds.RData"))
in.jplace <- read.jplace(file.path(parathaaDir, "output/20231203_kindom_fix_synth/merged_sub.jplace"))

in.tree.data <- resultData$tax_bestcuts
##Add unique labels for interior nodes:
in.tree.data <- in.tree.data %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
  select(!label) %>%
  rename(label=label_new)

load("output/Figures/synth_nomult_arc/V4V5_full_comparisons.RData")
nm <- "JOUE01000003.131499.133012"
pind <- 1
plotTree <- as.phylo(in.tree.data)
plc <- in.jplace@placements[which(in.jplace@placements$name==nm),]
plc <- plc %>% filter(like_weight_ratio > 0.5*max(like_weight_ratio))
plc

print(pind)
plotTree <- bind.tip(plotTree, tip.label=paste0("Placement_", pind), edge.length=plc$pendant_length[pind], where=plc$node[pind], position=plc$distal_length[pind])
plc1 <- which(plotTree$tip.label==paste0("Placement_", pind))
print(plc1)
test.tree <- tree_subset(plotTree, node=plc1, levels_back = 4)
if(nrow(as_tibble(test.tree))>50)
  test.tree <- tree_subset(plotTree, node=plc1, levels_back = 3)

labelInfo <- in.tree.data %>% select(label, Kingdom, Phylum, Class, Order, Family, Genus, Species)
plotTree2 <- left_join(as_tibble(test.tree), labelInfo)
truName <- compare.synth %>% filter(AccID==nm) %>% select(Species.silva) %>% as.character()
parName <- compare.synth %>% filter(AccID==nm) %>% select(Species.parathaa) %>% as.character()
plotData3 <- as.treedata(plotTree2) %>%
  mutate(Label1 = if_else(label==paste0("Placement_", pind), " Query", NA)) %>%
  mutate(Species = gsub(";", ";\n", Species)) 




t1 <- ggtree(plotData3, aes(color=Species)) + geom_tippoint(size=5) + geom_nodepoint(size=5) + 
  geom_tiplab(aes(label =Label1), size=6, show.legend=F) + 
  labs(title=paste0("Source: ", truName),  subtitle = paste0("Parathaa: ", parName)) +
  geom_treescale(fontsize = 6, linesize = 1, width = 0.005, x = .1) + theme(text=element_text(size=16)) +
  theme(legend.position = "bottom") + guides(color=guide_legend(nrow=3,byrow=TRUE)) + xlim(c(0,0.15)) + scale_color_manual(values =  c("Francisella noatunensis"="#F8766D", 
                                                                                                                                       "Francisella philomiragia" = "#C49A00", 
                                                                                                                                       "Piscirickettsia salmonis" = "#53B400",
                                                                                                                                       "Francisella noatunensis;\nFrancisella philomiragia" = "#00C094",
                                                                                                                                       "Legionella longbeachae"= "#00B6EB",
                                                                                                                                       "Legionella lytica" = "#A58AFF")) 
t1
ggsave(t1, filename = "output/Figures/Example2_v4v5.pdf", width = 8, height = 6, units="in")
