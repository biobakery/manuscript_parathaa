############################# 
#### Trees for Figure 1B ####
#############################


require(docopt)

'Usage:
  Plots.Figure.1b.R [-p <parathaa_PATH> --namedTreeV1V2 <named_tree_V1V2> --jV1V2 <jplace_V1V2> -o <output> --bV1V2 <benchmark_V1V2> --namedTreeV4V5 <named_tree_V4V5> --jV4V5 <jplace_V4V5> --bV4V5 <benchmark_V4V5>]
  
  
Options:
  -p directory where parathaa github repo is cloned
  --namedTreeV1V2 the named tree RData file for V1V2
  --jV1V2 the jplace file from V1V2 synthetic data testing
  --bV1V2 benchmarking V1V2 results
  --namedTreeV4V5 the named tree RData file for V1V2
  --jV4V5 the jplace file from V1V2 synthetic data testing
  --bV4V5 benchmarking V1V2 results
  -o output
  ]' -> doc


opts <- docopt(doc)

library(phytools)
library(ggtree)
parathaaDir <- opts$p

#V1V2
load(opts$namedTreeV1V2)
in.jplace <- read.jplace(opts$jV1V2)

in.tree.data <- resultData$tax_bestcuts
##Add unique labels for interior nodes:
in.tree.data <- in.tree.data %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
  select(!label) %>%
  rename(label=label_new)

load(opts$bV1V2)
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
ggsave(t1, filename = paste0(opts$o, "/Figures/Figure1B_V1V2.pdf"), width = 8, height = 6, units="in")


#V4V5
load(opts$namedTreeV4V5)
in.jplace <- read.jplace(opts$jV4V5)

in.tree.data <- resultData$tax_bestcuts
##Add unique labels for interior nodes:
in.tree.data <- in.tree.data %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
  select(!label) %>%
  rename(label=label_new)

load(opts$bV4V5)
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
ggsave(t1, filename = paste0(opts$o, "/Figures/Figure1B_V4V5.pdf"), width = 8, height = 6, units="in")
