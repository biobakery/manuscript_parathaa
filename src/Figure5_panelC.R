
library(phytools)
library(treeio)
library(ggtree)
library(dplyr)
library(RColorBrewer)
library(microbiome)

## ORAL - V4V5
parathaaDir <- "output/Oral_V4V5/"
in.tree <- read.tree(file.path(parathaaDir, "ref.refpkg/region_specific.tree"))
in.jplace <- read.jplace(file.path(parathaaDir, "merged_sub.jplace"))
merged_phyloseq <- readRDS(file.path(parathaaDir, "merged_phyloseq.RDS"))




tax.ids <- data.frame(tax_table(merged_phyloseq)) %>% filter(grepl(";", Genus)) %>% rownames()

#Most abundant genera with ambig assignments
merged_physeq_ge <- aggregate_taxa(merged_phyloseq, "Genus", verbose = FALSE)
tax.ids <- data.frame(tax_table(merged_physeq_ge)) %>% filter(grepl(";", Genus)) %>% rownames()
otu.amb <- otu_table(merged_physeq_ge)[tax.ids,]
sort(rowSums(otu.amb), decreasing = T)[1:10]

## Check out a few OTUs to see how they are assigned by both Parathaa and DADA2
rows1 <- data.frame(tax_table(merged_phyloseq)) %>% filter(Genus=="Enterobacter;Erwinia;Pantoea;Raoultella") %>% rownames()
rows1.1 <- substr(rows1,1,nchar(rows1)-1)
rowSums(otu_table(merged_phyloseq)[rows1,])
data.frame(tax_table(merged_phyloseq)) %>% filter(grepl(rows1.1[2], rownames(data.frame(tax_table(merged_phyloseq))))) 

rows1 <- data.frame(tax_table(merged_phyloseq)) %>% filter(Genus=="Haemophilus;Pasteurella") %>% rownames()
rows1.1 <- substr(rows1,1,nchar(rows1)-1)
rowSums(otu_table(merged_phyloseq)[rows1,])
data.frame(tax_table(merged_phyloseq)) %>% filter(grepl(rows1.1[1], rownames(data.frame(tax_table(merged_phyloseq))))) 

rows1 <- data.frame(tax_table(merged_phyloseq)) %>% filter(Genus=="Haemophilus;Mannheimia;Pasteurella") %>% rownames()
rows1.1 <- substr(rows1,1,nchar(rows1)-1)
rowSums(otu_table(merged_phyloseq)[rows1,])
data.frame(tax_table(merged_phyloseq)) %>% filter(grepl(rows1.1[12], rownames(data.frame(tax_table(merged_phyloseq))))) 

rows1 <- data.frame(tax_table(merged_phyloseq)) %>% filter(Genus=="Necropsobacter;Pasteurella") %>% rownames()
rows1.1 <- substr(rows1,1,nchar(rows1)-1)
rowSums(otu_table(merged_phyloseq)[rows1,])
data.frame(tax_table(merged_phyloseq)) %>% filter(grepl(rows1.1[5], rownames(data.frame(tax_table(merged_phyloseq))))) 

## Pasteurellaceae
tax.ids <- data.frame(tax_table(merged_physeq_ge)) %>% filter(grepl("Pasteurellaceae", Family)) %>% rownames()
otu.pas <- otu_table(merged_physeq_ge)[tax.ids,]
sort(rowSums(otu.pas), decreasing = T)

dada2_phy <- subset_samples(merged_phyloseq, Taxonomy_type=="DADA2")
tax.ids <- data.frame(tax_table(dada2_phy)) %>% filter(grepl("Pasteurellaceae", Family)) %>% rownames()
otu.pas <- otu_table(dada2_phy)[tax.ids,]
sort(rowSums(otu.pas), decreasing = T)[1:10]

print("Relative abundance in DADA2:")
sum(rowSums(otu.pas))/sum(rowSums(otu_table(dada2_phy)))


parathaa_phy <- subset_samples(merged_phyloseq, Taxonomy_type=="Parathaa")
tax.ids_phy <- data.frame(tax_table(parathaa_phy)) %>% filter(grepl("Pasteurellaceae", Family)) %>% rownames()
otu.pas_phy <- otu_table(parathaa_phy)[tax.ids_phy,]
sort(rowSums(otu.pas), decreasing = T)[1:10]

print("Relative abundance in Parathaa:")
sum(rowSums(otu.pas_phy))/sum(rowSums(otu_table(parathaa_phy)))


## Calculate genus relative abundances in full data
merged_physeq_ge_rel <- transform_sample_counts(merged_physeq_ge, function(x) x/sum(x))
## Subset full relabund dataset to Pasteurellaceae
merged_physeq_ge_pas <- subset_taxa(merged_physeq_ge_rel, Family=="Pasteurellaceae" )
## Print all ambiguous genus classifications
rownames(tax_table(merged_physeq_ge_pas))[grep(";", rownames(tax_table(merged_physeq_ge_pas)))]
## Subset to taxa included in ambiguous classifications
merged_physeq_ge_pas_amb <- subset_taxa(merged_physeq_ge_rel, grepl("Haemophilus", Genus) | grepl("Mannheimia", Genus) | grepl("Pasteurella", Genus) | grepl("Bibersteinia", Genus) | grepl("Necropsobacter", Genus) )

# Plots in full sample
#plot_composition(merged_physeq_ge_pas_amb)
#plot_composition(merged_physeq_ge_pas_amb, average_by = "Taxonomy_type")

# Extract 12 top relative abundance samples
head(sort(colSums(otu_table(merged_physeq_ge_pas_amb)), decreasing = T), n = 22)
for.plot.names <- names(head(sort(colSums(otu_table(merged_physeq_ge_pas_amb)), decreasing = T), n = 22))
for.plot.names <- substr(for.plot.names,1,nchar(for.plot.names)-1)

## Create plot

plot.dat <- subset_samples(merged_physeq_ge_pas_amb, 
                           sampleID %in% for.plot.names)

x <- plot.dat
abu <- abundances(x)

dfm <- psmelt(otu_table(abu, taxa_are_rows = TRUE))
names(dfm) <- c("Tax", "Sample", "Abundance")
dfm$Group <- meta(x)[["Taxonomy_type"]][match(as.character(dfm$Sample),
                                              sample_names(x))]
samp.names <- substr(sample_names(x),1,nchar(sample_names(x))-1)
sample.sort <- unique(rev(samp.names[order(abundances(x)["Haemophilus",])]))
dfm$Sample <- substr(dfm$Sample,1,nchar(dfm$Sample)-1) 
dfm$Sample <- factor(dfm$Sample, levels=sample.sort)
otu.sort <- taxa(x)
dfm$Tax <- factor(dfm$Tax, levels=otu.sort)

cols1 <- c(brewer.pal(10, "Spectral")[8],
           brewer.pal(10, "Spectral")[9],
           brewer.pal(10, "Spectral")[10],
           brewer.pal(10, "Spectral")[1],
           brewer.pal(10, "Spectral")[2],
           brewer.pal(10, "Spectral")[3],
           brewer.pal(10, "Spectral")[4],
           brewer.pal(10, "Spectral")[5],
           brewer.pal(10, "Spectral")[6])


Tax <- Sample <- Abundance <- NULL

# Provide barplot
dfm <- dfm %>% arrange(Tax)  # Show Taxs always in the same order
dfm$Tax <- factor(dfm$Tax, levels = c(unique(dfm$Tax[!grepl(";", dfm$Tax)]), unique(dfm$Tax[grepl(";", dfm$Tax)])))

p <- ggplot(dfm, aes(x=Sample, y=Abundance, fill=Tax)) +
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels=dfm$xlabel, breaks=dfm$Sample) +
  scale_discrete_manual(values=cols1, aesthetics="fill") 

# Name appropriately
p <- p + labs(y = "Relative Abundance")

# Rotate horizontal axis labels, and adjust
p <- p + theme(axis.text.x=element_text(angle=90, vjust=0.5,
                                        hjust=0))
p <- p + guides(fill=guide_legend(reverse=FALSE))


p <- p + facet_grid(Group~., drop = TRUE,
                    space = "free", scales = "free") 
p

ggsave(file.path(parathaaDir, "Pasteruellaceae_vignette.pdf"), width = 7, height = 3.5, units = "in") 

