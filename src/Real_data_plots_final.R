### load in the data
require(docopt)

'Usage:
  Final_PLOT.R [--OralV4V5 <directory oral V4V5 RDS> --MineV1V3 <directory Mine V1V3 RDS> --MineV4V5 <directory Mine V4V5 RDS> --output <output directory>]
  
Options:
  --OralV4V5 directory of oral data
  --MineV1V3 directory of mine v1v3 data
  --MineV4V5 directory of mine v4v5 data
  --output output directory
  ]' -> doc

opts <- docopt(doc)




get_variance_by_tool <- function(ps1.com, level, remove_unknown=FALSE){
  
  if(remove_unknown){
    ps1.com <- subset_taxa(ps1.com, eval(as.name(level)) != "unassigned")
  }
  metadata <- as(sample_data(ps1.com), "data.frame")
  
  ps1_Genus <- aggregate_taxa(ps1.com, level = level)
  ps1_Genus <- transform(ps1_Genus, transform = "compositional")

  res <- adonis2(distance(ps1_Genus, method="bray") ~ Taxonomy_type, data=metadata)
  
  return(res)
}

library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(vegan)
library(ggplot2)
library(patchwork)
library(pals)

message(opts$MineV4V5)
#opts$MineV4V5 <- "~/Desktop/Parathaa_fig3/MINE_V4V5/"
message(opts$OralV4V5)
#opts$OralV4V5 <- "~/Desktop/Parathaa_fig3/ORAL/"
message(opts$MineV1V3)
#opts$MineV1V3 <- "~/Desktop/Parathaa_fig3/MINE_V1V3/"
message(opts$output)
#opts$output <- "~/Desktop/Parathaa_fig3/Update_Sep_28/"

mine_taxa_v4v5 <- readRDS(file=file.path(opts$MineV4V5, "taxa_bars.RDS"))

### oral data
oral_brays <- readRDS(file=file.path(opts$OralV4V5, "bray_plots_sub.RDS"))
oral_average <- readRDS(file=file.path(opts$OralV4V5, "average_plots.RDS"))

color_hex_codes <- stepped2(n=10)

mine_v4v5_scatter <- readRDS(file=file.path(opts$MineV4V5, "scatter_genus_abundance.RDS"))
mine_v1v3_scatter <- readRDS(file=file.path(opts$MineV1V3, "scatter_genus_abundance.RDS"))

## fix up v1v3 scatter
mine_v1v3_keepers <- names(which(table(mine_v1v3_scatter$data$Phylum) > 3))

mine_v1v3_scatter$data <- mine_v1v3_scatter$data %>% mutate(filter_Phylum=case_when(
  Phylum %in% mine_v1v3_keepers ~ Phylum,
  !Phylum %in% mine_v1v3_keepers ~ "Other"
))

mine_v1v3_scatter$data$filter_Phylum[which(mine_v1v3_scatter$data$Phylum=="Unknown")] <- "Unknown"
mine_v1v3_scatter$data$filter_Phylum <- factor(mine_v1v3_scatter$data$filter_Phylum)
mine_v1v3_scatter$data$filter_Phylum <- relevel(mine_v1v3_scatter$data$filter_Phylum, "Other")
mine_v1v3_scatter$data$filter_Phylum <- relevel(mine_v1v3_scatter$data$filter_Phylum, "Unknown")

mine_v1v3_scatter <- mine_v1v3_scatter + scale_color_manual(values=c("grey", "#899499", color_hex_codes), name="Phylum")

## fix up v4v5 scatter
mine_v4v5_scatter$data$filter_Phylum[which(mine_v4v5_scatter$data$Phylum=="Unknown")] <- "Unknown"
mine_v4v5_scatter$data$filter_Phylum <- factor(mine_v4v5_scatter$data$filter_Phylum)
mine_v4v5_scatter$data$filter_Phylum <- relevel(mine_v4v5_scatter$data$filter_Phylum, "Other")
mine_v4v5_scatter$data$filter_Phylum <- relevel(mine_v4v5_scatter$data$filter_Phylum, "Unknown")

mine_v4v5_scatter <- mine_v4v5_scatter + scale_color_manual(values=c("grey", "#899499", color_hex_codes), name="Phylum")
mine_v4v5_scatter <- mine_v4v5_scatter + ggtitle("Sediment V4V5")
test <- ggplot_build(mine_v4v5_scatter)
test$data[[3]]$x <- -4

test2 <- patchwork::wrap_ggplot_grob(ggplot_gtable(test))


mine_scatters <- (test2 + ggtitle("Sediment V4V5")) | (mine_v1v3_scatter + ggtitle("Sediment V1V3"))
mine_scatters

ggsave(plot=mine_scatters, filename=file.path(opts$output, "mine_scatters.png"), height=6, width=13.8)


oral_plots <- (oral_brays[[1]] + coord_fixed() + theme(legend.position = "bottom")) | (oral_brays[[5]] + coord_fixed() & 
  theme(legend.position = "bottom"))


oral_plots <- (oral_brays[[1]] + coord_fixed() + theme(legend.position = "bottom")) | (plotList[[1]] + coord_fixed() & 
                                                                                         theme(legend.position = "bottom"))
ggsave(plot=oral_plots, filename=file.path(opts$output, "oral_brays.png"), height=6, width=13.8)

oral_average_test <- oral_average[[2]]
oral_average_test$data <- oral_average_test$data %>% filter(Tax!="Unknown")


oral_taxa_panel <- oral_average[[1]] | (oral_average_test) & theme(legend.justification = "left")

ggsave(plot=oral_taxa_panel, filename=file.path(opts$output, "oral_taxa_panel.png"), height=6, width=13.8)


### Mine supplemental
mine_taxa_v1v3 <- readRDS(file=file.path(opts$MineV1V3, "taxa_bars.RDS"))

mine_taxa_v1v3[[1]] <- mine_taxa_v1v3[[1]] + ggtitle("Sediment V1V3")

(mine_taxa_v1v3[[1]])/(mine_taxa_v4v5[[1]] + ggtitle("Sediment V4V5"))/(mine_taxa_v4v5[[2]] + ggtitle("Sediment V4V5"))
