### load in the data
require(docopt)

'Usage:
  real_data_plots.R [--OralV4V5 <directory oral V4V5 RDS> --MineV1V3 <directory Mine V1V3 RDS> --MineV4V5 <directory Mine V4V5 RDS> --output <output directory>]
  
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

message(opts$MineV4V5)
message(opts$OralV4V5)
message(opts$MineV1V3)
message(opts$output)

mine_taxa_v4v5 <- readRDS(file=file.path(opts$MineV4V5, "taxa_bars.RDS"))

### oral data
oral_brays <- readRDS(file=file.path(opts$OralV4V5, "bray_plots_sub.RDS"))
oral_average <- readRDS(file=file.path(opts$OralV4V5, "average_plots.RDS"))



## top
bray_panel <- oral_brays[[1]] + coord_fixed() + oral_brays[[5]] + coord_fixed() & theme(legend.position = "bottom")
bray_panel <- (bray_panel + plot_layout(guides="collect")) + plot_annotation(title="A)")

oral_taxa_panel <- oral_average[[1]] + oral_average[[2]] + plot_annotation(title="B)")
oral_taxa_panel

oral_panel <- bray_panel | oral_taxa_panel
oral_panel + plot_annotation(theme=theme(plot.background = element_rect(fill="grey")))


final_fig <- oral_panel / mine_taxa_v4v5[[1]] + plot_annotation(tag_levels = "A")
final_fig
ggsave(plot = final_fig, filename = file.path(opts$output, "Figure_4.pdf"), height=6, width=13.8)
ggsave(plot = final_fig, filename = file.path(opts$output, "Figure_4.png"), height=6, width=13.8)

### Supplemental for v13 sediment data

mine_v13_bray <- readRDS(file=file.path(opts$MineV1V3, "bray_plots.RDS"))
mine_v13_taxa <- readRDS(file=file.path(opts$MineV1V3, "taxa_bars.RDS"))

bray_panel <- wrap_plots(mine_v13_bray) & theme(legend.position = "bottom") & coord_fixed()
bray_panel <- bray_panel + plot_layout(guides="collect")

taxa_panel <- mine_v13_taxa[[1]]
taxa_panel


final <- bray_panel | taxa_panel
final
ggsave(plot=final, filename=file.path(opts$output, "Supplemental_Mine_V1V3.pdf"), height=6, width=13.8)
ggsave(plot=final, filename=file.path(opts$output, "Supplemental_Mine_V1V3.png"), height=6, width=13.8)

### V4V5 full supplemental
mine_v45_bray <- readRDS(file=file.path(opts$MineV4V5, "bray_plots.RDS"))

bray_panel <- wrap_plots(mine_v45_bray[-6]) & theme(legend.position = "bottom") & coord_fixed()
bray_panel <- bray_panel + plot_layout(guides="collect")

taxa_panel <- mine_taxa_v4v5[[2]]
taxa_panel


final <- bray_panel | taxa_panel
final
ggsave(plot=final, filename=file.path(opts$output, "Supplemental_mine_V4V5.pdf"), height=6, width=13.8)
ggsave(plot=final, filename=file.path(opts$output, "Supplemental_mine_V4V5.png"), height=6, width=13.8)
