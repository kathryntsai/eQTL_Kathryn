library(reshape2)
library(dplyr)
library(readxl)
library(gridExtra)
library(ggplot2)
library(ggthemes)
theme_set(theme_stata())

sc_eqtl <- read_xlsx("/Users/kathryntsai/OneDrive\ -\ Villanova\ University/College/2018-2019/Summer\ 2019/TFs_eQTLs_Research/RProjects/Project2/IMPACT/NIHMS76345-supplement-Supplementary_table_2_vanderWijst.xlsx", skip=1)

# How many SNP-gene associations exist in both single cell and bulk datasets?

g <-
  merge(
    sc_eqtl,
    franke_cis_data,
    by = c("SNP", "Gene"),
    all.x = F,
    all.y = F
  )

dim(g)
# 317 overlap

# Of these overlapping associations, are they more or less cell type specific than the associations that can 
# only be detected in the single cell data?

gg9 <- melt(g, 
           id.vars=c("SNP", "Gene"), 
           #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
           measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
      ) %>% 
      ggplot() + geom_histogram(alpha=.8, aes(x=-log(value), fill=variable)) + ggtitle("P-values of eQTLs, grouped by cell type") # possibly take out position="identity" 

gg10 <- melt(g, 
           id.vars=c("SNP", "Gene"), 
           #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
           measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of each cell type at each p-value") # possibly take out position="identity" 

gg11 <- melt(g, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(alpha=.8, aes(x=value, fill=variable)) + ggtitle("# eQTLs Significant at FDR < 0.05 by Cell Type") # possibly take out position="identity" 

gg12 <- melt(g, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of Cell Types eQTLs Significant at FDR < 0.05") # possibly take out position="identity" 


library(gridExtra)
grid.arrange(gg1, gg2, gg3, gg4)

# No PBMCs

gg13 <- melt(g, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(alpha=.8, aes(x=-log(value), fill=variable)) + ggtitle("P-values of eQTLs, grouped by cell type") # possibly take out position="identity" 

gg14 <- melt(g, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of each cell type at each p-value") # possibly take out position="identity" 

gg15 <- melt(g, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(alpha=.8, aes(x=value, fill=variable)) + ggtitle("# eQTLs Significant at FDR < 0.05 by Cell Type") # possibly take out position="identity" 

gg16 <- melt(g, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of Cell Types eQTLs Significant at FDR < 0.05") # possibly take out position="identity" 
grid.arrange(gg1, gg2, gg3, gg4)

#================

gg1 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(alpha=.8, aes(x=-log(value), fill=variable)) + ggtitle("P-values of eQTLs, grouped by cell type - all sc-eQTLs") # possibly take out position="identity" 

gg2 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of each cell type at each p-value - all sc-eQTLs") # possibly take out position="identity" 

gg3 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(alpha=.8, aes(x=value, fill=variable)) + ggtitle("# eQTLs Significant at FDR < 0.05 by Cell Type - all sc-eQTLs") # possibly take out position="identity" 

gg4 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of Cell Types eQTLs Significant at FDR < 0.05") # possibly take out position="identity" 


# No PBMCs

gg5 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(alpha=.8, aes(x=-log(value), fill=variable)) + ggtitle("P-values of eQTLs, grouped by cell type - all sc-eQTLs, no PBMCs") # possibly take out position="identity" 

gg6 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_histogram(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of each cell type at each p-value - all sc-eQTLs, no PBMCs") # possibly take out position="identity" 

gg7 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(alpha=.8, aes(x=value, fill=variable)) + ggtitle("# eQTLs Significant at FDR < 0.05 by Cell Type - all sc-eQTLs, no PBMCs") # possibly take out position="identity" 

gg8 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            measure.vars = c("T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "cMono__1")
            #measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "cMono")
) %>% 
  ggplot() + geom_bar(position="fill", alpha=.8, aes(x=value, fill=variable)) + ggtitle("Proportions of Cell Types eQTLs Significant at FDR < 0.05 - all sc-eQTLs, no PBMCs") # possibly take out position="identity" 
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, ncol=2)


### Making a Heat Map
library(pheatmap)
library(RColorBrewer)
library(viridis)

gg17 <- melt(sc_eqtl, 
            id.vars=c("SNP", "Gene"), 
            #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
            measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "Mono")
)

mat_col <- data.table(gg17$variable)
mat_colors <- list(group = brewer.pal(length(unique(gg17$variable)), "Set1"))
names(mat_colors$group) <- unique(gg17$variable)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat, n = 11)


draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

pheatmap(
  mat               = -log10(sc_eqtl[, c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "Mono")]),
  color             = inferno(length(mat_breaks) -1),
  breaks            = mat_breaks,
  border_color      = F,
  cellwidth = 15,
  cellheight= 3,
  legend = TRUE,
  show_colnames     = T,
  show_rownames     = F,
  treeheight_row = 0,
  treeheight_col = 0,
  annotation_names_col    =  data.table(group = unique(gg17$variable)),
  annotation_legend = TRUE,
  #annotation_names_row = data.table(group = paste0(gg17$SNP, "-", gg17$Gene)),
  annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 3,
  main              = "Quantile Color Scale with Log Values"
)

# https://slowkow.com/notes/pheatmap-tutorial/
  
# library(scales)
# library(plyr)
# gg17 <- melt(sc_eqtl, 
#              id.vars=c("SNP", "Gene"), 
#              #measure.vars = c("PBMC__1","T CD4+__1", "T CD8+__1", "NK__1", "B__1", "DC__1", "Mono__1", "cMono__1", "ncMono__1")
#              measure.vars = c("PBMC","T CD4+", "T CD8+", "NK", "B", "DC", "Mono")
# )
# gg17 <- ddply(gg17, .(variable), transform, rescale = rescale(-log(value)))
# gg17$eQTL <- paste0(gg17$SNP, gg17$Gene)
# base_size <- 9
# gg17_plot <- ggplot(gg17, aes(variable, eQTL)) + 
#   geom_tile(aes(fill = rescale), colour = "white") + 
#   #scale_fill_gradient(low = "white", high = "steelblue") + 
#   theme_grey(base_size = base_size) + 
#   labs(x = "", y = "") + 
#   scale_x_discrete(expand = c(0, 0)) +     
#   scale_y_discrete(expand = c(0, 0)) + 
#   theme_economist() + 
#   theme(axis.text.y=element_text(size=3), 
#         axis.text.x=element_text(size=8),
#         axis.title=element_text(size=8,face="bold"), 
#         legend.title = element_text(size = 3), 
#         legend.text = element_text(size = 3)
#        )
# theme(legend.position = "none", axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))


# Do we even need sc-eQTLs if we can use IMPACT to identify the cell type?

# Since Pfizer eQTLs are found in whole blood, can we use IMPACT to infer cell-type-specificity? Since eQTLgen eQTLs are found in 
# peripheral blood, can we use IMPACT to infer cell-type-specificity?

# After learning from the annotation analyses above, can we predict if a given SNP will participate in a cis or trans eQTL? 
# Consider 3 groups of SNPs: only cis, only trans, both cis and trans (with two separate genes, of course)

# What are the models of these datasets? Do all include expression PCs? Traditionally, expression PCs don’t take away signal from cis eQTLs but they might from trans eQTLs. Need to know data details, e.g. how many donors? How many samples? How many cells?

# proportion of pbmcs are different cells