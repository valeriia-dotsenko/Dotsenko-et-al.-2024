library(readxl)
library(DESeq2)
library(ggpubr)
library(scales)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library("kableExtra")
library(RGraphics)
library(gridExtra)
library("biomaRt")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(ggplotify)
library(reshape2)
library(plyr)
library(stats)
library(DescTools)
library(numform)
library("circlize")
library(ComplexHeatmap)
library(cowplot)
library(rlist)
library(export)



CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
CD <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")

#create theme for all plots
#create theme for all plots
mytheme <- theme_classic(base_family='sans') + theme(axis.text.x = element_text(size=6, color = "black"),
                                                     axis.text.y= element_text(size=6, color = "black"),
                                                     axis.title.y = element_text(size=7, color = "black"),
                                                     axis.title.x = element_text(size=7, color = "black"),
                                                     legend.text = element_text( size = 5, color = "black"),
                                                     legend.key.size = unit(5, 'mm'),
                                                     legend.position = "top",
                                                     legend.title = element_blank(),
                                                     plot.title = element_text(hjust = 0.5, size = 7,face = "bold"),
                                                     axis.line = element_line(colour = "black", linewidth=1/.pt),
                                                     panel.background = element_rect(fill = 'White'),
                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     strip.text.x = element_text(size = 5.5, colour = "black"))

palette_bio <- c("#00b347","#006729","#6500ff","#ff6500")
#-------------------------------------------------------------------------------panel A (PCA)-----------------------------------------------------------------------
set.seed(1)
vst_seobj <- vst(CD)
pca_obja <- plotPCA(vst_seobj, intgroup=c("label2"))
##### plotting prep
metadata <- colData(CD)
metadata$VHCrD_group <- NA
metadata[metadata$VHCrD >= 2.5,]$VHCrD_group <- "H (>= 2.5) "
metadata[metadata$VHCrD < 2.5 ,]$VHCrD_group <- "M (1.2-2.5) "
metadata[metadata$VHCrD < 1.2,]$VHCrD_group <- "L (<= 1.2) "

df <- pca_obja$data
df <- merge(df, metadata[,c("ID","VHCrD_group","label3")], by.x="name", by.y="ID")
df <- as.data.frame(df)
hulls <- plyr::ddply(df, "VHCrD_group", function(df) df[chull(df$PC1, df$PC2),])

#plotting
A <- ggplot()+
  geom_polygon(data = hulls,aes(PC1,PC2,group=VHCrD_group, fill=VHCrD_group), alpha = 0.25)+
  geom_point(data = df,aes(PC1, PC2,color=label2), size=1)+
  scale_color_manual(values=c("#00b347","#006729","#6500ff","#ff6500"))+
  scale_fill_manual(values=c("#ffd689","#ff89b2","#89edff"))+
  labs(x=sprintf("PC1: %s Variance", percent(pca_obja$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obja$plot_env$percentVar[2])) ) 
A <- A + mytheme + theme(legend.box="vertical", legend.margin=margin(),
                            plot.margin = margin(0, 1.5, 1, 1.5, "cm")) #top clockwise
A

#-------------------------------------------------------------------------------panel B (DEGs table)----------------------------------------------------------------
DEG1 <- CDres[[1]]
DEG2 <- CDres[[2]]
DEG3 <- CDres[[3]]

group_comp2 <- data.frame(`Comparisons` = names(CDres)[1:3],
                          `DEGs` = c(length(unique(DEG1[DEG1$Gene.type == "DE",]$ID)),length(unique(DEG2[DEG2$Gene.type == "DE",]$ID)),length(unique(DEG3[DEG3$Gene.type == "DE",]$ID))),
                          `Downregulated` = c(length(unique(DEG1[DEG1$Gene.type == "DE" & DEG1$log2FoldChange < 0.5,]$ID)) ,
                                              length(unique(DEG2[DEG2$Gene.type == "DE" & DEG2$log2FoldChange < 0.5,]$ID)),
                                              length(unique(DEG3[DEG3$Gene.type == "DE" & DEG3$log2FoldChange < 0.5,]$ID))),
                          `Upregulated` = c(length(unique(DEG1[DEG1$Gene.type == "DE" & DEG1$log2FoldChange > 0.5,]$ID)) ,
                                            length(unique(DEG2[DEG2$Gene.type == "DE" & DEG2$log2FoldChange > 0.5,]$ID)),
                                            length(unique(DEG3[DEG3$Gene.type == "DE" & DEG3$log2FoldChange > 0.5,]$ID))))
myt <- ttheme_default(base_size = 6, base_colour = "black", base_family = "",parse = T)
B <- tableGrob(group_comp2, rows = NULL,theme=myt)
as.ggplot(B)

#-----------------------------------------------------------------------------panel C (volcano plots)--------------------------------------------------------------------
DEG1 <- CDres[[1]]
DEG1$size = 0.5
DEG1$size[DEG1$Gene.type %in% c("DE")] = 0.75


#genes to highlight
data <- DEG1[DEG1$Gene.type=="DE", ] %>% top_n(5,log2FoldChange)
data <- rbind(data,DEG1[DEG1$Gene.type=="DE", ] %>% top_n(-5,log2FoldChange))
data <- rbind(data,DEG1[DEG1$Gene.type=="DE", ] %>% top_n(-5,padj))

C1 <- ggplot(DEG1, aes(x = log2FoldChange, y = -log10(padj), color=Gene.type)) + 
  geom_point(size=DEG1$size)+
  scale_color_manual(values=c("#64cdac","#d3d3d3"))+ 
  ylim(0,4)+ xlim(-2.5,2.5)+
  xlab("log2 fold change") + ylab("-log10(FDR) ")+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour="#858585", size=0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour="#858585", size=0.1) +
  ggtitle(unique(DEG1$Experiment))

C1 <- C1 + mytheme + theme(axis.title.x = element_blank())

DEG2 <- CDres[[2]]
DEG2$size = 0.5
DEG2$size[DEG2$Gene.type %in% c("DE")] = 0.75

#genes to highlight
data <- DEG2[DEG2$Gene.type=="DE", ] %>% top_n(5,log2FoldChange)

C2 <- ggplot(DEG2, aes(x = log2FoldChange, y = -log10(padj), color=Gene.type)) + 
  geom_point(size=DEG2$size)+
  scale_color_manual(values=c("#64cdac","#d3d3d3"))+ 
  xlab("log2FC") + ylab("-log10(FDR) ")+
  ylim(0,4)+ xlim(-2.5,2.5)+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour="#858585", size=0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour="#858585", size=0.1) +
  ggtitle(unique(DEG2$Experiment))+ 
  geom_text_repel(
    data = unique(data),
    aes(log2FoldChange, -log10(padj), label = external_gene_name),
    size = 2, # font size in the text labels
    colour="#858585",
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  )

C2 <- C2 + mytheme + theme(axis.title.y = element_blank())

DEG3 <- CDres[[3]]
DEG3$size = 0.5
DEG3$size[DEG3$Gene.type %in% c("DE")] = 0.75


C3 <- ggplot(DEG3, aes(x = log2FoldChange, y = -log10(padj), color=Gene.type)) + 
  geom_point(size=DEG3$size)+
  scale_color_manual(values=c("#64cdac","#d3d3d3"))+ 
  xlab("log2 fold change") + ylab("-log10(FDR) ")+
  ylim(0,4)+ xlim(-2.5,2.5)+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour="#858585", size=0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour="#858585", size=0.1) +
  ggtitle(unique(DEG3$Experiment))

C3 <- C3 + mytheme+ theme(axis.title.x = element_blank(), axis.title.y = element_blank())

C <- ggarrange(C1,C2,C3,
          ncol = 3, common.legend = T, align = "h")

C

#---------------------------------------------------------------------------panel D (Venn)-----------------------------------------------------------------------------------------------

# Chart all
D <- grobTree(venn.diagram(
  x = list(DEG1[DEG1$Gene.type=="DE",]$ID, 
           DEG3[DEG3$Gene.type=="DE",]$ID),
  category.names = c("PGCp VS GFDp" , "PGCp VS PGCd"),
  filename = NULL,
  fontfamily="sans", #counts font
  cex = 0.5, #counts size; 0.9 is 10pt
  cat.fontfamily =  "sans", #labels font
  cat.cex = 0.5, #labels size; 0.8 is appr 9pt
  #cat.default.pos = "outer", #label position
  cat.pos = c(285, 20),
  cat.dist = 0.01,
  col=c(alpha("#ff6500",0.5), alpha("#6500ff",0.5)),
  fill = c(alpha("#ff6500",0.3), alpha("#6500ff",0.3)),
  rotation.degree = 45,
  margin = 0.12))
D <- as.ggplot(D)
D

#--------------------------------------------------------------------------panel E()-------------------------------------------------------------------------------------------------
#correlation between the effect size of the comparison PGCp vs GFD and the effect size of the comparison PGCp and PGCd
df <- unique(DEG3[#DEG3$padj <= 0.05
  , c("ID", "log2FoldChange")])
df2 <- unique(DEG1[DEG1$ID %in% df$ID , c("ID", "log2FoldChange")])
df <- merge(df, df2, by = "ID")

diff1 <- df[df$log2FoldChange.x< 0 & df$log2FoldChange.y> 0,]
diff2 <- df[df$log2FoldChange.x> 0 & df$log2FoldChange.y< 0,]

G <- ggscatter(df, x = "log2FoldChange.y" , y = "log2FoldChange.x", size= 0.5,add = "reg.line",                                 # Add regression line
               conf.int = TRUE)+xlab("log2FC (PGCp VS GFDp)") + ylab("log2FC (PGCp VS PGCd)")+
  ylim(-1.5,2.5)+ xlim(-1.5,2.5) + mytheme


G
#----------------------------------------------------------------------------panel H (correlation heatmap)----------------------------------------------------------------------------------------------
normalized_counts <- counts(CD, normalized=TRUE) 


counts1<- normalized_counts[DEG1[DEG1$Gene.type == "DE",]$ID,]
counts3<- normalized_counts[DEG3[DEG3$Gene.type == "DE",]$ID,]

counts_comb <- rbind(counts1,counts3)
metadata <- metadata[order(metadata$label2, metadata$VHCrD),]
counts_comb <- counts_comb[,match(metadata$ID,colnames(counts_comb))]

cormat <- round(cor(counts_comb),2)


type1 <- metadata$label2


top = HeatmapAnnotation(Samples = anno_block(gp = gpar(fill = c("#00b347","#006729","#6500ff","#ff6500"), lty="blank")),
                        `VH:CrD` = anno_points(metadata$VHCrD, size = unit(0.5, "mm"), axis_param = list(side = "right", gp = gpar(fontsize = 3))),
                        #annotation_name_gp = gpar(fontsize = 3)
                        annotation_name_gp= gpar(fontsize = 5),
                        gap = unit(1, "mm"),
                        height = unit(0.5, "cm"))

col_type1 = list(`label` = c("GFDd" = "#00b347",
                             "GFDp" = "#006729",
                             "PGCd" = "#6500ff",
                             "PGCp" = "#ff6500")) 

row=rowAnnotation(#`VH:CrD` = anno_points(metadata$VHCrD, size = unit(0.5, "mm")),
                  #annotation_name_gp= gpar(fontsize = 3),
                  `label` = as.character(metadata$label2), col=col_type1, 
                  show_annotation_name = F,
                  annotation_legend_param = list(`label` = list(title_position = "topcenter",labels_gp = gpar(fontsize = 5),
                                                                title_gp = gpar(fontsize = 4),direction = "vertical",ncol = 1,#nrow = 1,
                                                                legend_width=unit(1,"cm"),grid_height = unit(0.2, "cm"))),
                  
                  gap = unit(0.5, "mm"),
                  #width= unit(2, "mm") 
                  simple_anno_size = unit(2, "mm"))
#simple_anno_size = unit(2, "mm"))

#create color scheme for heatmap
newcolors <- c(#"#67001F",
  "#B2182B",
  #"#D6604D", 
  #"#F4A582",
  #"#FDDBC7", 
  "#F7F7F7", #white 
  #"#D1E5F0", 
  #"#92C5DE",
  #"#4393C3", 
  "#2166AC"  #, 
  #"#053061"
)

ht_list = Heatmap(cormat, name = "Pearson\nCorrelation", top_annotation = top, left_annotation = row,
                  show_column_names = F,
                  show_row_names = F, 
                  col=newcolors,
                  column_order= metadata$ID,
                  row_order= metadata$ID,
                  cluster_rows = F,  
                  cluster_columns = F, 
                  column_split = factor(type1, levels = c("GFDd","GFDp", "PGCd","PGCp")),
                  column_gap = unit(0.5,"mm"),
                  row_split = factor(type1, levels = c("GFDd","GFDp", "PGCd","PGCp")),
                  row_title = NULL,
                  row_gap = unit(2,"mm"),
                  column_title_gp = gpar(fontsize = 0.5, col="White"),
                  heatmap_legend_param = list(
                    at = c(0:1),
                    labels = c(0:1),
                    title = "Pearson\nCorrelation",
                    labels_gp = gpar(fontsize = 5),
                    title_gp = gpar(fontsize = 5),direction = "horizontal",
                    legend_width=unit(1,"cm"),grid_height = unit(2, "mm"),
                    title_position = "topcenter"))#+   row
draw(ht_list,merge_legend = FALSE,heatmap_legend_side="bottom", annotation_legend_side="right", legend_grouping = "original" #, row_gap = unit(0.5, "mm"),auto_adjust = T
     )
H <- grid.grabExpr(draw(ht_list,merge_legend = FALSE,heatmap_legend_side="bottom", annotation_legend_side="right", legend_grouping = "original", row_gap = unit(0.5, "mm")))
#hm1

#----------------------------------------------------------------------------combine all together ---------------------------------------------------------------------------------
filler <- ggplot() + theme(panel.background	= element_rect(fill = "white"))
p0 <- ggarrange(as.ggplot(B),filler,
                ncol = 2,widths = c(1, 0.5))
p0
p1 <- ggarrange(A, p0, ncol = 2,   widths=c(1,1))
p1


#p2 <- ggarrange(C,filler,D, filler, ncol = 4,  align = "h", widths = c(3,0.75,1,0.75))
#p2

p3 <- ggarrange(D, G, H, ncol = 3,  align = "h", widths = c(1,1,1))
p3



p <- arrangeGrob(p1, C,p3,
                    ncol = 1, nrow = 3)
as_ggplot(p)


# Add labels to the arranged plots
p <- as_ggplot(p) +    
  draw_plot_label(label = c("a", "b", "c","d","e","f"), size = 11,
                  x = c(0, 0.5, 0,  0,   0.33,   0.66), 
                  y = c(1, 1,  0.66, 0.35, 0.35, 0.35)) # transform to a ggplot


p
#ggexport(p,filename = "./output_data/Figure_1_1.pdf") # Figure panels should be prepared at a minimum resolution of 300 dpi and saved at a maximum width of 180 mm (7,08661 inch).  
graph2ppt(file="./output_data/Figure_1.pptx", width=7.08, height=7.28)

#pdf(file = "./output_data/Figure_1.pdf",   # The directory you want to save the file in
#    width = 4, # The width of the plot in inches
#    height = 4) # The height of the plot in inches
