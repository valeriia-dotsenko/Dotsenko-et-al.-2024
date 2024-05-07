library(DESeq2)
library(dplyr)
library(scales)
library(ggpubr)
library(gridExtra)
library(ggplotify)
library(ggrepel)
library(VennDiagram)
library("circlize")
library(ComplexHeatmap)
library(cowplot)
library(export)

ORGraw <- readRDS("./input_data/hDuo ZED1277 DESeqDataSet object.rds")
ORGres <- readRDS('./input_data/hDuo ZED1277 DEG results.rds')


#create theme for all plots
mytheme <- theme_classic() + theme(axis.text.x = element_text(size=8, color = "black"),
                                   axis.text.y= element_text(size=8, color = "black"),
                                   axis.title.y = element_text(size=10, color = "black"),
                                   axis.title.x = element_text(size=10, color = "black"),
                                   legend.text = element_text( size = 6, color = "black"),
                                   legend.key.size = unit(5, 'mm'),
                                   legend.position = "top",
                                   legend.title = element_blank(),
                                   plot.title = element_text(hjust = 0.5, size = 10,face = "bold"),
                                   axis.line = element_line(colour = "black", size=1/.pt)
)
palette_ORG <- c("#FFBD35","#3FA796")
#-------------------------------------------------------------------------------panel A (PCA)-----------------------------------------------------------------------
set.seed(1)
vst_seobj <- vst(ORGraw)
pca_obja <- plotPCA(vst_seobj, intgroup=c("Treatment"))
##### plotting prep
metadata <- colData(ORGraw)

df <- pca_obja$data %>% filter(Treatment %in% c("mock", "50uM_ZED1227"))

#plotting
A <- ggplot()+
  geom_point(data = df,aes(PC1, PC2,color=Treatment), size=1)+
  scale_color_manual(values=palette_ORG)+
  labs(x=sprintf("PC1: %s Variance", percent(pca_obja$plot_env$percentVar[1])),
       y=sprintf("PC2: %s Variance", percent(pca_obja$plot_env$percentVar[2])) ) 
A <- A + mytheme + theme(legend.box="vertical", legend.margin=margin(),
                         plot.margin = margin(0, 1.5, 1, 1.5, "cm")) #top clockwise
A
#-------------------------------------------------------------------------------panel B (DEGs table)----------------------------------------------------------------
DEG1 <- ORGres[[1]]


group_comp2 <- data.frame(`Comparisons` = "50uM ZED1227 VS mock",
                          `DEGs` = length(unique(DEG1[DEG1$Gene.type == "DE",]$ID)),
                          `Downregulated` = length(unique(DEG1[DEG1$Gene.type == "DE" & DEG1$log2FoldChange < 0.5,]$ID)),
                          `Upregulated` = length(unique(DEG1[DEG1$Gene.type == "DE" & DEG1$log2FoldChange > 0.5,]$ID)))
myt <- ttheme_default(base_size = 6, base_colour = "black", base_family = "",parse = T)
B <- tableGrob(group_comp2, rows = NULL,theme=myt)
as.ggplot(B)


#-----------------------------------------------------------------------------panel C (volcano plots)--------------------------------------------------------------------
DEG1$size = 0.5
DEG1$size[DEG1$Gene.type %in% c("DE")] = 0.75


#genes to highlight
data <- DEG1[DEG1$Gene.type=="DE", ] %>% top_n(5,log2FoldChange)
data <- rbind(data,DEG1[DEG1$Gene.type=="DE", ] %>% top_n(-5,log2FoldChange))
data <- rbind(data,DEG1[DEG1$Gene.type=="DE", ] %>% top_n(-5,padj))

C1 <- ggplot(DEG1, aes(x = log2FoldChange, y = -log10(padj), color=Gene.type)) + 
  geom_point(size=DEG1$size)+
  scale_color_manual(values=c("#64cdac","#d3d3d3"))+ 
  #ylim(0,4)+ xlim(-2.5,2.5)+
  xlab("log2 fold change") + ylab("-log10(FDR) ")+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour="#858585", linewidth=0.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour="#858585", linewidth=0.1) +
  #ggtitle(unique(DEG1$Experiment))+ 
  geom_text_repel(
    data = unique(data),
    aes(log2FoldChange, -log10(padj), label = external_gene_name),
    size = 2, # font size in the text labels
    colour="#858585",
    box.padding = unit(0.1, "lines"),
    point.padding = unit(0.1, "lines")
  )

C1 <- C1 + mytheme + theme(axis.title.x = element_blank())
C1

#-------------------------------------------------------------------------D Venn diagram-------------------------------------------------------------------------
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')


p1 <- grobTree(venn.diagram(
  x = list(
    CDres[["PGCd VS GFDd"]][CDres[["PGCd VS GFDd"]]$Gene.type=="DE",]$ID,
    CDres[["PGCp VS PGCd"]][CDres[["PGCp VS PGCd"]]$Gene.type=="DE",]$ID,
    DEG1[DEG1$Gene.type=="DE",]$ID),
  category.names = c(#"PGCp VS GFD",
    "PGCd VS GFDd","PGCp VS PGCd", "50uM ZED1227 VS mock"),
  filename = NULL,
  fontfamily="sans", #counts font
  cex = 0.6, #counts size; 0.9 is 10pt
  cat.fontfamily =  "sans", #labels font
  cat.cex = 0.5, #labels size; 0.8 is appr 9pt
  #cat.default.pos = "outer", #label position
  #cat.pos = c(200, 20),
  cat.dist = 0.01,
  #col=c(alpha("#ff6500",0.5), alpha("#8267BE",0.5)),
  #fill = c(alpha("#ff6500",0.3), alpha("#8267BE",0.3)),
  rotation.degree = 45,
  margin = 0.12))

D <- ggarrange(as.ggplot(p1))
D  

#--------------------------------------------------------------------------panel E (heatmap)------------------------------------------------------------------------------

normalized_counts <- counts(ORGraw, normalized=TRUE)
DEG1 <- DEG1[DEG1$Gene.type=="DE",]
counts1<- normalized_counts[DEG1$ID,]
#sort the data by logFC
DEG1 <- DEG1[order(DEG1$log2FoldChange, decreasing=T),]
counts1<-counts1[match(rownames(counts1),DEG1$ID),]

#change id to names
rownames(counts1) <- DEG1$external_gene_name[match(rownames(counts1), DEG1$ID)]
  
#scale the data
m.data_sc<-t(scale(t(counts1)))
  
#sort columns
metadata1 <- metadata[metadata$label %in% c("M", "Z"),]
m.data_sc <- m.data_sc[,metadata1$Sample_Name]
#create annotation 1(Sample type)
type1<- metadata1$label

top = HeatmapAnnotation(Samples = anno_block(gp = gpar(fill = palette_ORG, lty="blank"),
                                               labels = unique(type1), 
                                               labels_gp = gpar(col = "white", fontsize = 10)),
                          gap = unit(2, "mm"),
                          annotation_height = unit(0.3, "cm"))
  #height = unit(0.3, "cm"))
  
  #create annotation 2 (log2FC)
  col_fun = colorRamp2(c(max(DEG1$log2FoldChange), 
                         mean(DEG1$log2FoldChange), min(DEG1$log2FoldChange)), 
                       c("#053061", "#F7F7F7","#67001F" ))
  
  #create color scheme for heatmap
newcolors <- c(#"#67001F", #red
    "#B2182B",
    #"#D6604D", 
    #"#F4A582",
    "#FDDBC7", 
    "#F7F7F7", #white 
    "#D1E5F0", 
    #"#92C5DE",
    "#4393C3", 
    #"#2166AC", 
    "#053061") #blue
ht_list = Heatmap(m.data_sc, 
                    name = "Row Z-Score of Expression",
                    top_annotation = top, #ha1,
                    show_column_names = FALSE, 
                    col=newcolors,
                    column_title_gp = gpar(fontsize = 1, col="White"),
                    column_order = metadata1$Sample_Name,
                    row_names_side = "left",
                    #row_order = order,
                    cluster_rows = T, 
                    show_row_dend = F,
                    row_names_gp = gpar(fontsize = 6), 
                    column_split = factor(type1, levels = unique(type1)),
                    column_gap = unit(0.5,"mm"),
                    heatmap_legend_param = list(
                      at = c(-1:4),
                      labels = c(-1:4),
                      title = "Row Z-Score",
                      labels_gp = gpar(fontsize = 6),
                      title_gp = gpar(fontsize = 6),
                      direction = "horizontal",
                      legend_width=unit(2.1,"cm"),
                      grid_height = unit(0.2, "cm"),
                      title_position = "topcenter"))#+row
  
E <- grid.grabExpr( draw(ht_list, merge_legend = F, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", auto_adjust = FALSE))



#--------------------------------------------------------------------------Combine----------------------------------------------------------------------------------
#p0 <- ggarrange(as.ggplot(B),D,
#                ncol = 2,widths = c(1, 0.5))
#p0

p0 <- ggarrange(B, C1, nrow = 2)
p0
p1 <- ggarrange(A, p0,D, E, ncol = 2, nrow = 2,align = "h", widths = c(1,1))
p1

p <- p1 +    
  draw_plot_label(label = c("a", "b", "c","d","e"), size = 11,
                  x = c(0, 0.5, 0.5,  0.0,   0.5), 
                  y = c(1, 1,  0.75, 0.5, 0.5)) # transform to a ggplot

p
graph2ppt(file="./output_data/Figure_S5.pptx", width=7.08, height=7.28)


