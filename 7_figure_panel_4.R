library(readxl)
library(VennDiagram)
library(reshape2)
library(rlist)
library(numform)
library(ggpubr)
library(mGSZ)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(dplyr)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(rstatix)
library(export)


human =  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")


#open dadplyr#open data sets for comparisons
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
CDraw <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")

metadata <- colData(CDraw)


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
                                                     strip.text.x = element_text(size = 5.5, colour = "black"),
                                                     axis.title.x.top = element_text(color='#4daf4a',size=7), 
                                                     axis.text.x.top = element_text(color='#4daf4a',size=6))

palette_bio <- c("#00b347","#006729","#6500ff","#ff6500")
palette_3 <- c("#00b347","#6500ff","#ff6500")
#-----------------------------------------------------------------------panel A PPARG expression --------------------------------------------------------------------------
IDs <-c('ENSG00000198959', 'ENSG00000132170', 'ENSG00000007171') # TGM2 & PPARG & iNOS

#prepare data
normalized_counts <- counts(CDraw, normalized=TRUE)

melt<-melt(normalized_counts[rownames(normalized_counts) %in% IDs,])
#melt$variable <- rownames(melt)
melt$log10_Expr <- log10(melt$value+1)
melt <-merge(melt, as.data.frame(metadata[,c("ID","label2", "IELs")]) ,by.x="Var2", by.y="ID")
group_comp <- data.frame(`Sample Group` = c("PGCp","PGCd","PGCp"),
                         `Control Group` = c("GFDp","GFDd","PGCd"),
                         `Group Comparisons` = c("PGCp VS GFDp","PGCd VS GFDd","PGCp VS PGCd"))

anno<- data.frame()
for(j in IDs){
for(i in 1:nrow(group_comp)){
  DEG1 <- CDres[[i]]
  anno1 <-subset(DEG1, DEG1$ID %in% j, select=c(external_gene_name,baseMean,padj)) 
  anno1$group1 <- group_comp$`Sample.Group`[i]
  anno1$group2 <- group_comp$`Control.Group`[i]
  anno <- rbind(anno,anno1)
}}

anno %>% 
  mutate(p_new = ifelse(padj > 0.01, c(paste("italic('P')~`=", f_num(padj,2), "`")), padj))%>% 
  mutate(p_new = ifelse(padj < 0.01, c(paste("italic('P')~`=", f_num(padj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(padj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->anno

anno$y_pos <- c(1,1,1,1.5,1.6,1.7,2,1.8,1.9)

A <- ggboxplot(as.data.frame(melt[melt$Var1 == "ENSG00000132170",]), x = "label2", y = "log10_Expr", outlier.colour = NA, order = c("GFDd", "GFDp", "PGCd","PGCp"),#c("GFD", "PGCd","PGCp"), #unique(metadata$label3), #c("GFDd", "GFDp", "PGCd","PGCp")
          ylab = "PPARG mRNA \n log10(Expression, counts)",# yscale="log2",
          palette = palette_bio, # c("#ffd700","#ffad00", "#ff934d","#ff6500"),
          fill="label2" #, #facet.by = c("Name"),ncol = 3,scales = "free_y",
          # title = "TGM2 mRNA expression"
)+
  #facet_grid(Marker~Name, scales="free", space="free_x")+
  geom_jitter(size=0.5, alpha=0.5)+
  #scale_y_continuous(limits=c(0.7, 3.0))
  labs(fill = element_blank())+
  #theme(strip.text.x = element_text(size = 6))+
  geom_signif(data=anno[anno$external_gene_name == "PPARG",],
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3) + mytheme +
  theme(axis.title.x = element_blank(),
        axis.text.y=element_blank(), legend.position='none')+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

A

#-----------------------------------------------------------------------panel B PPAG GSZ---------------------------------------------------------


SampleGSZ <- CDres[["SampleGSZ"]]
groupGSZ <- CDres[["groupGSZ"]]

df3 <- merge(SampleGSZ[,c("Samples","KEGG_PPAR_SIGNALING_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION",
                            "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "BIOCARTA_NO2IL12_PATHWAY")], as.data.frame(metadata[,c("ID","label2", "label3", "IELs")]), by.x = "Samples", by.y = "ID")


groupGSZ1 <- groupGSZ[groupGSZ$gene.sets == "KEGG_PPAR_SIGNALING_PATHWAY",]
groupGSZ1$y_pos <- c(10,11,12)

B <- ggboxplot(df3, x = "label3", y = "KEGG_PPAR_SIGNALING_PATHWAY", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
               ylab = "KEGG PPAR SIGNALING \n PATHWAY gene set GSZ",
               palette = palette_3, fill="label3") +
  geom_jitter(size=0.5, alpha=0.5)+
  scale_x_discrete(labels=c("GFDd+p","PGCd","PGCp"))+
  geom_hline(yintercept=c(mean(SampleGSZ$KEGG_PPAR_SIGNALING_PATHWAY)), linetype="dashed", color = "grey60")+
  geom_signif(data=groupGSZ1,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position= y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)+ mytheme +
  theme(axis.title.x = element_blank(),
        #axis.text.y=element_blank(), 
        legend.position='none')+
   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
B


#-------------------------------------------------------------------------panel C TGM2 VS PPARG---------------------------------------------
df <- as.data.frame(t(normalized_counts[rownames(normalized_counts) %in% IDs,]))
df$ID <- rownames(df)
df <- merge(df, as.data.frame(metadata[,c("ID","label2", "IELs")]),by.x="ID", by.y="ID") 


cor.test.result <- cor.test(df$ENSG00000198959, 
                            df$ENSG00000132170 , method = "pearson") #alternative hypothesis: true correlation is not equal to 0, two-tailed test

correlations <- data.frame(correlation_coefficient = cor.test.result$estimate, p.value = cor.test.result$p.value)

correlations <- correlations %>% mutate_all(as.numeric) %>%
  mutate(label = ifelse(`p.value` > 0.01, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `=", f_num(`p.value`,2), "`")), `p.value`))%>% 
  mutate(label = ifelse(`p.value` < 0.01, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `=", f_num(`p.value`,3), "`")), label)) %>%
  mutate(label = ifelse(`p.value` < 0.001, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `", "<.001", "`")),label))%>% 
  mutate(p_full = ifelse(`p.value` > 0.01, round(p.value,2), `p.value`))%>% 
  mutate(p_full = ifelse(`p.value` < 0.01, round(p.value,3), p_full)) %>%
  mutate(p_full = ifelse(`p.value` < 0.001, formatC(p.value, format = "e", digits = 2),p_full))

correlations$x <- 20
correlations$y <- 40

C1 <- ggscatter(as.data.frame(df), x = "ENSG00000132170", y = "ENSG00000198959",color = "label2",palette = palette_bio, size = 0.5) + #c(PGCd = "#6500ff", PGCp = "#ff6500")
  labs(x = "PPARG mRNA \n log10(Expression, counts)", y= "TGM2 mRNA \n log10(Expression, counts)", )+
  geom_smooth(method = "lm", se = F, color = "black", linewidth = 1/.pt)+
  geom_text(
    data    = correlations,
    mapping = aes(x = x, y = y, label = label),
    parse=TRUE, size=5/.pt)+
  mytheme + theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(), 
                  #legend.position = "right", 
                  legend.key.size = unit(2, 'mm')) + guides(color=guide_legend(nrow = 2, byrow=T))
#+   stat_cor(method = "pearson", label.x = 10, label.y = 30)
C1


cor.test.result <- cor.test(df$IELs, 
                            df$ENSG00000132170 , method = "pearson")

correlations <- data.frame(correlation_coefficient = cor.test.result$estimate, p.value = cor.test.result$p.value)

correlations <- correlations %>% mutate_all(as.numeric) %>%
  mutate(label = ifelse(`p.value` > 0.01, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `=", f_num(`p.value`,2), "`")), `p.value`))%>% 
  mutate(label = ifelse(`p.value` < 0.01, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `=", f_num(`p.value`,3), "`")), label)) %>%
  mutate(label = ifelse(`p.value` < 0.001, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `", "< .001", "`")),label))%>% 
  mutate(p_full = ifelse(`p.value` > 0.01, round(p.value,2), `p.value`))%>% 
  mutate(p_full = ifelse(`p.value` < 0.01, round(p.value,3), p_full)) %>%
  mutate(p_full = ifelse(`p.value` < 0.001, formatC(p.value, format = "e", digits = 2),p_full))

correlations$x <- 20
correlations$y <- 60

C2 <- ggscatter(as.data.frame(df), x = "ENSG00000132170", y = "IELs",color = "label2",palette = palette_bio, size = 0.5) + #c(PGCd = "#6500ff", PGCp = "#ff6500")
  labs(x = "PPARG mRNA \n log10(Expression, counts)", y= "IEL density", )+
  geom_smooth(method = "lm", se = F, color = "black", linewidth = 1/.pt)+
  geom_text(
    data    = correlations,
    mapping = aes(x = x, y = y, label = label),
    parse=TRUE, size=5/.pt)+
  mytheme+theme(#legend.position = "right", 
    legend.key.size = unit(2, 'mm'))  + guides(color=guide_legend(nrow = 2, byrow=T))
C2

C <- ggarrange(C1,C2, ncol = 1, align = "v", common.legend = T, legend = "top")
C
#-----------------------------------------------------------------------panel D NOS2 expression --------------------------------------------------------------------------

D <- ggboxplot(as.data.frame(melt[melt$Var1 == "ENSG00000007171",]), x = "label2", y = "log10_Expr", outlier.colour = NA, order = c("GFDd", "GFDp", "PGCd","PGCp"),#c("GFD", "PGCd","PGCp"), #unique(metadata$label3), #c("GFDd", "GFDp", "PGCd","PGCp")
               ylab = "NOS2 mRNA \n log10(Expression, counts)",# yscale="log2",
               palette = palette_bio, # c("#ffd700","#ffad00", "#ff934d","#ff6500"),
               fill="label2" #, #facet.by = c("Name"),ncol = 3,scales = "free_y",
               # title = "TGM2 mRNA expression"
)+
  #facet_grid(Marker~Name, scales="free", space="free_x")+
  geom_jitter(size=0.5, alpha=0.5)+
  #scale_y_continuous(limits=c(0.7, 3.0))+
  geom_signif(data=anno[anno$external_gene_name == "NOS2",],
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3) + mytheme +
  theme(axis.title.x = element_blank(),
        axis.text.y=element_blank(), legend.position='none')

D

#-----------------------------------------------------------------------panel E  GSZ---------------------------------------------------------

groupGSZ2 <- groupGSZ[groupGSZ$gene.sets %in% c("BIOCARTA_NO2IL12_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                                                 "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION"),]
groupGSZ2$gene.sets <- factor(groupGSZ2$gene.sets, levels = c("BIOCARTA_NO2IL12_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                                                              "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION"))
groupGSZ2 <- groupGSZ2[order(groupGSZ2$gene.sets),]
groupGSZ2$y_pos <- c(-1.2, -0.5, 0.2,
                     6,8,10,
                     32, 36, 40,
                     9,11,13
                     )



df4 <- melt(df3)
df4 <- df4[df4$variable %in% c("BIOCARTA_NO2IL12_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                               "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION"),]
colnames(df4)[4] <- "gene.sets"
df4$gene.sets <- factor(df4$gene.sets, levels = c("BIOCARTA_NO2IL12_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY",
                                                  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION"))





groupGSZ2$gene.sets2 <- gsub("_", " ",groupGSZ2$gene.sets )
df4$gene.sets2 <- gsub("_", " ",df4$gene.sets )

E <- ggboxplot(df4, x = "label3", y = "value", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
               ylab =  "GSZ",
               palette = palette_3, fill="label3" #, facet.by = 'gene.sets', ncol =  1,scales = "free_y"
) +
  facet_wrap(~ gene.sets2, 
             labeller = label_wrap_gen(multi_line = TRUE), ncol = 1, scales = "free_y")+
  geom_jitter(size=0.5, alpha=0.5)+
  scale_x_discrete(labels=c("GFDd+p","PGCd","PGCp"))+
  geom_signif(data=groupGSZ2,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position= y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)+ mytheme +
  theme(legend.position = "none", axis.title.x = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

E



#-----------------------------------------------------------------------panel Fheatmap--------------------------------------------------------

markers <- readxl::read_xlsx(path = "./input_data/41467_2022_32691_MOESM4_ESM.xlsx", sheet = 1) #data from Atlasy et al, 2022
markers_T <- readxl::read_xlsx(path = "./input_data/41467_2022_32691_MOESM4_ESM.xlsx", sheet = 4) 
markers_CD4 <- readxl::read_xlsx(path = "./input_data/41467_2022_32691_MOESM4_ESM.xlsx", sheet = 6)
markers_MF <- readxl::read_xlsx(path = "./input_data/41467_2022_32691_MOESM4_ESM.xlsx", sheet = 2)



clusters <- data.frame(cluster = c(0,1,2,3,4),
                       Cell.Type = c("PCs","T cells" ,"MCs","mf/DC", "B cells"))
markers <- merge(markers, clusters)
clusters_T <- data.frame(cluster = c(0,1,2,3,4),
                         Cell.Type = c("CD4+ T-cells","CD8+ T-cells" ,"NK","ILC", "proliferating cells"))
markers_T <- merge(markers_T, clusters_T)
clusters_CD4 <- data.frame(cluster = c(0,1),
                           Cell.Type = c("C0 CD4+ T-cells","C1 CD4+ T-cells"))
markers_CD4 <- merge(markers_CD4, clusters_CD4)
clusters_MF <- data.frame(cluster = c(0,1,2,3),
                          Cell.Type = c("MFs","Inf-MF" ,"cDC2","cDC1"))
markers_MF <- merge(markers_MF, clusters_MF)

markers <- rbind(markers,markers_T,markers_CD4, markers_MF)

#select top up and down-reg in each group
top_genes <- data.frame()
for(i in unique(markers$Cell.Type)){
  df <- markers[markers$Cell.Type == i,]
  df <- df %>% arrange(p_val_adj) 
  df1 <- df[df$avg_logFC > 0, ]
  genes1 <- df1[c(1:20), c("Cell.Type","gene")]
  df2 <- df[df$avg_logFC < 0, ]
  genes2 <- df2[c(1:20), c("Cell.Type","gene")]
  top_genes <- rbind(top_genes, genes1, genes2)
}

#searchAttributes(mart = human, pattern = "type")
annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                     filters = "external_gene_name", values = top_genes$gene, mart = human)
top_genes <- merge(top_genes, annot_human, by.x = "gene", by.y = "external_gene_name")

top_genes_M <- top_genes[top_genes$ensembl_gene_id %in% rownames(normalized_counts) &
                         top_genes$Cell.Type %in% c("CD4+ T-cells","CD8+ T-cells","PCs","NK", "Inf-MF" ),]

openxlsx::write.xlsx(top_genes_M, file = "./output_data/Immune_Marker_genes_Figure_4.xlsx")

#top_genes_M <- top_genes_M[1:50,]
counts1<- normalized_counts[top_genes_M$ensembl_gene_id,]

#rename rows
rownames(counts1) <- top_genes_M$gene[match(rownames(counts1), top_genes_M$ensembl_gene_id)]

#define the order of rows if needed
counts1<-counts1[match(rownames(counts1),top_genes_M$gene),]

#scale the data
m.data_sc<-t(scale(t(counts1)))

#sort columns
#metadata <- colData(CD)
metadata <- metadata[order(metadata$label3, metadata$IELs),]
m.data_sc <- m.data_sc[,metadata$ID]

#create annotation 1(Sample type)
type1<- metadata$label2
type3 <- metadata$label
type2 <- top_genes_M$Cell.Type
n = length(unique(top_genes_M$Cell.Type))
palette.cell <- hcl.colors(palette = "Dark2", n = n)
#palette.cell <- c("#72E2AD", "#B7F1B2")

#create annotation 1(Sample type)
type1<- metadata$label2
type3 <- metadata$label

top = HeatmapAnnotation( `VH:CrD` = anno_points(metadata$VHCrD, size = unit(0.5, "mm"), axis_param = list(side = "left", gp = gpar(fontsize = 3))),
                         #annotation_name_gp = gpar(fontsize = 3),
                         `IELs` = anno_points(metadata$IELs, size = unit(0.5, "mm"), axis_param = list(side = "left", gp = gpar(fontsize = 3))),
                         annotation_name_gp = gpar(fontsize = 5),
                         Samples = anno_block(gp = gpar(fill = palette_bio, lty="blank"),
                                              labels = c("GFDd","GFDp", "PGCd","PGCp"), 
                                              labels_gp = gpar(col = "white", fontsize = 5, face = "bold")),
                         gap = unit(1, "mm"),
                         height = unit(1, "cm"))
#height = unit(0.3, "cm"))

row2 = rowAnnotation(foo = anno_block(gp = gpar(fill = palette.cell, lty="blank"), labels = c("CD4+ T-cells","CD8+ T-cells","PCs","NK", "Inf-MF" ), 
                                      labels_gp = gpar(col = "white", fontsize = 5, face = "bold")),
                     width = unit(0.3, "cm"))

col_fun = colorRamp2(c(min(m.data_sc), quantile(m.data_sc, 0.25), 0, quantile(m.data_sc, 0.75),  max(m.data_sc)), 
                     c("#B2182B","#F4A582"  ,"#F7F7F7","#92C5DE", "#2166AC"))

#col_fun = colorRamp2(c(min(m.data_sc),  0,   max(m.data_sc)), 
#                     c("#B2182B" ,"#F7F7F7", "#2166AC"))

ht_list = Heatmap(m.data_sc, 
                  name = "Row Z-Score of Expression",
                  top_annotation = top, #ha1,
                  show_column_names = FALSE, 
                  col=col_fun,
                  column_title_gp = gpar(fontsize = 1, col="White"),
                  column_order= metadata$ID,
                  row_names_side = "right",
                  show_row_names = FALSE,
                  #row_order = order,
                  cluster_rows = T, 
                  show_row_dend = F,
                  row_names_gp = gpar(fontsize = 5), 
                  column_split = factor(type1, levels = c("GFDd","GFDp", "PGCd","PGCp")),
                  column_gap = unit(0.5,"mm"),
                  row_split = factor(type2, levels=c("CD4+ T-cells","CD8+ T-cells","PCs","NK", "Inf-MF" )),
                  row_title = NULL,
                  #row_title_gp = gpar(fontsize = 6),
                  row_gap = unit(0.5,"mm"),
                  cluster_row_slices = FALSE,
                  heatmap_legend_param = list(
                    at = c(ceiling(min(m.data_sc)):floor(max(m.data_sc))),
                    labels = c(ceiling(min(m.data_sc)):floor(max(m.data_sc))),
                    title = "Row Z-Score",
                    labels_gp = gpar(fontsize = 4),
                    title_gp = gpar(fontsize = 4), direction = "vertical",
                    legend_width=unit(2,"cm"),grid_height = unit(0.2, "cm"),
                    title_position = "topcenter"),
                  left_annotation = row2)

Fh <- grid.grabExpr(draw(ht_list, merge_legend = F, heatmap_legend_side = "right", auto_adjust = FALSE))
as.ggplot(Fh)

#----------------------------------------------------------------------------combine all together---------------------------------------------------------------------------------
p1 <- ggarrange(A,B, ncol=1, align = "v")
p2 <- ggarrange(p1, C, ncol = 2)
p3 <- ggarrange(D, E, ncol = 1, heights = 1:5)

p <- arrangeGrob(p2,p3, Fh,
                 ncol = 3, nrow = 4, layout_matrix = rbind(c(1,1,2),c(1,1,2),c(3,3,2),c(3,3,2)))

p <- as_ggplot(p) +                        # transform to a ggplot
  draw_plot_label(label = c("a", "b", "c","d", "e", "f"), size = 11,
                  x = c(0, 0,   0.33,0.66, 0.66, 0 ), 
                  y = c(1, 0.75, 1,1, 0.66, 0.5))
p



graph2ppt(file="./output_data/Figure_4.pptx", width=7.08, height=7.28)

