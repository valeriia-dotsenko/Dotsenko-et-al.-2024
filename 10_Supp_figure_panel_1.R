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



#prepare data
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
CD <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")
metadata <- as.data.frame(colData(CD))

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


#--------------------------------------------------------------------------panel A ( Violin)-------------------------------------------------------------------------------------------------
DEG1 <- CDres[[1]]
DEG2 <- CDres[[2]]
DEG3 <- CDres[[3]]

common2 <- Reduce(function(...) merge(..., all = TRUE, by = "external_gene_name"),
                  list(DEG1[DEG1$Gene.type == "DE",c("external_gene_name","log2FoldChange")],
                       DEG2[DEG2$Gene.type == "DE",c("external_gene_name","log2FoldChange")],
                       DEG3[DEG3$Gene.type == "DE",c("external_gene_name","log2FoldChange")]))
colnames(common2) <- c("Gene","PGCp_VS_GFD","PGCd_VS_GFD","PGCp_VS_PGCd")
common2 <- unique(as.data.frame(common2))
rownames(common2) <- common2$Gene
common2[is.na(common2)] <-  "non-DE"
common2[common2 > 0 & common2 != "non-DE"] <-  "DE UP"
common2[common2 < 0 & common2 != "non-DE"] <-  "DE DOWN"
common2$Gene <- rownames(common2)

#DESeq2 size factor normalization
normalized_counts <- counts(CD, normalized=TRUE)

genes56_DOWN <- DEG1[DEG1$external_gene_name %in% common2[common2$PGCp_VS_GFD %in% c( "DE DOWN") & common2$PGCp_VS_PGCd %in% c( "DE DOWN"),]$Gene,]$ID
counts56_DOWN <- normalized_counts[genes56_DOWN,]
#change id to names
tmp <- DEG1[DEG1$ID %in% rownames(counts56_DOWN), c("ID","external_gene_name")]

rownames(counts56_DOWN) <- tmp$external_gene_name[match(rownames(counts56_DOWN), tmp$ID)]

#melt data for plotting
data1_melt <- melt(counts56_DOWN) 
data1_melt <- merge(data1_melt, metadata, by.x = "Var2",by.y = "ID")
library(plyr)
data1_mean <- ddply(data1_melt, .(Var1,label3), plyr::summarize,  `mean(value)`=mean(value))
data1_mean$`log10_value` <- log10(data1_mean$`mean(value)`)
data1_mean$`log10_value`[data1_mean$log10_value=="-Inf"] <- (-0.01)

#Kruskal-Wallis test followed by the Dunn's post hoc test
kr.test <- kruskal.test(log10_value~label3, data = data1_mean)$p.value
data1_mean$label3 = factor(data1_mean$label3,levels=c("GFD",  "PGCd", "PGCp"))
l <- DunnTest(log10_value~label3, data = data1_mean,method="fdr") #c("holm", "hochberg", "hommel", "bonferroni", "BH","BY", "fdr", "none"))
anno_df4 <- data.frame(Reduce(rbind, l[1]))[1:length(levels(data1_mean$label3)),]
anno_df4$comparison <- row.names(anno_df4)
anno_df4 <- cbind(anno_df4,as.data.frame(matrix(unlist(strsplit(anno_df4$comparison,"-")), ncol=2)))
colnames(anno_df4) <- c("mean.rank.diff","p.adj","comparison","group1","group2")
anno_df4 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new)) %>%
  mutate(p_full = ifelse(`p.adj` > 0.01, round(p.adj,2), `p.adj`))%>% 
  mutate(p_full = ifelse(`p.adj` < 0.01, round(p.adj,3), p_full)) %>%
  mutate(p_full = ifelse(`p.adj` < 0.001, formatC(p.adj, format = "e", digits = 2),p_full))->anno_df4
anno_df4$y_pos <- max(data1_mean$log10_value)*seq(from=1.05, to=1.2, length.out=nrow(anno_df4))

E1 <- ggplot(data1_mean, aes(x = factor(label3), y = log10_value, fill=label3)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.shape=NA)+
  scale_fill_manual(values=alpha(c("#00b347","#6500ff","#ff6500"),0.8))+
  geom_jitter(size=0.5, fill="black", alpha=0.5)+
  scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(-0.1, 4.5))+
  annotate(geom="text", x=2, y=-0.1, label="downgulated genes",size=6/.pt)+
  #annotate("text", x=2, y=4.5, label = "Kruskal-Wallis test",size=6/.pt)+
  #annotate("text", x=0, y=4.4, label = paste("italic(P) == ", f_num(kr.test,2)), parse = TRUE,size=6/.pt)+
  ylab("log10(Expression, counts) ")+
  geom_signif(annotations = anno_df4$p_new,
              y_position = anno_df4$y_pos, xmin=anno_df4$group1, xmax=anno_df4$group2,textsize=6/.pt, manual= F, parse=T, size=0.3, tip_length = 0)


E1 <- E1 + mytheme + theme(axis.title.x = element_blank())

genes56_UP <- DEG1[DEG1$external_gene_name %in% common2[common2$PGCp_VS_GFD %in% c( "DE UP") & common2$PGCp_VS_PGCd %in% c( "DE UP"),]$Gene,]$ID
counts56_UP <- normalized_counts[genes56_UP,]
#change id to names
tmp <- DEG1[DEG1$ID %in% rownames(counts56_UP), c("ID","external_gene_name")]
rownames(counts56_UP) <- tmp$external_gene_name[match(rownames(counts56_UP), tmp$ID)]
#melt data for plotting
data1_melt <- melt(counts56_UP) 
data1_melt <- merge(data1_melt, metadata, by.x = "Var2",by.y = "ID")
data1_mean <- ddply(data1_melt, .(Var1,label3), plyr::summarize,  `mean(value)`=mean(value))
data1_mean$`log10_value` <- log10(data1_mean$`mean(value)`)
data1_mean$`log10_value`[data1_mean$log10_value=="-Inf"] <- (-0.01)

kr.test <- kruskal.test(log10_value~label3, data = data1_mean)$p.value
data1_mean$label3 = factor(data1_mean$label3,levels=c("GFD",  "PGCd", "PGCp"))
l <- DunnTest(log10_value~label3, data = data1_mean,method="fdr") #c("holm", "hochberg", "hommel", "bonferroni", "BH","BY", "fdr", "none"))
anno_df5 <- data.frame(Reduce(rbind, l[1]))[1:length(levels(data1_mean$label3)),]
anno_df5$comparison <- row.names(anno_df5)
anno_df5 <- cbind(anno_df5,as.data.frame(matrix(unlist(strsplit(anno_df5$comparison,"-")), ncol=2)))
colnames(anno_df5) <- c("mean.rank.diff","p.adj","comparison","group1","group2")
anno_df5 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))%>%
  mutate(p_full = ifelse(`p.adj` > 0.01, round(p.adj,2), `p.adj`))%>% 
  mutate(p_full = ifelse(`p.adj` < 0.01, round(p.adj,3), p_full)) %>%
  mutate(p_full = ifelse(`p.adj` < 0.001, formatC(p.adj, format = "e", digits = 2),p_full))->anno_df5
#anno_df4 <- merge(anno_df4, data1_mean %>%
#group_by(label3) %>%
#  dplyr::summarise(mean = mean(log10_value)), by.x= "group1", by.y = "label3")
anno_df5$y_pos <- max(data1_mean$log10_value)*seq(from=1.05, to=1.2, length.out=nrow(anno_df4))

#anno_df4 %>% rstatix::add_y_position(data1_mean, formula = label3 ~ cancer_group, fun="max")

E2 <- ggplot(data1_mean, aes(x = factor(label3), y = log10_value, fill=label3)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.shape=NA)+
  scale_fill_manual(values=alpha(c("#00b347","#6500ff","#ff6500"),0.8))+
  geom_jitter(size=0.5, fill="black", alpha=0.5)+
  scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(-0.1, 4.5))+
  annotate(geom="text", x=2, y=-0.1, label="upregulated genes",size=6/.pt)+
  #annotate("text", x=2, y=4.5, label = "Kruskal-Wallis test",size=6/.pt)+
  #annotate("text", x=0.6, y=4.4, label = paste("italic(P) == ", f_num(kr.test,2)), parse = TRUE,size=6/.pt)+
  geom_signif(annotations = anno_df5$p_new,y_position = anno_df5$y_pos, 
              xmin=anno_df5$group1, xmax=anno_df5$group2, textsize=6/.pt, manual= F, parse=T, size=0.3, tip_length = 0)


E2 <- E2 + mytheme + theme(axis.title.y = element_blank(),axis.title.x = element_blank())

E <- ggarrange(E1,E2,ncol = 2, nrow = 1, #font.label = list(size = 16),
                   common.legend = TRUE, legend="none")
E

#--------------------------------------------------------------------------panel F ( Violin)-------------------------------------------------------------------------------------------------
genes124_DOWN <- DEG3[DEG3$external_gene_name %in% common2[common2$PGCp_VS_GFD %in% c("non-DE") & common2$PGCp_VS_PGCd %in% c("DE DOWN"),]$Gene,]$ID
counts124_DOWN <- normalized_counts[genes124_DOWN,]
#change id to names
tmp <- DEG1[DEG1$ID %in% rownames(counts124_DOWN), c("ID","external_gene_name")]
rownames(counts124_DOWN) <- tmp$external_gene_name[match(rownames(counts124_DOWN), tmp$ID)]
#melt data for plotting
data1_melt <- melt(counts124_DOWN) 
data1_melt <- merge(data1_melt, metadata, by.x = "Var2",by.y = "ID")
library(plyr)
data1_mean <- ddply(data1_melt, .(Var1,label3), plyr::summarize,  `mean(value)`=mean(value))
data1_mean$`log10_value` <- log10(data1_mean$`mean(value)`)
data1_mean$`log10_value`[data1_mean$log10_value=="-Inf"] <- (-0.01)

#Kruskal-Wallis test followed by the Dunn's post hoc test
kr.test <- kruskal.test(log10_value~label3, data = data1_mean)$p.value
data1_mean$label3 = factor(data1_mean$label3,levels=c("GFD",  "PGCd", "PGCp"))
l <- DunnTest(log10_value~label3, data = data1_mean,method="fdr") #c("holm", "hochberg", "hommel", "bonferroni", "BH","BY", "fdr", "none"))
anno_df4 <- data.frame(Reduce(rbind, l[1]))[1:length(levels(data1_mean$label3)),]
anno_df4$comparison <- row.names(anno_df4)
anno_df4 <- cbind(anno_df4,as.data.frame(matrix(unlist(strsplit(anno_df4$comparison,"-")), ncol=2)))
colnames(anno_df4) <- c("mean.rank.diff","p.adj","comparison","group1","group2")
anno_df4 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))%>%
  mutate(p_full = ifelse(`p.adj` > 0.01, round(p.adj,2), `p.adj`))%>% 
  mutate(p_full = ifelse(`p.adj` < 0.01, round(p.adj,3), p_full)) %>%
  mutate(p_full = ifelse(`p.adj` < 0.001, formatC(p.adj, format = "e", digits = 2),p_full))->anno_df4
anno_df4$y_pos <- max(data1_mean$log10_value)*seq(from=1.05, to=1.2, length.out=nrow(anno_df4))

F1 <- ggplot(data1_mean, aes(x = factor(label3), y = log10_value, fill=label3)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.shape=NA)+
  scale_fill_manual(values=alpha(c("#00b347","#6500ff","#ff6500"),0.8))+
  geom_jitter(size=0.5, fill="black", alpha=0.5)+
  scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(-0.1, 5.0))+
  annotate(geom="text", x=2, y=-0.1, label="downgulated genes",size=6/.pt)+
  #annotate("text", x=2, y=4.8, label = "Kruskal-Wallis test",size=6/.pt)+
  #annotate("text", x=2, y=4.7, label = paste("italic(P) == ", f_num(kr.test,2)), parse = TRUE,size=6/.pt)+
  ylab("log10(Expression, counts) ")+
  geom_signif(annotations = anno_df4$p_new,
              y_position = anno_df4$y_pos, xmin=anno_df4$group1, xmax=anno_df4$group2,textsize=6/.pt, manual= F, parse=T, size=0.3, tip_length = 0)


F1 <- F1 + mytheme + theme(axis.title.x = element_blank())

genes124_UP <- DEG3[DEG3$external_gene_name %in% common2[common2$PGCp_VS_GFD %in% c("non-DE") & common2$PGCp_VS_PGCd %in% c("DE UP"),]$Gene,]$ID
counts124_UP <- normalized_counts[genes124_UP,]
#change id to names
tmp <- DEG1[DEG1$ID %in% rownames(counts124_UP), c("ID","external_gene_name")]
rownames(counts124_UP) <- tmp$external_gene_name[match(rownames(counts124_UP), tmp$ID)]
#melt data for plotting
data1_melt <- melt(counts124_UP) 
data1_melt <- merge(data1_melt, metadata, by.x = "Var2",by.y = "ID")
data1_mean <- ddply(data1_melt, .(Var1,label3), plyr::summarize,  `mean(value)`=mean(value))
data1_mean$`log10_value` <- log10(data1_mean$`mean(value)`)
data1_mean$`log10_value`[data1_mean$log10_value=="-Inf"] <- (-0.01)

kr.test <- kruskal.test(log10_value~label3, data = data1_mean)$p.value
data1_mean$label3 = factor(data1_mean$label3,levels=c("GFD",  "PGCd", "PGCp"))
l <- DunnTest(log10_value~label3, data = data1_mean,method="fdr") #c("holm", "hochberg", "hommel", "bonferroni", "BH","BY", "fdr", "none"))
anno_df4 <- data.frame(Reduce(rbind, l[1]))[1:length(levels(data1_mean$label3)),]
anno_df4$comparison <- row.names(anno_df4)
anno_df4 <- cbind(anno_df4,as.data.frame(matrix(unlist(strsplit(anno_df4$comparison,"-")), ncol=2)))
colnames(anno_df4) <- c("mean.rank.diff","p.adj","comparison","group1","group2")
anno_df4 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))%>%
  mutate(p_full = ifelse(`p.adj` > 0.01, round(p.adj,2), `p.adj`))%>% 
  mutate(p_full = ifelse(`p.adj` < 0.01, round(p.adj,3), p_full)) %>%
  mutate(p_full = ifelse(`p.adj` < 0.001, formatC(p.adj, format = "e", digits = 2),p_full))->anno_df4
#anno_df4 <- merge(anno_df4, data1_mean %>%
#group_by(label3) %>%
#  dplyr::summarise(mean = mean(log10_value)), by.x= "group1", by.y = "label3")
anno_df4$y_pos <- max(data1_mean$log10_value)*seq(from=1.05, to=1.2, length.out=nrow(anno_df4))

#anno_df4 %>% rstatix::add_y_position(data1_mean, formula = label3 ~ cancer_group, fun="max")

F2 <- ggplot(data1_mean, aes(x = factor(label3), y = log10_value, fill=label3)) +
  geom_violin() +
  geom_boxplot(width=0.1,outlier.shape=NA)+
  scale_fill_manual(values=alpha(c("#00b347","#6500ff","#ff6500"),0.8))+
  geom_jitter(size=0.5, fill="black", alpha=0.5)+
  scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(-0.1, 5))+
  annotate(geom="text", x=2, y=-0.1, label="upregulated genes",size=6/.pt)+
  #annotate("text", x=2, y=4.8, label = "Kruskal-Wallis test",size=6/.pt)+
  #annotate("text", x=2, y=4.7, label = paste("italic(P) == ", f_num(kr.test,2)), parse = TRUE,size=6/.pt)+
  geom_signif(annotations = anno_df4$p_new,y_position = anno_df4$y_pos, 
              xmin=anno_df4$group1, xmax=anno_df4$group2, textsize=6/.pt, manual= F, parse=T, size=0.3, tip_length = 0)


F2 <- F2 + mytheme + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
Fpl <- ggarrange(F1,F2,ncol = 2, nrow = 1, #font.label = list(size = 16),
               common.legend = TRUE, legend="none")
Fpl



p1 <- ggarrange(E,Fpl,
                ncol = 2, common.legend = T, align = "h")
p1


p <- p1 +    
  draw_plot_label(label = c("a", "b"), size = 11,
                                 x = c(0, 0.5), 
                               y = c(1, 1)) # transform to a ggplot
p

graph2ppt(file="./output_data/Figure_S1.pptx", width=4.76, height=5.11)
