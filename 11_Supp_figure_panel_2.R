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


library(ggplot2)
library(enrichR) #library(devtools) #install_github("wjawaid/enrichR")
library(stringr)
library(ggpubr)
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


DEG1 <- CDres[[1]]
DEG2 <- CDres[[2]]
DEG3 <- CDres[[3]]

#---------------------------------------------------------------------------panel D---------------------------------------------------------------------------------------------
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE




genes <- DEG1[DEG1$Gene.type=="DE",]$external_gene_name

genes1 <- DEG3[DEG3$Gene.type=="DE",]$external_gene_name 

dbs_selected <- c("TRRUST_Transcription_Factors_2019")#"ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X") #"TRRUST_Transcription_Factors_2019") "ENCODE_TF_ChIP-seq_2015"
enriched <- enrichr(genes, dbs_selected)
df <- enriched[[1]]
#df[df$P.value< 0.05,]

TF_toshow <- c("STAT1","RELA","IRF1","HNF4A","CREB5","PPARG", "SREBF1")

#str_split(df$Term, " ", simplify = T)
df <- cbind(df, str_split(df$Term, " ", simplify = T))
df <- df[df$`2`=="human",]
df

df_selected <- df[df$`1` %in% TF_toshow,]
plotting_df <- df_selected[, c("1","P.value", "Overlap")]
plotting_df <- cbind(plotting_df, Freq = as.numeric(str_split(plotting_df$Overlap, "/", simplify = T)[,1]))
plotting_df$`-log10(p-value)` <- -log10(plotting_df$P.value)


p <- plotting_df %>% 
  ggplot(aes(x = reorder(`1`,`-log10(p-value)`), y = Freq, fill=c("#86C5D8"))) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts")+ scale_fill_manual(values= c("#86C5D8"))+
  geom_text(aes(label=Freq,y=Freq/2, x=`1`), vjust = 0.5, size=6/.pt)
p

f <- sum(abs(range(na.omit(plotting_df$Freq))))/(diff(range(plotting_df$`-log10(p-value)`))+1)*0.5 #factor for top y-axis (coordinates are flipped) scaling
plotting_df$y <- f*(plotting_df$`-log10(p-value)`-median(plotting_df$`-log10(p-value)`))
plotting_df$group <- "1"
#highlifgt non-significant points
highlight_df <- plotting_df %>% filter(P.value > 0.05)
#add scatter plot on top of barplot


bp <-  p + 
  geom_line(data=plotting_df, aes(x=reorder(`1`,`-log10(p-value)`),y=`-log10(p-value)`, group = group), color="#09e8b5",size = 0.5) +
  geom_point(data=plotting_df, aes(x=reorder(`1`,`-log10(p-value)`), y=`-log10(p-value)`, group = group), shape=21,color="black", fill="#09e8b5", size=3) +
  theme_bw()+ 
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "-log10(p-value)")) +
  ggtitle("PGCp VS GFD")

bp <- bp + mytheme + theme(axis.title.x.top = element_text(color='#4daf4a',size=7), 
                           axis.text.x.top = element_text(color='#4daf4a',size=6), 
                           axis.title.y = element_blank(),#text(size=16),
                           strip.text.x = element_text(size = 6),
                           legend.position = "none")
bp



enriched1 <- enrichr(genes1, dbs_selected)
df1 <- enriched1[[1]]

#str_split(df1$Term, " ", simplify = T)
df1 <- cbind(df1, str_split(df1$Term, " ", simplify = T))
df1 <- df1[df1$`2`=="human",]
head(df1)

df1_selected <- df1[df1$`1` %in% TF_toshow,]
plotting_df1 <- df1_selected[, c("1","P.value", "Overlap")]
plotting_df1 <- cbind(plotting_df1, Freq = as.numeric(str_split(plotting_df1$Overlap, "/", simplify = T)[,1]))
plotting_df1$`-log10(p-value)` <- -log10(plotting_df1$P.value)
plotting_df <- plotting_df[order(plotting_df$`-log10(p-value)`, decreasing = F),]
plotting_df1$`1` <- factor(plotting_df1$`1`, levels=unique(plotting_df$`1`))

p <- plotting_df1 %>% 
  ggplot(aes(x = `1`, y = Freq, fill=c("#86C5D8"))) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts") + scale_fill_manual(values= c("#86C5D8"))+
  geom_text(aes(label=Freq,y=Freq/2, x=`1`), vjust = 0.5, size=6/.pt)+
  scale_x_discrete(position = "top")
p

f <- sum(abs(range(na.omit(plotting_df1$Freq))))/(diff(range(plotting_df1$`-log10(p-value)`))+1)*0.5 #factor for top y-axis (coordinates are flipped) scaling
plotting_df1$y <- f*(plotting_df1$`-log10(p-value)`-median(plotting_df1$`-log10(p-value)`))
plotting_df1$group <- "1"
#highlifgt non-significant points
highlight_df1 <- plotting_df1[plotting_df1$P.value > 0.05,]
#add scatter plot on top of barplot


bp2 <-  p + 
  geom_line(data=plotting_df1, aes(x=`1`,y=`-log10(p-value)`, group = group), color="#09e8b5",size = 0.5) +
  geom_point(data=plotting_df1, aes(x=`1`, y=`-log10(p-value)`, group = group), shape=21,color="black", fill="#09e8b5", size=3) +
  geom_point(data=highlight_df1, aes(x=`1`,y=`-log10(p-value)`, group=group), color='#8c8f8e', size=3) +
  theme_bw()+ 
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "-log10(p-value)")) +
  ggtitle("PGCp VS PGCd")

bp2 <- bp2 + mytheme + theme(axis.title.x.top = element_text(color='#4daf4a',size=7), 
                             axis.text.x.top = element_text(color='#4daf4a',size=6), 
                             axis.title.y = element_blank(),#text(size=16),
                             strip.text.x = element_text(size = 6),
                             legend.position = "none")
bp2
D <- ggarrange(bp, bp2, ncol = 2)
D

graph2ppt(file="./output_data/Figure_S2.pptx", width=4.76, height=5.11)


ENRICHR_res <- list(df,df1)
names(ENRICHR_res) <- c("PGCp_VS_GFDp_ENRICHR_results", "PGCp_VS_PGCd_ENRICHR_results")
openxlsx::write.xlsx(ENRICHR_res, file = "./output_data/ALL_ENRICHR_results_CeD_biopsies.xlsx")
