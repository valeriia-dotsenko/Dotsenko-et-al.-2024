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
library(export)
library(dplyr)


#open data sets for comparisons
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
CDraw <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")

metadata <- colData(CDraw)

ORG <- readRDS('./input_data/hDuo ZED1277 DEG results.rds')

DEG1 <- CDres[[1]]
DEG_IFNg <- ORG[[2]]

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
#----------------------------------------------------------------------------panel A Venn ------------------------------------------------------------------------------------------------
p1 <- grobTree(venn.diagram(
  x = list(DEG1[DEG1$Gene.type=="DE",]$ID,
           DEG_IFNg[DEG_IFNg$Gene.type=="DE",]$ID),
  category.names = c("PGCp VS GFDp", "IFNg VS mock"),
  filename = NULL,
  fontfamily="sans", #counts font
  cex = 0.6, #counts size; 0.9 is 10pt
  cat.fontfamily =  "sans", #labels font
  cat.cex = 0.5, #labels size; 0.8 is appr 9pt
  #cat.default.pos = "outer", #label position
  cat.pos = c(200, 20),
  cat.dist = 0.01,
  col=c(alpha("#ff6500",0.5), alpha("#8267BE",0.5)),
  fill = c(alpha("#ff6500",0.3), alpha("#8267BE",0.3)),
  rotation.degree = 45,
  margin = 0.12))

A <- ggarrange(as.ggplot(p1))
A  


#--------------------------------------------------------------------------panel B ------------------------------------------------------------------------------------------------------------------------------

SampleGSZ <- CDres[["SampleGSZ"]]
groupGSZ <- CDres[["groupGSZ"]]


GSZ_melt <- melt(SampleGSZ)
colnames(GSZ_melt)[2] <- "gene.sets"
GSZ_melt <- merge(GSZ_melt,metadata, by.x="Samples", by.y = "ID")
GSZ_melt_panelB <- as.data.frame(GSZ_melt[GSZ_melt$gene.sets == "epithelial_IFNg_responce_CeD",])
GSZ_melt_panelB$value_norm <- GSZ_melt_panelB$value-mean(GSZ_melt_panelB$value)

writexl::write_xlsx(GSZ_melt_panelB[, c("Samples","gene.sets","value"  )], path = "./input_data/epithelial_IFNg_responce_CeD_GSZ.xlsx")

groupGSZ <- groupGSZ %>% filter(gene.sets == "epithelial_IFNg_responce_CeD")
groupGSZ %>%  #p-full for figure caption
  mutate(p_full = ifelse(`pvalue` > 0.01, round(pvalue,2), `pvalue`))%>% 
  mutate(p_full = ifelse(`pvalue` < 0.01, round(pvalue,3), p_full)) %>%
  mutate(p_full = ifelse(`pvalue` < 0.001, formatC(pvalue, format = "e", digits = 2),p_full))->groupGSZ

df <-   as.data.frame(GSZ_melt_panelB) %>% group_by(gene.sets, label3) %>% 
  dplyr::summarise(max=max(value_norm))
groupGSZ_panelB <- merge(groupGSZ,df, by.x=c("gene.sets", "group1"), by.y=c("gene.sets", "label3"))

j="Epithelial responce to IFNg"
bp<-ggboxplot(GSZ_melt_panelB, x = "label3", y = "value_norm", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
              palette = c("#00b347","#6500ff","#ff6500"),
              fill="label3")+
  labs(y = paste0(j," \n gene set GSZ"))+
  scale_x_discrete(labels=c("GFDd+p","PGCd","PGCp"))+
  geom_jitter(size=0.5, alpha=0.5)+
  geom_hline(yintercept=c(mean(GSZ_melt_panelB$value_norm)), linetype="dashed", color = "grey60")+
  mytheme+
  geom_signif(data=groupGSZ_panelB,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=max*c(1.7)),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)

B <- bp + mytheme+
  theme( legend.position = "none", axis.title.x = element_blank())
B

#-----------------------------------------------------------------------panel C---------------------------------------------------------------------
# barplot
#j=names(mglist)[6]
#GSZ_melt1 <- as.data.frame(GSZ_melt[GSZ_melt$gene.sets==j,])
GSZ_melt_panelB <- GSZ_melt_panelB[order(GSZ_melt_panelB$label2, GSZ_melt_panelB$VHCrD),]

plot <- ggplot(GSZ_melt_panelB, aes(x = Samples, y = value_norm,
                              fill = label2))+
  geom_hline(yintercept=c(sd(GSZ_melt_panelB$value_norm), -sd(GSZ_melt_panelB$value_norm)), linetype="dashed", color = "grey")+
  geom_bar(stat="identity")+
  scale_fill_manual(values=palette_bio)+
  geom_point(aes(x=Samples, y=VHCrD,group=1),stat="identity",color="black",size=0.3)+
  #geom_line(aes(x=Samples, y=VHCrD,group=1),stat="identity",color="red",size=0.5)+
  labs(x="Sample",
       y= paste0(j," \n gene set GSZ"))+
  rremove("x.text")+
  scale_x_discrete(limits=GSZ_melt_panelB$Samples)+
  
  scale_y_continuous(sec.axis=sec_axis(~.*1,name="VH:CrD"))+
  mytheme + theme(axis.text.x = element_blank())
plot
#onoff plot  
GSZ_melt_panelB$x <- seq(0, 1, length.out = nrow(GSZ_melt_panelB))
GSZ_melt_panelB$y <- rep(0, nrow(GSZ_melt_panelB))
GSZ_melt_panelB %>% 
  mutate(onoff = ifelse(GSZ_melt_panelB$value_norm > sd(GSZ_melt_panelB$value_norm), 1, 0)) %>% 
  mutate(onoff = ifelse(GSZ_melt_panelB$value_norm < -sd(GSZ_melt_panelB$value_norm), -1, onoff))%>% 
  mutate(col = ifelse(onoff == 1, "#98a148", '#ffffcc')) %>% 
  mutate(col = ifelse(onoff == -1, "#e74c3c", col))%>% 
  mutate(status = ifelse(onoff == 1, "on", '-')) %>% 
  mutate(status = ifelse(onoff == -1, "off", status)) -> GSZ_melt_panelB

tmp <- GSZ_melt_panelB[, c("Samples","gene.sets","value_norm", "status")]
#write.table(tmp, file=paste("./", j, " Sample GSZ.txt",sep=""))

onoff <- ggplot(GSZ_melt_panelB, aes(x, y, fill= as.factor(status))) + 
  geom_tile()+
  #scale_fill_manual(labels = unique(GSZ_melt_panelB$status), values = unique(col))+
  scale_fill_manual(values=c("on"="#98a148", "off"="#e74c3c", "-"='#ffffcc'))+
  geom_vline(xintercept=c(0.5,0.8), linetype="dashed", color = "black")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.background= element_rect(fill = NA),
        legend.title=element_blank(), legend.position = "bottom",
        legend.text = element_text( size = 6, color = "black"),
        legend.key.size = unit(3, 'mm'))

onoff

C1 <- ggarrange(plot, onoff, heights = c(2, 1),
                ncol = 1, nrow = 2, align = "v")
C1


#-------------------------------------------------------------------------panel D TG2 barplot #1------------------------------------------------------------------------------
IDs <-c('ENSG00000198959') # TGM2
#prepare data
normalized_counts <- counts(CDraw, normalized=TRUE)

melt<-melt(normalized_counts[rownames(normalized_counts) %in% IDs,])
melt$variable <- rownames(melt)
melt$log10_Expr <- log10(melt$value+1)
melt <-merge(melt, metadata[,c("ID","label2")],by.x="variable", by.y="ID") 
group_comp <- data.frame(`Sample Group` = c("PGCp","PGCd","PGCp"),
                         `Control Group` = c("GFDp","GFDd","PGCd"),
                         `Group Comparisons` = c("PGCp VS GFDp","PGCd VS GFDd","PGCp VS PGCd"))

anno<- data.frame()
for(i in 1:nrow(group_comp)){
  DEG1 <- CDres[[i]]
  anno1 <-subset(DEG1, DEG1$ID %in% IDs, select=c(external_gene_name,baseMean,padj)) 
  anno1$group1 <- group_comp$`Sample.Group`[i]
  anno1$group2 <- group_comp$`Control.Group`[i]
  anno <- rbind(anno,anno1)
}

anno %>% 
  mutate(p_new = ifelse(padj > 0.01, c(paste("italic('P')~`=", f_num(padj,2), "`")), padj))%>% 
  mutate(p_new = ifelse(padj < 0.01, c(paste("italic('P')~`=", f_num(padj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(padj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->anno

anno$y_pos <- max(melt$log10_Expr)*seq(from=1.1, to=1.5, length.out = nrow(group_comp))
#anno[anno$group2 %in% c("GFDp","GFDd"),]$group2 <- "GFD"


bp<-ggboxplot(as.data.frame(melt), x = "label2", y = "log10_Expr", outlier.colour = NA, order = c("GFDd", "GFDp", "PGCd","PGCp"),#c("GFD", "PGCd","PGCp"), #unique(metadata$label3), #c("GFDd", "GFDp", "PGCd","PGCp")
              ylab = "TGM2 mRNA \n log10(Expression, counts)",# yscale="log2",
              palette = palette_bio, # c("#ffd700","#ffad00", "#ff934d","#ff6500"),
              fill="label2" #, #facet.by = c("Name"),ncol = 3,scales = "free_y",
              # title = "TGM2 mRNA expression"
)+
  #facet_grid(Marker~Name, scales="free", space="free_x")+
  geom_jitter(size=0.5, alpha=0.5)+
  #scale_y_continuous(limits=c(0.7, 3.0))+
  theme(axis.title.x = element_blank(),
        axis.text.y=element_blank())+
  labs(fill = element_blank())+
  theme(strip.text.x = element_text(size = 6))+
  geom_signif(data=anno,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)


D <- bp + mytheme + theme(axis.title.x = element_blank(), legend.position = "none")


D

#---------------------------------------------------------------------------panel F----------------------------------------------------------------------------
Orgraw <- readRDS( "./input_data/hDuo ZED1277 DESeqDataSet object.rds")
normalized_counts <- counts(Orgraw, normalized=TRUE)
metadata2 <- colData(Orgraw)

melt2<-melt(normalized_counts[rownames(normalized_counts) %in% IDs,])
melt2$variable <- rownames(melt2)
melt2$log10_Expr <- log10(melt2$value+1)
melt2 <-merge(melt2, metadata2[,c("Sample_Name","label")],by.x="variable", by.y="Sample_Name") 
#group_comp <- data.frame(`Sample Group` = c("Z","M","M","I"),
#                         `Control Group` = c("M","I","I+Z","I+Z"),
#                         `Group Comparisons` = c("Z VS M","I VS M","M VS I+Z", "I VS I+Z"))
group_comp <- data.frame(`Sample Group` = c("M"),`Control Group` = c("I"),`Group Comparisons` = c('I_VS_M'))
anno<- data.frame()
for(i in 1:nrow(group_comp)){
  DEG1 <- ORG[[group_comp$Group.Comparisons[i]]]
  anno1 <-subset(DEG1, DEG1$ID %in% IDs, select=c(external_gene_name,baseMean,padj)) 
  anno1$group1 <- group_comp$`Sample.Group`[i]
  anno1$group2 <- group_comp$`Control.Group`[i]
  anno <- rbind(anno,anno1)
}

anno %>% 
  mutate(p_new = ifelse(padj > 0.01, c(paste("italic('P')~`=", f_num(padj,2), "`")), padj))%>% 
  mutate(p_new = ifelse(padj < 0.01, c(paste("italic('P')~`=", f_num(padj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(padj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->anno

anno %>% 
  mutate(p_full = ifelse(`padj` > 0.01, round(padj,2), `padj`))%>% 
  mutate(p_full = ifelse(`padj` < 0.01, round(padj,3), p_full)) %>%
  mutate(p_full = ifelse(`padj` < 0.001, formatC(padj, format = "e", digits = 2),p_full))->anno

anno$y_pos <- max(melt2$log10_Expr)*seq(from=1.05, to=1.2, length.out = nrow(group_comp))
#palette = c("#FFBD35","#3FA796","#8267BE","#502064")
palette2 = c("#FFBD35","#8267BE")

bp2<-ggboxplot(as.data.frame(melt2[melt2$label %in% c("I","M"),]), x = "label", y = "log10_Expr", outlier.colour = NA, order = c("M","I"),#c("M", "Z","I","I+Z"), #unique(metadata$label3), #c("GFDd", "GFDp", "PGCd","PGCp")
              ylab = "TGM2 mRNA \n log10(Expression, counts)",# yscale="log2",
              palette = palette2,
              fill="label" #, #facet.by = c("Name"),ncol = 3,scales = "free_y",
              # title = "TGM2 mRNA expression"
)+
  #facet_grid(Marker~Name, scales="free", space="free_x")+
  geom_jitter(size=0.5, alpha=0.5)+
  #scale_y_continuous(limits=c(0.7, 3.0))+
  theme(axis.title.x = element_blank(),
        axis.text.y=element_blank())+
  labs(fill = element_blank())+
  theme(strip.text.x = element_text(size = 6))+
  geom_signif(data=anno,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)


E <- bp2 + mytheme + theme(axis.title.x = element_blank(), legend.position = "none")

E

#-----------------------------------------------------------------------panel E--------------------------------------------------------------------------------
GSZ_melt_panelB <- merge(GSZ_melt_panelB, as.data.frame(melt)[, -2], by.x = c("Samples", "label2"), by.y = c("variable", "label2"))
GSZ_melt_panelB_challenge <- GSZ_melt_panelB[GSZ_melt_panelB$Timepoint == "challenge",]

cor.test.result <- cor.test(GSZ_melt_panelB_challenge$log10_Expr, 
                            GSZ_melt_panelB_challenge$value_norm , method = "pearson")

correlations <- data.frame(correlation_coefficient = cor.test.result$estimate, p.value = cor.test.result$p.value)

correlations <- correlations %>% mutate_all(as.numeric) %>%
  mutate(label = ifelse(`p.value` > 0.01, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `=", f_num(`p.value`,2), "`")), `p.value`))%>% 
  mutate(label = ifelse(`p.value` < 0.01, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `=", f_num(`p.value`,3), "`")), label)) %>%
  mutate(label = ifelse(`p.value` < 0.001, c(paste("R==", round(correlation_coefficient, 2), "~';'~italic(P)~ `", "< .001", "`")),label)) %>% 
  mutate(p_full = ifelse(`p.value` > 0.01, round(p.value,2), `p.value`))%>% 
  mutate(p_full = ifelse(`p.value` < 0.01, round(p.value,3), p_full)) %>%
  mutate(p_full = ifelse(`p.value` < 0.001, formatC(p.value, format = "e", digits = 2),p_full))
correlations$x <- -1
correlations$y <- 1.7

Fp <- ggscatter(GSZ_melt_panelB_challenge, x = "value_norm", y = "log10_Expr",color = "label2", palette = c(PGCd = "#6500ff", PGCp = "#ff6500"), size = 0.5) + 
  labs(x = paste0(j, " \n gene set GSZ"), y= "TGM2 mRNA \n log10(Expression, counts)", )+
  #grids(linetype = "dashed")+
  geom_text(
    data    = correlations,
    mapping = aes(x = x, y = y, label = label),
    parse=TRUE, size=5/.pt)+
  geom_smooth(method="lm",formula=y~x, se = F, color = "black", linewidth = 1/.pt)+
  #stat_cor(method = "pearson", label.x = 1.7, label.y = 0.7)+
  mytheme
Fp

#----------------------------------------------------------------------------combine all together ver2 ---------------------------------------------------------------------------------


p1 <- ggarrange(A,B, ncol=2, align = "h", widths = c(1,2))
p1
p2 <- ggarrange(p1,C1, ncol=1, heights = c(2,1))
p2

p3 <- ggarrange(D,E,Fp, ncol=1, align = "v")
p3
p <- ggarrange(p2,p3, ncol=2, widths = c(2.5,1))
p

p <- arrangeGrob(p, 
                 ncol = 1, nrow = 1) 
#as_ggplot(p)
p <- as_ggplot(p) +                        # transform to a ggplot
  draw_plot_label(label = c("a", "b", "c","d", "e","f"), size = 11,
                  x = c(0, 0.22, 0,   0.75, 0.75, 0.75), 
                  y = c(1, 1,    0.33, 1,    0.66,0.33))
p
graph2ppt(file="./output_data/Figure_3.pptx", width=7.08, height=7.28)



