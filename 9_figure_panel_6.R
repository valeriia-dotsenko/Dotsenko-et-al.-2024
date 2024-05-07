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
library(ReactomePA)
library(topGO)
library(org.Hs.eg.db) #annotation package
library(GOxploreR)
library(rlist)
library(matrixStats)
library(mGSZ)
library(forcats)
library(export)


human =  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

CDraw <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")
normalized_counts <- counts(CDraw, normalized=TRUE)
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
                                                     strip.text.x = element_text(size = 5.5, colour = "black"))

palette_bio <- c("#00b347","#006729","#6500ff","#ff6500")


load("./input_data/MODEL1.rda")
model1<-c("ATP8B2", "PLA2R1", "PDIA3", "TM4SF4")
annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype", "entrezgene_id"),
                     filters = "external_gene_name", values = model1, mart = human)
#=====================================================old model1===================================================
model1_id <- annot_human[annot_human$external_gene_name %in% model1,]$ensembl_gene_id

#-------------------------------------------------------------------------model1 prediction--------------------------------------------------
modeldata <- as.data.frame(normalized_counts[model1_id,])
modeldata$id <- rownames(modeldata)
modeldata <- merge(modeldata, annot_human[,c("ensembl_gene_id","external_gene_name")], by.x="id", by.y="ensembl_gene_id")
rownames(modeldata) <- modeldata$external_gene_name
modeldata <- as.data.frame(t(modeldata[,-c(1,117)]))

modeldata$pred <- predict(fit1,modeldata)
modeldata$Sample <- row.names(modeldata)
modeldata <- as.data.frame(merge(modeldata,metadata, by.x="Sample",by.y="ID"))

#-------------------------------------------------------------------------scatter plots-------------------------------------------------------

modeldata_melt <- melt(modeldata[, c("Sample","ATP8B2","PLA2R1","PDIA3","TM4SF4", "VHCrD","label3")], id.vars=c("Sample","VHCrD","label3"))
modeldata_melt[modeldata_melt$label3 == "GFD",]$label3 <- "GFDd+p"

#color condition, regression line
A <- ggscatter(modeldata_melt, x = "VHCrD", y = "value", color = "label3",
          palette=c("#00b347","#6500ff","#ff6500"),size = 1, alpha = 1)+
  facet_wrap(~variable, scales="free_y",ncol=2) +
  labs(y = expression(paste("Gene expression, UMI counts")),x="VH:CrD")+
  stat_smooth(method = "lm") + mytheme
A

fit <- lm(VHCrD ~ pred - 1, data = modeldata)
n <- length(modeldata$VHCrD) # the sample size
RSMD<-(1/(n-1)*sum((modeldata$pred-modeldata$VHCrD)^2))^0.5

text="VH:CrD"
eqn <- sprintf("italic(R)^2 ~ '=' ~ %.2g *',' ~~ RSMD == %.2f",
               summary(fit)$r.squared, RSMD)
B <- ggscatter(modeldata, x = "pred", y = "VHCrD",
          color = "black", shape = 21, size = 1)+ # Points color, shape and size
  labs(y = paste("observed",text,sep=" "
  ), x=paste("predicted",text,sep=" "
  ))+
  geom_abline(intercept = 0, slope = 1,linetype="dashed", linewidth=0.5, color="red")+
  grids(linetype = "dashed")+
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = eqn, parse = TRUE,
    hjust = 1.1, vjust = -.5, size = 6/.pt)+
  geom_smooth(method="lm",formula=y~0+x) + mytheme
B

modeldata <- modeldata[modeldata$newID!="401021",]
my_comparisons <- list( c("PGCp", "GFDp"), c("PGCd", "GFDd"), c("PGCd", "PGCp"))
j="VHCrD"
stat_comp <- data.frame()
for(i in 1:length(my_comparisons)){
  group1 <- my_comparisons[[i]][1]
  group2 <- my_comparisons[[i]][2]
  N1 <- length(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
  N2 <- length(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
  
  if (N1==N2){
    p.value <- t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$p.value
    est <- t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$estimate
    est_CI <- paste(round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$conf.int[1],2), "to", 
                    round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$conf.int[2],2), sep=" ")
    mean1 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    sd1 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    
    mean2 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    sd2 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    
    tmp <- cbind(group1,group2,N1,mean1,sd1,N2,mean2,sd2,p.value, est, est_CI)
    stat_comp <- rbind(stat_comp,tmp)
  } else{
    p.value <- t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$p.value
    est <- -diff(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$estimate)
    est_CI <- paste(round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$conf.int[1],2), "to", 
                    round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$conf.int[2],2), sep=" ")
    mean1 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    sd1 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    
    mean2 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    sd2 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    tmp <- cbind(group1,group2,N1,mean1,sd1,N2,mean2,sd2,p.value, est, est_CI)
    stat_comp <- rbind(stat_comp,tmp)
  }
  
}


stat_comp %>% mutate_at(c(3:10), as.numeric) %>% 
  mutate(p_new = ifelse( p.value > 0.01, c(paste("italic('P')~`=", f_num( p.value,2), "`")),  p.value))%>% 
  mutate(p_new = ifelse( p.value < 0.01, c(paste("italic('P')~`=", f_num( p.value,3), "`")), p_new))%>%
  mutate(p_new = ifelse( p.value < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))%>%
  mutate(p_new1 = ifelse( p.value > 0.01, c(f_num( p.value,2)),  p.value))%>% 
  mutate(p_new1 = ifelse( p.value < 0.01, c(f_num( p.value,3)), p_new1))%>%
  mutate(p_new1 = ifelse( p.value < 0.001, c("<.001"), p_new1))%>%
  mutate(p_full = ifelse(`p.value` > 0.01, round(p.value,2), `p.value`))%>% 
  mutate(p_full = ifelse(`p.value` < 0.01, round(p.value,3), p_full)) %>%
  mutate(p_full = ifelse(`p.value` < 0.001, formatC(p.value, format = "e", digits = 2),p_full))->stat_comp

stat_comp$mean11 <- paste(round(stat_comp$mean1,2),"$\\pm$",round(stat_comp$sd1,2), sep=" ")
stat_comp$mean21 <- paste(round(stat_comp$mean2,2),"$\\pm$",round(stat_comp$sd2,2), sep=" ")
stat_comp$est1 <- paste(round(stat_comp$est,2)," \n (",stat_comp$est_CI,")", sep="")

#3. Plot the boxplot
p1 <- ggplot(modeldata, aes(x=label2, y=VHCrD, group=label2, fill=label2)) + 
  geom_boxplot(outlier.colour=NA)+
  geom_line(aes(x = label2, group = Pair),linetype="longdash",# color="grey",
            alpha=0.5,size=1/.pt)+
  geom_signif(annotations = stat_comp[3,]$p_new,
              y_position = c( 4.5), xmin=stat_comp[3,]$group1, xmax=stat_comp[3,]$group2,textsize=6/.pt, manual= F, parse=T, size=0.3)+
  labs(y = "VH:CrD observed")+
  ggtitle("Histomorphometry \nmeasured VH:CrD") +
  scale_x_discrete(limits = c("GFDd","PGCd","GFDp","PGCp"), labels=c("GFDd","PGCd","GFDp","PGCp"))+
  scale_fill_manual(values=alpha(palette_bio,1))+
  scale_y_continuous(expand = c(.2,0)) 
p1 <- p1 + geom_point(aes(x = label2), size=0.5)+
            mytheme+
              theme(axis.title.x = element_blank(),legend.position="none",
                plot.margin = margin(0, 0, 0.1, 0.4, "cm"))
p1
j="pred"
stat_comp2 <- data.frame()
for(i in 1:length(my_comparisons)){
  group1 <- my_comparisons[[i]][1]
  group2 <- my_comparisons[[i]][2]
  N1 <- length(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
  N2 <- length(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
  
  if (N1==N2){
    p.value <- t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$p.value
    est <- t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$estimate
    est_CI <- paste(round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$conf.int[1],2), "to", 
                    round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=T)$conf.int[2],2), sep=" ")
    mean1 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    sd1 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    
    mean2 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    sd2 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    
    tmp <- cbind(group1,group2,N1,mean1,sd1,N2,mean2,sd2,p.value, est, est_CI)
    stat_comp2 <- rbind(stat_comp2,tmp)
  } else{
    p.value <- t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$p.value
    est <- -diff(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$estimate)
    est_CI <- paste(round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$conf.int[1],2), "to", 
                    round(t.test(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j], modeldata[modeldata$label2==my_comparisons[[i]][2],][,j], paired=F)$conf.int[2],2), sep=" ")
    mean1 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    sd1 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][1],][,j])
    
    mean2 <- mean(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    sd2 <- sd(modeldata[modeldata$label2==my_comparisons[[i]][2],][,j])
    tmp <- cbind(group1,group2,N1,mean1,sd1,N2,mean2,sd2,p.value, est, est_CI)
    stat_comp2 <- rbind(stat_comp2,tmp)
  }
  
}


stat_comp2 %>% mutate_at(c(3:10), as.numeric) %>% 
  mutate(p_new = ifelse( p.value > 0.01, c(paste("italic('P')~`=", f_num( p.value,2), "`")),  p.value))%>% 
  mutate(p_new = ifelse( p.value < 0.01, c(paste("italic('P')~`=", f_num( p.value,3), "`")), p_new))%>%
  mutate(p_new = ifelse( p.value < 0.001, c(paste("italic('P')~`", "<.001", "`")), p_new))%>%
  mutate(p_new1 = ifelse( p.value > 0.01, c(f_num( p.value,2)),  p.value))%>% 
  mutate(p_new1 = ifelse( p.value < 0.01, c(f_num( p.value,3)), p_new1))%>%
  mutate(p_new1 = ifelse( p.value < 0.001, c("<.001"), p_new1))%>%
  mutate(p_full = ifelse(`p.value` > 0.01, round(p.value,2), `p.value`))%>% 
  mutate(p_full = ifelse(`p.value` < 0.01, round(p.value,3), p_full)) %>%
  mutate(p_full = ifelse(`p.value` < 0.001, formatC(p.value, format = "e", digits = 2),p_full)) -> stat_comp2

stat_comp2$mean11 <- paste(round(stat_comp2$mean1,2),"$\\pm$",round(stat_comp2$sd1,2), sep=" ")
stat_comp2$mean21 <- paste(round(stat_comp2$mean2,2),"$\\pm$",round(stat_comp2$sd2,2), sep=" ")
stat_comp2$est1 <- paste(round(stat_comp2$est,2)," \n (",stat_comp2$est_CI,")", sep="")

#3. Plot the boxplot
p2 <- ggplot(modeldata, aes(x=label2, y=pred, group=label2, fill=label2)) + 
  geom_boxplot(outlier.colour=NA)+
  geom_line(aes(x = label2, group = Pair),linetype="longdash", #color="grey",
            alpha=0.5,size=1/.pt)+
  geom_signif(annotations = stat_comp2[3,]$p_new,
              y_position = c(5), xmin=stat_comp2[3,]$group1, xmax=stat_comp2[3,]$group2,textsize=6/.pt, manual= F, parse=T, size=0.3)+
            labs(y = "VH:CrD predicted")+
  ggtitle("Molecular Histomorphometry \nRegression model based on RNA expression") +   
  scale_x_discrete(limits = c("GFDd","PGCd","GFDp","PGCp"), labels=c("GFDd","PGCd","GFDp","PGCp"))+
  scale_fill_manual(values=alpha(palette_bio,0.6))+
  scale_y_continuous(expand = c(.2,0)) 
p2 <- p2 + geom_point(aes(x = label2), size=0.5)+
            mytheme+
                theme(axis.title.x = element_blank(),legend.position="none",
                                  plot.margin = margin(0, 0, 0.1, 0.4, "cm"))
p2

C <- ggarrange(p1,p2,ncol = 2, nrow = 1, 
                   common.legend = TRUE,legend="none")
C

#----------------------------------------------------------------------------combine all together---------------------------------------------------------------------------------
p0 <- ggarrange(A,B,
                ncol = 2,widths = c(1, 0.5))
p0

p1 <- ggarrange(p0,C,
                nrow = 2, common.legend = T, align = "v", heights=c(1,1))
p1


# Add labels to the arranged plots
p <- p1 +                        
  draw_plot_label(label = c("a", "b", "c"), size = 11,
                  x = c(0, 0.66,0), 
                  y = c(1, 1, 0.5))
p
graph2ppt(file="./output_data/Figure_6.pptx", width=7.08, height=7.28)


