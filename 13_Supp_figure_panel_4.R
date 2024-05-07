library(readxl)
library(DESeq2)
library(VennDiagram)
library(reshape2)
library(rlist)
library(numform)
library(ggpubr)
library(mGSZ)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(plotrix)
library(plyr)
library(dplyr)
library(export)
#open data sets for comparisons
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




#--------------------------------------------------------------------------panel A ------------------------------------------------------------------------------------------------------------------------------

#retrieve Sample GSZ

SampleGSZ <- CDres[["SampleGSZ"]]
GSZ_melt <- melt(SampleGSZ)
colnames(GSZ_melt)[2] <- "gene.sets"

j="REACTOME_INTERFERON_GAMMA_SIGNALING"
GSZ_melt_panelA <- as.data.frame(GSZ_melt[GSZ_melt$gene.sets == j,])
GSZ_melt_panelA <- as.data.frame(merge(GSZ_melt_panelA,metadata, by.x="Samples", by.y = "ID"))
groupGSZ <- CDres[["groupGSZ"]]

df <-   as.data.frame(GSZ_melt_panelA) %>% group_by(gene.sets, label3) %>% 
  dplyr::summarise(max=max(value))
groupGSZ_panelB <- merge(groupGSZ,df, by.x=c("gene.sets", "group1"), by.y=c("gene.sets", "label3"))


bp<-ggboxplot(GSZ_melt_panelA, x = "label3", y = "value", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
              palette = c("#00b347","#6500ff","#ff6500"),
              fill="label3")+
  labs(y = paste0(j," \n gene set GSZ"))+
  scale_x_discrete(labels=c("GFDd+p","PGCd","PGCp"))+
  geom_jitter(size=0.5, alpha=0.5)+
  geom_hline(yintercept=c(mean(GSZ_melt_panelA$value)), linetype="dashed", color = "grey60")+
  mytheme+
  geom_signif(data=groupGSZ_panelB,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=max*c(1.2,1.3,1.4)),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)

A <- bp + mytheme+
  theme( legend.position = "none", axis.title.x = element_blank())
A

#--------------------------------------------------------------------------panel B ------------------------------------------------------------------------------------------------------------------------------

#retrieve Sample GSZ

j="REACTOME_INTERLEUKIN_21_SIGNALING"
GSZ_melt_panelB <- as.data.frame(GSZ_melt[GSZ_melt$gene.sets == j,])
GSZ_melt_panelB <- as.data.frame(merge(GSZ_melt_panelB,metadata, by.x="Samples", by.y = "ID"))
groupGSZ <- CDres[["groupGSZ"]]

df <-   as.data.frame(GSZ_melt_panelB) %>% group_by(gene.sets, label3) %>% 
  dplyr::summarise(max=max(value))
groupGSZ_panelB <- merge(groupGSZ,df, by.x=c("gene.sets", "group1"), by.y=c("gene.sets", "label3"))


bp<-ggboxplot(GSZ_melt_panelB, x = "label3", y = "value", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
              palette = c("#00b347","#6500ff","#ff6500"),
              fill="label3")+
  labs(y = paste0(j," \n gene set GSZ"))+
  scale_x_discrete(labels=c("GFDd+p","PGCd","PGCp"))+
  geom_jitter(size=0.5, alpha=0.5)+
  geom_hline(yintercept=c(mean(GSZ_melt_panelB$value)), linetype="dashed", color = "grey60")+
  mytheme+
  geom_signif(data=groupGSZ_panelB,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=max*c(1.2,1.3,1.4)),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)

B <- bp + mytheme+
  theme( legend.position = "none", axis.title.x = element_blank())
B

#---------------------------------------------------------------------------panelC-------------------------------------------------------------------------------------------------
#open all files in folder
files_names <- list.files(path="./input_data/Ct_files", pattern=NULL, all.files=FALSE, full.names=FALSE)

CT.data <- list()
for(i in 1:length(files_names)){
  Ctdata<-readxl::read_xlsx(path = paste0("./input_data/Ct_files/", files_names[i]), skip = 18)
  CT.data[[i]] <- Ctdata
}

#open metadata
meta <- readxl::read_xlsx("./input_data/ZED1227-IFNg dose-responce in organoids.xlsx", sheet = "samples meta")

#check tech replicates
#function for removing outliers
goodCt <- function(data, treshhold) {
  mean <-ddply(data, .(Target, Sample), summarize,
               `mean.Ct` = mean(na.omit(as.numeric(as.character(Cq)))),`sd.Ct` = sd(na.omit(Cq)))
  badCt_genes <- na.omit(mean[mean$sd.Ct >treshhold,])
  if (nrow(badCt_genes)!= 0){
    badCt <- data.frame()
    for(i in 1:dim(badCt_genes)[1]){
      tmp <- Ctdata[Ctdata$Sample == badCt_genes$Sample[i],]
      tmp2 <- tmp[tmp$Target == badCt_genes$Target[i],]
      a <-cbind(1,2, abs(tmp2$Cq[1]-tmp2$Cq[2]))
      b <-  cbind(1,3, abs(tmp2$Cq[1]-tmp2$Cq[3]))
      c <-  cbind(2,3, abs(tmp2$Cq[2]-tmp2$Cq[3]))
      tmp3 <- rbind(a,b,c)
      mdiff <- subset(tmp3, tmp3[,3] == min(tmp3[,3]))
      tmp4 <- tmp2[c(mdiff[,1],mdiff[,2]),] 
      diff <- anti_join(tmp2,tmp4)
      badCt <- rbind(badCt,diff)
    }
    
    Ctdata_goodCt <- anti_join(data,badCt)
  } else {
    Ctdata_goodCt <- data
  }
  
  return(Ctdata_goodCt)
}
#add metadata and define reference
reference.sample <- "0uM ZED1227 + 0 U/mL hIFNg"
reference.gene <-"GAPDH"

#ddCt calculations
ddCt_comb <- data.frame()
for(i in 1:length(CT.data)){
  Ctdata <- CT.data[[i]]
  Ctdata <- goodCt(data=Ctdata, treshhold = 0.25)
  Ctdata <- merge(Ctdata, meta, by = "Sample")
  #calculate means for all samples
  mean <-plyr::ddply(Ctdata, .(Target, Sample,`Bio rep`,Treatment), summarize,
               `mean.Ct` = mean(na.omit(as.numeric(as.character(Cq)))),
               `sd.Ct` = sd(na.omit(Cq)))
  #calculate dCt + '2^(-dCt)'
  dCtdata<-data.frame() ##calculation of dCt
  for (n in as.character(unique(mean$Sample))){
    tmp3<-mean[mean$Sample==n,]
    tmp3$dCt<-tmp3$mean.Ct - tmp3$mean.Ct[tmp3$Target == reference.gene]
    tmp3$dCt.error<-sqrt((tmp3$sd.Ct)^2+
                           (tmp3$sd.Ct[tmp3$Target == reference.gene])^2)
    tmp3$'2^(-dCt)' <- 2^(-tmp3$dCt)
    tmp3$conf_int1<-2^-(tmp3$dCt+tmp3$dCt.error)
    tmp3$conf_int2<-2^-(tmp3$dCt-tmp3$dCt.error)
    dCtdata<-rbind(dCtdata,tmp3)
  }
  ## calculate ddCt
  ddCt_results<-data.frame()
  for(i in unique(mean$`Bio rep`)){
    tmp5 <- dCtdata[dCtdata$`Bio rep`==i,]
    for (e in as.character(unique(dCtdata$Target))){
      tmp4<-tmp5[tmp5$Target==e,]
      #tmp4$p.value <- t.test(`2^(-dCt)` ~ medium,data=tmp4,var.equal=TRUE)$p.value
      if (reference.sample %in%tmp4$Treatment ==T){
        tmp4$ddCt<- tmp4$dCt - tmp4$dCt[tmp4$Treatment == reference.sample]
        tmp4$ddCt.error<-tmp4$dCt.error
        tmp4$rate<-2^-(tmp4$ddCt)
        tmp4$conf_int3<-2^-(tmp4$ddCt+tmp4$ddCt.error)
        tmp4$conf_int4<-2^-(tmp4$ddCt-tmp4$ddCt.error)
        ddCt_results<-rbind(ddCt_results,tmp4)
      } else {
        tmp4$ddCt<-NA
        tmp4$ddCt.error<-NA
        tmp4$rate<-NA
        tmp4$conf_int3<-NA
        tmp4$conf_int4<-NA
        ddCt_results<-rbind(ddCt_results,tmp4)
        
      }
    }} 
  ddCt_results %>% 
    mutate(FC = ifelse(rate < 1, round(-1/rate,2), round(rate,2)))%>% 
    mutate(conf_int5 = ifelse(rate < 1, -1/conf_int3, conf_int3))%>% 
    mutate(conf_int6 = ifelse(rate < 1, -1/conf_int4, conf_int4))->ddCt_results
  ddCt_comb <- rbind(ddCt_comb, ddCt_results)
  
}

#writexl::write_xlsx(ddCt_comb, path = "ddCt_results.xlsx")

#plotting
ddCt_comb <- merge(ddCt_comb, meta[, c("Sample","Zed1227_uM","IFNg_U")], by = "Sample")

plotdata <- ddCt_comb[ddCt_comb$Target=="TGM2",] 
plotdata$`Bio rep` <- factor(plotdata$`Bio rep`,levels = c("hDUO04", "hDUO02", "hDUO05"))
plotdata$rate <- round(plotdata$rate,1)


P <- ggbarplot(plotdata, x = "Bio rep", y ='FC', fill="Treatment",
               label = TRUE,  lab.pos = "out", legend = "none",
               #palette = c("#b6afd6", "#6e60ae", "#4d4379")
)+
  facet_grid(`Zed1227_uM` ~ `IFNg_U`,labeller = label_both)+
  geom_hline(yintercept = 1, colour="red", linetype="dashed")+
  geom_errorbar(aes(ymin = conf_int5, ymax = conf_int6), width = 0.5)+
  xlab(NULL) + ylab('TGM2 Fold Change vs mock treatment') + xlab("Biological replicate")

P <- ggpar(P,xtickslab.rt = 45,main = "INFg treatment for 24h")

#P


zed.labs <- c("0\u03bcM ZED1227", "5\u03bcM ZED1227","25\u03bcM ZED1227","50\u03bcM ZED1227")
names(zed.labs) <- c(0,5,25,50)
#plotdata <- 
  
  plotdata  %>% group_by(IFNg_U, Zed1227_uM) %>% select(IFNg_U, Zed1227_uM, rate, Sample)%>% 
    mutate(avg=mean(rate,na.rm=T), min=avg-std.error(rate,na.rm=T),max=avg+std.error(rate,na.rm=T)) -> df
  
Cpl <-  ggplot(unique(df), aes(x=as.factor(IFNg_U), y=rate, fill = factor(IFNg_U))) + 
    geom_jitter(size=0.5) + 
    geom_errorbar(aes(ymin = min, ymax = max), width = 0.5,size=0.2)+
    geom_errorbar(aes(ymin = min, ymax = avg), width = 0.5,size=0.2) + 
    facet_wrap(. ~ Zed1227_uM, ncol = 2, labeller = labeller( Zed1227_uM = zed.labs))+
    ylab("TGM2 Relative Expression") +
    xlab("IFNg, U/mL") +
    scale_y_sqrt() + mytheme +
  theme(legend.position = "none")
Cpl   



#=============================================================================panel D ===============================================
df_act <- read_xlsx(path="./input_data/Caco2_colorimetric_assay_results.xlsx", sheet = 1)
df_act_test <- df_act[df_act$Type == "test", c("Sample", "Name","nmol")]
df_act_test_mean <- df_act_test %>% group_by(Sample, Name) %>% summarise(mean_nmol = mean(nmol))

df_prot <- read_xlsx(path="./input_data/Caco2_colorimetric_assay_results.xlsx", sheet = 2)
df_prot_test <- df_prot[df_prot$Type == "test", c("Content","Protein conc, ug/mL")]

df_prot_test_mean <- df_prot_test %>% group_by(Content) %>% summarise(mean_prot_conc = mean(`Protein conc, ug/mL`)) 

df_comb <- merge(df_act_test_mean, df_prot_test_mean, by.x = "Name", by.y = "Content")
df_comb$prot_well_ug <- df_comb$mean_prot_conc*50/1000
df_comb$activity <- (df_comb$mean_nmol*1.5)/(df_comb$prot_well_ug/1000)

palette = c(M="#FFBD35",`I+Z`="#502064",I="#8267BE",Z="#3FA796")

D <- ggbarplot(df_comb, x = "Sample", y = "activity", color = "Sample", palette = palette,
               add = c("mean_se"))+
  geom_jitter(size=0.5, alpha=0.5)+
  scale_x_discrete(limits = c("M","Z","I","I+Z"))+
  scale_y_continuous(breaks = seq(0, 800, by = 200), limits=c(0, 950))+
  labs(y = "Transglutaminase activity, \n nmole hydroxamate/120min/mg)")+ mytheme+
  theme(axis.title.x = element_blank())

D

#----------------------------------------------------------------------------combine all together---------------------------------------------------------------------------------
fplot <- arrangeGrob(A, B, Cpl, D,
                     ncol = 2, nrow = 2, 
                     layout_matrix = rbind(c(1,2),
                                           c(3,4)))
p <- as_ggplot(fplot) +                        # transform to a ggplot
  draw_plot_label(label = c("a", "b", "c","d"), size = 11,
                  x = c(0, 0.5,0,0.5), 
                  y = c(1, 1, 0.5, 0.5))
p
#setwd("D:/OneDrive/OneDrive - TUNI.fi/ISE/RNA seq CD new project/our data")

graph2ppt(file="./output_data/Figure_S4.pptx", width=7.08, height=7.28)
