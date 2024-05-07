library(rlist)
library(DESeq2)
library(biomaRt)
library(reshape2)
library(mGSZ)
library(ggpubr)
library(dplyr)
library(numform)
library(rstatix)
library(readxl)
library(export)
library(cowplot)

#open data sets for comparisons
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
CDraw <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")

#metadata <- colData(CDraw)


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
                                                     strip.text.x = element_text(size = 5, colour = "black"))
palette3 <- c("#404040","#a0b6f7","#f2f261", "#404040","#a0b6f7","#f2f261")

#-----------------------------------------------------------------------------panelA-------------------------------------------------------------------------------------------------------------
metadata <- read_xlsx(path = "./input_data/metadata-comb.xlsx")
metadata_wide <- reshape2::dcast(as.data.frame(metadata[,c("newID", "Timepoint", "VH", "HLA_Genotype_Group", "Treatment" )]), newID + HLA_Genotype_Group +Treatment ~ Timepoint, value.var = "VH")

metadata_wide$baseline <- as.numeric(metadata_wide$baseline)
metadata_wide$challenge <- as.numeric(metadata_wide$challenge)
metadata_wide <- metadata_wide[metadata_wide$HLA_Genotype_Group != "Not identified", ]
metadata_wide <- na.omit(metadata_wide)


metadata_wide1 <- metadata_wide[metadata_wide$Treatment == "drug",]

res.aov2 <- metadata_wide1 %>% anova_test(challenge ~ baseline+HLA_Genotype_Group)
get_anova_table(res.aov2)
get_test_label(res.aov2,  detailed = TRUE)



gg <- ggplot(metadata_wide1, aes(x = baseline,  y = challenge,color = HLA_Genotype_Group)) +
  geom_point(aes(color=HLA_Genotype_Group), size=0.5) +
  geom_smooth(method = "lm",  mapping = aes(y = predict(lm(challenge ~ baseline + HLA_Genotype_Group,data = metadata_wide1), metadata_wide1)), linewidth=0.5) +
  scale_color_manual(values = palette3) +
  ylab("VH in PGCd group") +
  xlab("VH in GFDd group") +  mytheme +
  theme(legend.position = "right") +
  guides(color=guide_legend(override.aes=list(fill=NA))) #removes legend grey background
#gg

# Pairwise comparisons
pwc2 <- metadata_wide1 %>% 
  emmeans_test(  #Performs pairwise comparisons between groups using the estimated marginal means
    challenge ~ HLA_Genotype_Group, covariate = baseline,
    p.adjust.method = "bonferroni", detailed = T
  )
#pwc2 %>% filter(p.adj < 0.05)
pwc2

pwc2 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->pwc2
pwc2$contrast <- paste0(pwc2$group1, " - ", pwc2$group2)


gg_effect <- ggplot(data = pwc2, aes(y = contrast,estimate)) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),width = 0, color = "black") + # confidence level of effect
  geom_point(size = 0.5) + # estimate of effect
  geom_vline(xintercept = 0, linetype = 2) + # draw a line at effect = 0
  annotate(geom = "text", label = pwc2$p_new, size = 5/.pt, manual=TRUE, parse=T, y = 1:3, x = round(max(pwc2$conf.high)*4,2)) + # p-value. The y coordinates are set by eye
  ylab("Contrast") +
  coord_cartesian(xlim = c(round(min(pwc2$conf.low)*1.05,2),round(max(pwc2$conf.high)*5,2))) +
  scale_x_continuous(position="top") +
  theme_pubr() + mytheme + theme(axis.title.y = element_blank())

gg_effect

A <- plot_grid(gg_effect,
                gg ,
                nrow=2,
                align = "v",
                axis = "rl",
                rel_heights = c(0.4, 1)
)
A

#-----------------------------------------------------------------------------panelB-------------------------------------------------------------------------------------------------------------
metadata_wide <- reshape2::dcast(as.data.frame(metadata[,c("newID", "Timepoint", "CrD", "HLA_Genotype_Group", "Treatment" )]), newID + HLA_Genotype_Group +Treatment ~ Timepoint, value.var = "CrD")

metadata_wide$baseline <- as.numeric(metadata_wide$baseline)
metadata_wide$challenge <- as.numeric(metadata_wide$challenge)
metadata_wide <- metadata_wide[metadata_wide$HLA_Genotype_Group != "Not identified", ]
metadata_wide <- na.omit(metadata_wide)


metadata_wide1 <- metadata_wide[metadata_wide$Treatment == "drug",]

res.aov2 <- metadata_wide1 %>% anova_test(challenge ~ baseline+HLA_Genotype_Group)
get_anova_table(res.aov2)
get_test_label(res.aov2,  detailed = TRUE)



gg <- ggplot(metadata_wide1, aes(x = baseline,  y = challenge,color = HLA_Genotype_Group)) +
  geom_point(aes(color=HLA_Genotype_Group), size=0.5) +
  geom_smooth(method = "lm",  mapping = aes(y = predict(lm(challenge ~ baseline + HLA_Genotype_Group,data = metadata_wide1), metadata_wide1)), linewidth=0.5) +
  scale_color_manual(values = palette3) +
  ylab("CrD in PGCd group") +
  xlab("CrD in GFDd group") +  mytheme +
  theme(legend.position = "right") +
  guides(color=guide_legend(override.aes=list(fill=NA))) #removes legend grey background
#gg

# Pairwise comparisons
pwc2 <- metadata_wide1 %>% 
  emmeans_test(  #Performs pairwise comparisons between groups using the estimated marginal means
    challenge ~ HLA_Genotype_Group, covariate = baseline,
    p.adjust.method = "bonferroni", detailed = T
  )
#pwc2 %>% filter(p.adj < 0.05)
pwc2

pwc2 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->pwc2
pwc2$contrast <- paste0(pwc2$group1, " - ", pwc2$group2)


gg_effect <- ggplot(data = pwc2, aes(y = contrast,estimate)) +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),width = 0, color = "black") + # confidence level of effect
  geom_point(size = 0.5) + # estimate of effect
  geom_vline(xintercept = 0, linetype = 2) + # draw a line at effect = 0
  annotate(geom = "text", label = pwc2$p_new, size = 5/.pt, manual=TRUE, parse=T, y = 1:3, x = round(max(pwc2$conf.high)*4,2)) + # p-value. The y coordinates are set by eye
  ylab("Contrast") +
  coord_cartesian(xlim = c(round(min(pwc2$conf.low)*1.05,2),round(max(pwc2$conf.high)*5,2))) +
  scale_x_continuous(position="top") +
  theme_pubr() + mytheme + theme(axis.title.y = element_blank())

gg_effect

B <- plot_grid(gg_effect,
                gg ,
                nrow=2,
                align = "v",
                axis = "rl",
                rel_heights = c(0.4, 1)
)
B
#---------------------------------------------------------------------------panel C--------------------------------------------------------------------------------------------------

#retrieve Sample GSZ
gene.sets <- c("Transit amplifying cells", "Mature enterocytes", "Immune cells", "Duodenal transporters", "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA",
               "REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS","REACTOME_METABOLISM_OF_LIPIDS","REACTOME_PHOSPHOLIPID_METABOLISM")

SampleGSZ <- CDres[["SampleGSZ"]]
GSZ_melt <- reshape2::melt(SampleGSZ)
colnames(GSZ_melt)[2] <- "gene.sets"
GSZ_melt <- as.data.frame(GSZ_melt[GSZ_melt$gene.sets %in% gene.sets,])
GSZ_melt$gene.sets <- factor(as.character(GSZ_melt$gene.sets), levels=gene.sets)
GSZ_melt <- as.data.frame(merge(GSZ_melt,metadata, by.x="Samples", by.y = "ID"))

GSZ_melt$gene.sets <- gsub("_", " ",GSZ_melt$gene.sets )

Cp <- ggplot(GSZ_melt[GSZ_melt$Timepoint == "challenge" & GSZ_melt$HLA_Genotype_Group != "Not identified",], aes(x = label3, y = value)) + 
  geom_jitter(aes(color = HLA_Genotype_Group), size = 0.5,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + ylab("GSZ") +
  facet_wrap(. ~ gene.sets, scales="free", ncol = 4, labeller = label_wrap_gen(multi_line = TRUE)) +
  stat_summary(aes(color = HLA_Genotype_Group), fun.data="mean_sdl", fun.args = list(mult=1), size = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values =  palette3)+ 
  mytheme + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Cp

#-----------------------------------------------------------------------------combine------------------------------------------------------------------
p0 <- ggarrange(A,B,  ncol = 2,nrow = 1, align = "h", common.legend = TRUE)
p0
p1 <- ggarrange(p0,Cp,
                nrow = 2,  heights = c(1,2))
p1
p <- p1 +                        
  draw_plot_label(label = c("a", "b","c"), size = 11,
                  x = c(0, 0.5,0), 
                  y = c(1,1,0.66))
p

graph2ppt(file="./output_data/Extended_data_Figure_2.pptx", width=7.08, height=7.28)
