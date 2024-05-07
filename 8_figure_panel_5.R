library(ggplot2)
library(reshape2)
library(dplyr)
library(rstatix)
library(DESeq2)
library(emmeans)
library(ggpubr)
library(data.table)
library(cowplot)
library(numform)
library(ggplotify)
library(ggthemes)
library(readxl)
library(lazyWeave)
library(gridExtra)
library(export)



CDraw <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')

metadata <- read_xlsx("./input_data/metadata-comb.xlsx")
IFNg_GSZ <- read_xlsx("./input_data/epithelial_IFNg_responce_CeD_GSZ.xlsx")
colnames(IFNg_GSZ)[3] <- "epithelial_IFNg_responce"

metadata <- merge(metadata, IFNg_GSZ[, c(1,3)], by.x = "ID", by.y = "Samples", all =T)

#create theme for all plots
mytheme <- theme_classic(base_family='sans') + theme(axis.text.x = element_text(size=6, color = "black"),
                                                     axis.text.y= element_text(size=6, color = "black"),
                                                     axis.title.y = element_text(size=7, color = "black"),
                                                     axis.title.x = element_text(size=7, color = "black"),
                                                     legend.text = element_text( size = 5, color = "black"),
                                                     legend.key.size = unit(3, 'mm'),
                                                     legend.position = "top",
                                                     legend.title = element_blank(),
                                                     plot.title = element_text(hjust = 0.5, size = 7,face = "bold"),
                                                     axis.line = element_line(colour = "black", linewidth=1/.pt),
                                                     panel.background = element_rect(fill = 'White'),
                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     strip.text.x = element_text(size = 5.5, colour = "black"),
                                                     axis.title.x.top = element_text(color="black",size=7), 
                                                     axis.text.x.top = element_text(color='black',size=6))

palette_bio <- c("#00b347","#006729","#6500ff","#ff6500")

#------------------------------------------------------------------------------Table 2 ---------------------------------------------------------------------------------------------------
table2_1 <- metadata[metadata$Timepoint == "baseline" & metadata$Treatment == "drug",] %>% 
  group_by(HLA_Genotype_Group, HLA_Genotype, Treatment) %>%
  dplyr::summarize(n = n())

sum <- sum(table2_1$n)

table2_1 <- table2_1 %>%
  mutate(freq = round(n/sum*100, 1))

table2_2 <- metadata[metadata$Timepoint == "baseline" & metadata$Treatment == "placebo",] %>% 
  group_by(HLA_Genotype_Group, HLA_Genotype, Treatment) %>%
  dplyr::summarize(n = n()) 
sum <- sum(table2_2$n)

table2_2 <- table2_2 %>%
  mutate(freq = round(n/sum*100, 1))

table2 <- rbind(table2_1,table2_2)
writexl::write_xlsx(table2, path = "./output_data/table_2_Distribution_HLA-DQ.xlsx")

metadata[metadata$Timepoint == "baseline",] %>% 
  group_by(HLA_Genotype_Group, Treatment) %>%
  dplyr::summarize(n = n())

#------------------------------------------------------------------------------Table 3 ---------------------------------------------------------------------------------------------------

metadata_wide <- reshape2::dcast(as.data.frame(metadata[,c("newID", "Timepoint", "VHCrD", "HLA_Genotype_Group", "Treatment" )]), newID + HLA_Genotype_Group +Treatment ~ Timepoint, value.var = "VHCrD")
metadata_wide <- metadata_wide %>% arrange(newID) %>% filter(HLA_Genotype_Group != "Not identified") %>% mutate_at(c("baseline", "challenge"), as.numeric)

table3 <- metadata_wide %>% 
  group_by(HLA_Genotype_Group, Treatment) %>%
  dplyr::summarize(n = n(),mean_baseline = mean(baseline), sd_baseline = round(sd(baseline),2), mean_challenge = mean(challenge), sd_challenge = round(sd(challenge),2),
                   p.val= round(t.test(challenge,baseline, paired=T)$p.value,3), #alternative hypothesis: true mean difference is not equal to 0, two-tailed
                   est = round(t.test(challenge,baseline, paired=T)$estimate,2), 
                   est_CI = paste(round(t.test(challenge,baseline, paired=T)$conf.int[1],2), "to", 
                                  round(t.test(challenge,baseline, paired=T)$conf.int[2],2), sep=" "))

writexl::write_xlsx(table3, path = "./output_data/table_3_changes_VHCrD.xlsx")


#------------------------------------------------------------------------------repeated measures ANOVA------------------------------------------------------------------------------------
metadata <- metadata %>% convert_as_factor(newID,Pair, Treatment, Timepoint) %>% mutate_at(c('VHCrD'), as.numeric)
#metadata1 <- metadata[metadata$HLA_Genotype_Group != "Not identified" & metadata$Treatment == "drug",]

res.aov <- metadata[metadata$HLA_Genotype_Group != "Not identified",] %>%  group_by(Treatment) %>%
  anova_test(dv = VHCrD, wid = newID,within = c(Timepoint), between = c(HLA_Genotype_Group))

get_anova_table(res.aov)


# Pairwise comparisons (the same as one.way2)
pwc <- metadata[metadata$HLA_Genotype_Group != "Not identified",] %>%
  group_by(HLA_Genotype_Group, Treatment) %>%
  pairwise_t_test(
    VHCrD ~ Timepoint, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc %>% arrange(Treatment)


# Effect of HLA_Genotype_Group at each time point in placebo group
one.way <- metadata[metadata$HLA_Genotype_Group != "Not identified" & metadata$Treatment == "placebo",] %>%
  group_by(Timepoint) %>%
  anova_test(dv = VHCrD, wid = newID, between = HLA_Genotype_Group) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

# Effect of time point at each HLA_Genotype_Group in placebo group
one.way2 <- metadata[metadata$HLA_Genotype_Group != "Not identified" & metadata$Treatment == "placebo",] %>%
  group_by(HLA_Genotype_Group) %>%
  anova_test(dv = VHCrD, wid = newID, within = c(Timepoint)) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#DRUG non-significant two-way interaction
# comparisons for Timepoint variable
metadata[metadata$HLA_Genotype_Group != "Not identified" & metadata$Treatment == "drug",] %>%
  pairwise_t_test(
    VHCrD ~ Timepoint, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )
# comparisons for HLA_Genotype_Group
metadata[metadata$HLA_Genotype_Group != "Not identified" & metadata$Treatment == "drug",] %>%
  group_by(HLA_Genotype_Group) %>%
  pairwise_t_test(
    VHCrD ~ Timepoint, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )

#---------------------------------------------------------------------------panelA------------------------------------------------------------------------
means <- metadata[metadata$HLA_Genotype_Group != "Not identified",]%>%
  group_by(Treatment, HLA_Genotype_Group,Timepoint) %>%
  get_summary_stats(VHCrD, type = "mean_sd")

palette3 <- c("#404040","#a0b6f7","#f2f261")

A <- ggplot(means, aes(x = Timepoint, y = mean, color = HLA_Genotype_Group))+
  geom_point(size=1)+
  geom_line(aes(x = Timepoint, y = mean, group  = HLA_Genotype_Group), linewidth = 0.2)+
  facet_wrap("Treatment")+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = HLA_Genotype_Group), 
                width = 0.05, linewidth = 0.05)+
  scale_colour_manual(values = palette3)+
  ylab("VH:CrD")+ mytheme+
  theme(legend.position="top",
        legend.title=element_blank(),
        axis.title.x =element_blank())
A


#---------------------------------------------------------------------------panel B-------------------------------------------------------------------------

i="VH:CrD"

#Computation Two-way ANCOVA
res.aov <- metadata_wide %>% anova_test(challenge ~ baseline + Treatment+HLA_Genotype_Group)
get_anova_table(res.aov)
get_test_label(res.aov,  detailed = TRUE)

#Post-hoc test
#Simple main effect analyses for treatment
# Effect of treatment at each level of exercise
metadata_wide %>%
  group_by(HLA_Genotype_Group) %>%
  anova_test(challenge ~ baseline + Treatment)

# Pairwise comparisons
pwc1 <- metadata_wide %>% 
  group_by(HLA_Genotype_Group) %>%
  emmeans_test(  #Performs pairwise comparisons between groups using the estimated marginal means
    challenge ~ Treatment, covariate = baseline,
    p.adjust.method = "bonferroni"
  )
pwc1 %>% filter(p.adj < 0.05)


#Report
# Line plot
lp <- ggscatter(#ggline(
  get_emmeans(pwc1), x = "HLA_Genotype_Group", y = "emmean", 
  color = "Treatment", palette = c("#6500ff","#ff6500")) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, color = Treatment), 
    width = 0.05
  )+
  xlab("HLA-DQ genotype group") +
  ylab(paste(i, "at PGC \n estimated marginal means", sep=" "))+
  theme(legend.position="top",
        legend.title=element_blank())+ mytheme

lp

# Comparisons between treatment group at each exercise level
pwc1 <- pwc1 %>% rstatix::add_xy_position(x = "HLA_Genotype_Group", fun="mean_sd")
#pwc1$y.position <- 2*seq(from=1, to=1.3, length.out=nrow(pwc1))
pwc1.filtered <- pwc1 %>% filter(p.adj < 0.05)

library(numform)
pwc1.filtered %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->pwc1.filtered


B <- lp + 
  geom_signif(data=pwc1.filtered,
              aes(xmin=xmin, xmax=xmax, annotations=p_new, y_position=y.position),
              textsize = 6/.pt, 
              manual=TRUE, 
              parse=T, size=0.3, tip_length = 0) 
B <- B +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

#---------------------------------------------------------------------------panel C-------------------------------------------------------------------------
metadata_wide1 <- metadata_wide[metadata_wide$Treatment == "drug",]

res.aov2 <- metadata_wide1 %>% anova_test(challenge ~ baseline+HLA_Genotype_Group)
get_anova_table(res.aov2)
get_test_label(res.aov2,  detailed = TRUE)



gg <- ggplot(metadata_wide1, aes(x = baseline,  y = challenge,color = HLA_Genotype_Group)) +
  geom_point(aes(color=HLA_Genotype_Group), size=0.5) +
  geom_smooth(method = "lm",  mapping = aes(y = predict(lm(challenge ~ baseline + HLA_Genotype_Group,data = metadata_wide1), metadata_wide1)), linewidth=0.5) +
  scale_color_manual(values = palette3) +
  ylab("VH:CrD in PGCd group") +
  xlab("VH:CrD in GFDd group") +  mytheme +
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
   annotate(geom = "text", label = pwc2$p_new, size = 6/.pt, manual=TRUE, parse=T, y = 1:3, x = round(max(pwc2$conf.high)*4,2)) + # p-value. The y coordinates are set by eye
  ylab("Contrast") +
  coord_cartesian(xlim = c(round(min(pwc2$conf.low)*1.05,2),round(max(pwc2$conf.high)*5,2))) +
  scale_x_continuous(position="top") +
  theme_pubr() + mytheme + theme(axis.title.y = element_blank())

gg_effect

Cp <- plot_grid(gg_effect,
                gg ,
                nrow=2,
                align = "v",
                axis = "rl",
                rel_heights = c(0.4, 1)
)
Cp

#---------------------------------------------------------------------------- E panel (ind genes)--------------------------------------------------------------------------------

IDs <-c('ENSG00000084674',
        'ENSG00000110244',
        'ENSG00000118137',
        'ENSG00000169903',
        'ENSG00000106541',
        'ENSG00000148773',
        'ENSG00000179583',
        'ENSG00000132170',
        'ENSG00000007171',
        'ENSG00000271503',
        "ENSG00000198959",
        "ENSG00000115415",
        "ENSG00000117228",
        "ENSG00000117724") 


#prepare data
normalized_counts <- counts(CDraw, normalized=TRUE)
anno <- as.data.frame(CDres[1])

melt<-reshape2::melt(normalized_counts[rownames(normalized_counts) %in% IDs,])
melt$variable <- rownames(melt)
melt$log10_Expr <- log10(melt$value+1)
#melt$Var1 <- factor(melt$Var1,levels=IDs)
melt <-as.data.frame(merge(melt, metadata[,c("ID","label3","HLA_Genotype_Group")],by.x="Var2", by.y="ID"))
melt <- merge(melt, anno[,c(1,10)],by.x="Var1", by.y="PGCp.VS.GFDp.ID")
melt <- melt[melt$HLA_Genotype_Group != "Not identified" & melt$label3 != "GFD",]
melt <- melt[!(melt$PGCp.VS.GFDp.external_gene_name %in% c("APOA4","CCL5")),] # %in% c("NOS2","APOA4","CCL5")),] c("NOS2","APOA4","CCL5","APOA1","TM4SF4","AGR2","CENPF","GBP1","TGM2"))
melt$PGCp.VS.GFDp.external_gene_name <- factor(melt$PGCp.VS.GFDp.external_gene_name, levels =  c("APOB","APOA1","TM4SF4",
                                                                                                "AGR2","MKI67","CENPF",
                                                                                                "STAT1","GBP1","TGM2",
                                                                                                "CIITA","PPARG","NOS2")) #c("APOB","MKI67","STAT1"))

E <- ggplot(melt, aes(x = label3, y = value)) + 
  geom_jitter(aes(color = HLA_Genotype_Group), size = 0.5,position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) + ylab("mRNA Expression, counts") +
  facet_wrap(. ~ PGCp.VS.GFDp.external_gene_name, scales="free", ncol = 3) +
  stat_summary(aes(color = HLA_Genotype_Group), fun.data="mean_sdl", fun.args = list(mult=1), size = 0.4, position = position_dodge(0.8))+
  scale_color_manual(values =  palette3)+ 
  mytheme + theme(axis.title.x = element_blank())
E

#---------------------------------------------------------------------------panel D-------------------------------------------------------------------------


metadata_wide2 <- reshape2::dcast(as.data.frame(metadata[,c("newID", "Timepoint", "epithelial_IFNg_responce", "HLA_Genotype_Group", "Treatment" )]), newID + HLA_Genotype_Group +Treatment ~ Timepoint, value.var = "epithelial_IFNg_responce")

metadata_wide2$baseline <- as.numeric(metadata_wide2$baseline)
metadata_wide2$challenge <- as.numeric(metadata_wide2$challenge)
metadata_wide2 <- metadata_wide2[metadata_wide2$HLA_Genotype_Group != "Not identified", ]
metadata_wide2 <- na.omit(metadata_wide2)



#Two-way ANCOVA
res.aov <- metadata_wide2 %>% anova_test(challenge ~ baseline + Treatment + HLA_Genotype_Group)
get_anova_table(res.aov)


#Post-hoc test

# Pairwise comparisons
pwc1 <- metadata_wide2 %>% 
  group_by(HLA_Genotype_Group) %>%
  emmeans_test(
    challenge ~ Treatment, covariate = baseline,
    p.adjust.method = "bonferroni"
  )
pwc1 %>% filter(p.adj < 0.05)

#Report
# Line plot
lp <- ggscatter(
  get_emmeans(pwc1), x = "HLA_Genotype_Group", y = "emmean", 
  color = "Treatment", palette = c("#6500ff","#ff6500")) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high, color = Treatment), 
    width = 0.1
  )+
  xlab("DQ genotype group") +
  ylab("Epithelial responce to IFNg GSZ \n estimated marginal means")+
  theme(legend.position="top",
        legend.title=element_blank())+ mytheme

#lp

# Comparisons between treatment group at each exercise level
pwc1 <- pwc1 %>% rstatix::add_xy_position(x = "HLA_Genotype_Group", fun="mean_sd")
#pwc1$y.position <- 2*seq(from=1, to=1.3, length.out=nrow(pwc1))
#pwc1.filtered <- pwc1 %>% filter(p.adj < 0.05)

library(numform)
pwc1 %>% 
  mutate(p_new = ifelse(p.adj > 0.01, c(paste("italic('P')~`=", f_num(p.adj,2), "`")), p.adj))%>% 
  mutate(p_new = ifelse(p.adj < 0.01, c(paste("italic('P')~`=", f_num(p.adj,3), "`")), p_new)) %>%
  mutate(p_new = ifelse(p.adj < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))%>% 
  mutate(p_full = ifelse(`p.adj` > 0.01, round(p.adj,2), `p.adj`))%>% 
  mutate(p_full = ifelse(`p.adj` < 0.01, round(p.adj,3), p_full)) %>%
  mutate(p_full = ifelse(`p.adj` < 0.001, formatC(p.adj, format = "e", digits = 2),p_full))->pwc1


D <- lp + 
  geom_signif(data=pwc1,
              aes(xmin=xmin, xmax=xmax, annotations=p_new, y_position=y.position),
              textsize = 6/.pt, 
              manual=TRUE, 
              parse=T, size=0.3, tip_length = 0) #+
D <- D +   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

#----------------------------------------------------------------------------combine all together---------------------------------------------------------------------------------
p0 <- ggarrange(A,B,Cp,  ncol = 1,nrow = 3)
p0
p0_2 <- ggarrange(D,E,  ncol = 1,nrow = 2, heights = c(1,2))
p0_2
p1 <- ggarrange(p0,p0_2,
                ncol = 2, align = "hv", heights = c(1,1))
p1
p <- p1 +                        
  draw_plot_label(label = c("a", "b","c", "d", "e"), size = 11,
                  x = c(0, 0,0,0.5,0.5), 
                  y = c(1,0.66,0.33,1,0.66))
p
graph2ppt(file="./output_data/Figure_5.pptx", width=7.08, height=7.28)



