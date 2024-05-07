library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(ggplotify)
library(dplyr)
library(gridExtra)
library(cowplot)
library(mGSZ)
library(numform)
library(rlist)
library(export)





CD <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
metadata <-colData(CD)

#create theme for all plots
mytheme <- theme_classic(base_family='sans') + theme(axis.text.x = element_text(size=6, color = "black"),
                                                     axis.text.y= element_text(size=6, color = "black"),
                                                     axis.title.y = element_text(size=7, color = "black"),
                                                     axis.title.x = element_text(size=7, color = "black"),
                                                     legend.text = element_text( size = 5, color = "black"),
                                                     legend.key.size = unit(5, 'mm'),
                                                     legend.position = "top",
                                                     legend.title = element_blank(),
                                                     plot.title = element_text(hjust = 0.5, size = 5,face = "bold"),
                                                     axis.line = element_line(colour = "black", linewidth=1/.pt),
                                                     panel.background = element_rect(fill = 'White'),
                                                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     strip.text.x = element_text(size = 5.5, colour = "black"),
                                                     axis.title.x.top = element_text(color='#4daf4a',size=6), 
                                                     axis.text.x.top = element_text(color='#4daf4a',size=5))

palette_bio <- c("#00b347","#006729","#6500ff","#ff6500")
#-------------------------------------------------------------------------------panel A (heatmap)-----------------------------------------------------------------------
results1 <- CDres[["PGCp VS GFDp"]]

terms_to_show2 <- c("regulation of innate immune response", "defense response to other organism", "regulation of inflammatory response",
                    "response to type I interferon", "type II interferon-mediated signaling pa...", "cellular response to organic substance")
goEnrichment1 <- CDres[["PGCp_VS_GFDp_ALL_GO_BP_terms"]]
goEnrichment1 <- goEnrichment1[goEnrichment1$Term %in% terms_to_show2,]

GOplottingobject <- function (terms, category, ID, term,adj_pval,genes, sep,FC, identifier, log2FC) {
  names(FC)[names(FC) == identifier] <- "Name"
  names(FC)[names(FC) == log2FC] <- "logFC"
  names(terms)[names(terms) == category] <- 'category'
  names(terms)[names(terms) == ID] <- 'ID'
  names(terms)[names(terms)== term] <- 'term'
  names(terms)[names(terms) == adj_pval] <- 'adj_pval'
  names(terms)[names(terms) == genes] <- 'genes'
  tgenes <- strsplit(as.vector(terms$genes), sep)
  if (length(tgenes[[1]]) == 1) 
  tgenes <- strsplit(as.vector(terms$genes), ",")
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) FC$logFC[match(x, 
                                                             FC$Name)])
  if (class(logFC) == "factor") {
    logFC <- gsub(",", ".", gsub("\\.", 
                                 "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1
  zsc <- c()
  for (c in 1:length(count)) {
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 
                                                   1, -1))
    zsc <- c(zsc, sum(value)/sqrt(count[c]))
    s <- e + 1
  }
  if (is.null(terms$id)) {
    df <- data.frame(category = rep(as.character(terms$category), 
                                    count), term = rep(as.character(terms$term), count), 
                     ID = rep(as.character(terms$ID), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), 
                     logFC = logFC, adj_pval = rep(terms$adj_pval, count), 
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  else {
    df <- data.frame(category = rep(as.character(terms$category), 
                                    count), ID = rep(as.character(terms$id), count), 
                     term = rep(as.character(terms$term), count), count = rep(count, 
                                                                              count), 
                     ID = rep(as.character(terms$ID), count),genes = as.character(unlist(tgenes)), 
                     logFC = logFC, adj_pval = rep(terms$adj_pval, count), 
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

circa <- GOplottingobject(terms=goEnrichment1, category="Category", ID="ID", term="Term",adj_pval="adj_pval",genes="genes", sep=", ",
                          FC=results1, identifier="ID", log2FC="log2FoldChange")
circa <- circa[order(circa$ID),]
circa <- circa[circa$genes %in% results1[results1$Gene.type == "DE",]$ID,]
circa <- circa[, c("genes","ID","term")]
notannot <- data.frame(genes = results1[!(results1$ID %in% circa$genes) & results1$Gene.type == "DE",]$ID, 
                       ID = rep("OBP", nrow(results1[!(results1$ID %in% circa$genes) & results1$Gene.type == "DE",])),
                       term = rep("other biological processes", nrow(results1[!(results1$ID %in% circa$genes) & results1$Gene.type == "DE",])))
circa <- rbind(circa, notannot)
length(unique(circa$genes)) #95
circa <- merge(circa, results1[,c("ID","external_gene_name")], by.x = "genes", by.y = "ID")
IDs <- c("GO:0034340",
         "GO:0045088", 
         "GO:0050727", 
         "GO:0060333",
         "GO:0071310",
         "GO:0098542",
         "OBP" )


circa <- circa[order(match(circa$ID,IDs)),]
#circa <- circa[circa$genes %in% results1[results1$Gene.type == "DE",]$external_gene_name,]

#DESeq2 size factor normalization
ddsobj <- estimateSizeFactors(CD)
normalized_counts <- counts(ddsobj, normalized=TRUE) 
counts1<- normalized_counts[ circa$genes,]


#define the order of rows if needed
counts1<-counts1[match(rownames(counts1),circa$genes),]

#scale the data
m.data_sc<-t(scale(t(counts1)))

#sort columns
metadata <- colData(CD)
metadata <- metadata[order(metadata$label3, metadata$VHCrD),]
m.data_sc <- m.data_sc[,metadata$ID]

#create annotation 1(Sample type)
type1<- metadata$label2
type3 <- metadata$label
type2 <- circa$ID


top = HeatmapAnnotation(`VH:CrD` = anno_points(metadata$VHCrD, size = unit(0.5, "mm"), axis_param = list(side = "right", gp = gpar(fontsize = 3))),
                        annotation_name_gp = gpar(fontsize = 5),
                        
                        Samples = anno_block(gp = gpar(fill = palette_bio, lty="blank"),
                                             labels = c("GFDd","GFDp", "PGCd","PGCp"), 
                                             labels_gp = gpar(col = "white", fontsize = 5)),
                        
                        gap = unit(1, "mm"),
                        height = unit(0.6, "cm"))

row2 = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#3b565a", "#729ca3","#36b8ea", "#52a68a","#cd8162","#ff6666","#67001F"), lty="blank")), width = unit(0.3, "cm"))


termdf<-unique(circa[,2:3])
termdf$col <- c("#3b565a", "#729ca3","#36b8ea", "#52a68a","#cd8162","#ff6666","#67001F")
termdf$fontsize <- c(rep(4, 7))
termdf$fontfamily <- c(rep("sans", 7))
rownames(termdf) <- termdf$ID
termdf <- termdf[,-1]
colnames(termdf)[1] <- "text"

termlist <- list()
for(i in 1:nrow(termdf)){
  termlist[[i]] <- termdf[i,]
}

names(termlist) <- rownames(termdf)
GOnames = rowAnnotation(textbox = anno_textbox(align_to = type2, termlist, by = "anno_block", background_gp = gpar(fill = "#fcfcfc", col = "white")))


#create color scheme for heatmap
col_fun = colorRamp2(c(min(m.data_sc), quantile(m.data_sc, 0.25), 0, quantile(m.data_sc, 0.75),  max(m.data_sc)), 
                     c("#B2182B","#F4A582" ,"#F7F7F7","#92C5DE", "#2166AC"))

ht_list = Heatmap(m.data_sc, name = "Row Z-Score of Expression", top_annotation = top, #right_annotation = row,
                  left_annotation = row2 ,
                  right_annotation = GOnames,
                  show_column_names = F, col=col_fun,
                  column_order= metadata$ID,
                  column_title_gp = gpar(fontsize = 1, col="White"),
                  #row_names_side = "left",
                  cluster_rows = T,# row_dend_reorder = F, 
                  show_row_dend = FALSE,
                  column_split = factor(type1, levels = c("GFDd","GFDp", "PGCd","PGCp")),
                  column_gap = unit(0.5,"mm"),
                  row_split = factor(type2, levels=unique(circa$ID)),
                  row_gap = unit(0.5,"mm"),
                  cluster_row_slices = FALSE, 
                  show_row_names = F,
                  row_title_gp = gpar(fontsize = 5, direction = "horizontal"),
                  row_title_rot = 0,
                  #col_names_gp = gpar(fontsize = 30),
                  heatmap_legend_param = list(
                    at = c(ceiling(min(m.data_sc)):floor(max(m.data_sc))),
                    labels = c(ceiling(min(m.data_sc)):floor(max(m.data_sc))),
                    title = "Row Z-Score",
                    labels_gp = gpar(fontsize = 5),
                    title_gp = gpar(fontsize = 5), direction = "vertical",
                    legend_width=unit(2,"cm"),grid_height = unit(0.2, "cm"),
                    title_position = "topcenter"))#+row #+row2

#hm1 <- grid.grabExpr(draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(1, "mm"),auto_adjust = FALSE))
#draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", row_gap = unit(0.5, "mm"),auto_adjust = FALSE)

A <- grid.grabExpr(draw(ht_list, merge_legend = F, heatmap_legend_side = "right", 
                        annotation_legend_side = "right", auto_adjust = FALSE))

as.ggplot(A)
#terms table
#termdf<-unique(circa[,2:3])
#colnames(termdf) <- c("ID", "term")
#stable.p <- ggtexttable(termdf, rows = NULL,theme = ttheme(
#  colnames.style = colnames_style(color = "black", size=5),
#  tbody.style = tbody_style(color = "black", size=5)
#))
#stable.p

#-----------------------------------------------------------------------------panel B (reactome barplot)---------------------------------------------------------------------------------
Reactome1 <- CDres[["PGCp_VS_GFDp_ALL_Reactome_terms"]]
Reactome2 <- CDres[["PGCp_VS_PGCd_ALL_Reactome_terms"]]
common <- c("Interferon Signaling","Interferon alpha/beta signaling","Interferon gamma signaling","Growth hormone receptor signaling", "PPARA activates gene expression",                                        
            "Metabolism of vitamins and cofactors")
Reactome2 <- Reactome2[Reactome2$Description %in% common,]
Reactome1 <- Reactome1[Reactome1$Description %in% common,]


#plot #1
Reactome1$Category <- "R"
circa <- GOplottingobject(terms=Reactome1, category="Category", ID="ID", term="Description",adj_pval="p.adjust",genes="geneID", sep="/",
                          FC=results1, identifier="external_gene_name", log2FC="log2FoldChange")
circa <- circa[order(circa$ID),]
circa$direction = NA
circa$direction[circa$logFC>0] = 'up'
circa$direction[circa$logFC<0] = 'down'
circa$`-log10(adj.p-value)` <- (-log10(as.numeric(as.character(circa$adj_pval))))
circa %>% arrange(desc(`-log10(adj.p-value)`)) -> circa
circa$group <- rep(1, dim(circa)[1])


# prepare data for plotting
plotting_df <-   circa %>% group_by(term, direction, `-log10(adj.p-value)`) %>% 
  dplyr::summarise(Freq = n()) %>% mutate(Freq = if_else(direction == "down", -Freq, Freq))%>%
  arrange(desc(`-log10(adj.p-value)`))
plotting_df$term1 <- c( "- Interferon Signaling",
                        "- Interferon alpha/beta \nsignaling",
                        "- Interferon gamma \nsignaling",
                        "- Growth hormone receptor\n signaling",
                       "- Growth hormone receptor\n signaling",
                       "- PPARA activates \ngene expression",
                       "- PPARA activates \ngene expression",
                        "- Metabolism of vitamins \nand cofactors",
                       "- Metabolism of vitamins \nand cofactors")

circa <- merge(circa, plotting_df[,c(1,5)], by.x = "term", by.y="term")
circa %>% group_by(term, `-log10(adj.p-value)`, direction) %>% filter(direction == 'up') %>% dplyr::summarize(counts=mean(count),p=n()) -> circ2.counts


# side-by-side plot
p <- plotting_df %>% 
  ggplot(aes(x = reorder(term1,`-log10(adj.p-value)`), y = Freq, group = direction, fill = direction)) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts")+
  scale_fill_manual(values= c("#b21878","#18b252"),
                    name="log2FC",
                    breaks=c("down", "up"),
                    labels=c("down", "up"))  #+
  #geom_text(aes(label=abs(Freq),y=Freq/1.6, x=term1), vjust = 0.5, size=5/.pt)
p
f <- sum(abs(range(na.omit(plotting_df$Freq))))/(diff(range(plotting_df$`-log10(adj.p-value)`))+1)*2.2 #factor for top y-axis (coordinates are flipped) scaling
circa$y <- f*(circa$`-log10(adj.p-value)`-1.2*median(circa$`-log10(adj.p-value)`))
#highlifgt non-significant points
highlight_df <- circa %>% filter(adj_pval > 0.05)

#add scatter plot on top of barplot


bp <-  p + 
  geom_line(data=circa, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group), color="#09e8b5",linewidth = 0.5) +
  geom_point(data=circa, aes(x=reorder(term1,`-log10(adj.p-value)`), y=y, group=group), shape=21, fill="#09e8b5",color="#09e8b5", size=1) + #color="black"
  geom_point(data=highlight_df, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group), color='#8c8f8e', size=1)+
  # scale_y_continuous(sec.axis = sec_axis(~.*0.02+median(circ$`-log10(adj.p-value)`), name = "log10")) +
  theme_bw()+ 
  #geom_text(data=circ2.counts, aes(x = reorder(term1,`-log10(adj.p-value)`), y=20,label = p),size=6/.pt) +
  scale_y_continuous(sec.axis = sec_axis(~.*1/f+1.2*median(circa$`-log10(adj.p-value)`), name = "-log10(FDR)"))+
  ggtitle("PGCp VS GFDp")

bp <- bp + mytheme + theme(axis.title.y = element_blank(), legend.position = c(0.8, 0.3), legend.key.size = unit(0.2, 'cm'),
                           legend.key = element_rect(colour = "transparent", fill = "white"),
                           legend.background  = element_rect(colour = "transparent", fill = "white"))
bp


#plot #2
Reactome2$Category <- "R"
circa2 <- GOplottingobject(terms=Reactome2, category="Category", ID="ID", term="Description",adj_pval="p.adjust",genes="geneID", sep="/",
                          FC=results1, identifier="external_gene_name", log2FC="log2FoldChange")
circa2 <- circa2[order(circa2$ID),]
circa2$direction = NA
circa2$direction[circa2$logFC>0] = 'up'
circa2$direction[circa2$logFC<0] = 'down'
circa2$`-log10(adj.p-value)` <- (-log10(as.numeric(as.character(circa2$adj_pval))))
circa2 %>% arrange(desc(`-log10(adj.p-value)`)) -> circa2
circa2$group <- rep(1, dim(circa2)[1])


# prepare data for plotting
plotting_df1 <-   circa2 %>% group_by(term, direction, `-log10(adj.p-value)`) %>% 
  dplyr::summarise(Freq = n()) %>% mutate(Freq = if_else(direction == "down", -Freq, Freq))%>%
  arrange(match(term, plotting_df$term))
plotting_df1$term1 <- c( "- Interferon Signaling",
                        "- Interferon alpha/beta \nsignaling",
                        "- Interferon gamma \nsignaling",
                        "- Growth hormone receptor\n signaling",
                        "- Growth hormone receptor\n signaling",
                        #"- PPARA activates \ngene expression",
                        "- PPARA activates \ngene expression",
                        "- Metabolism of vitamins \nand cofactors",
                        "- Metabolism of vitamins \nand cofactors")

#plotting_df <- plotting_df[order(plotting_df$`-log10(adj.p-value)`, decreasing = F),]
plotting_df1$term1 <- factor(plotting_df1$term1, levels=unique(plotting_df[order(plotting_df$`-log10(adj.p-value)`, decreasing = F),]$term1))

# side-by-side plot
p <- plotting_df1 %>% 
  ggplot(aes(x = term1, y = Freq, group = direction, fill = direction)) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts")+
  scale_fill_manual(values= c("#b21878","#18b252"),
                    name="log2FC",
                    breaks=c("down", "up"),
                    labels=c("down", "up"))+
  #geom_text(aes(label=abs(Freq),y=Freq/1.6, x=term1), vjust = 0.5, size=5/.pt)+
  scale_x_discrete(position = "top")
f <- sum(abs(range(na.omit(plotting_df$Freq))))/(diff(range(plotting_df$`-log10(adj.p-value)`))+1)*2.2 #factor for top y-axis (coordinates are flipped) scaling

circa2 <- merge(circa2, plotting_df[,c(1,5)], by.x = "term", by.y="term")
circa2 %>% group_by(term1, `-log10(adj.p-value)`, direction) %>% filter(direction == 'up') %>% dplyr::summarize(counts=mean(count),p=n()) -> circ2.counts
circa2$y <- f*(circa2$`-log10(adj.p-value)`-1.2*median(circa2$`-log10(adj.p-value)`))
circa2$adj_pval <- round(circa2$adj_pval,2)
#highlifgt non-significant points
highlight_df <- circa2 %>% filter(adj_pval > 0.05)

#add scatter plot on top of barplot


bp2 <-  p + 
  geom_line(data=circa2, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group), color="#09e8b5",linewidth = 0.5) +
  geom_point(data=circa2, aes(x=reorder(term1,`-log10(adj.p-value)`), y=y, group=group), shape=21, fill="#09e8b5",color="#09e8b5", size=1) +
  geom_point(data=highlight_df, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group), color='#8c8f8e', size=1)+
  # scale_y_continuous(sec.axis = sec_axis(~.*0.02+median(circ$`-log10(adj.p-value)`), name = "log10")) +
  theme_bw()+ 
  #geom_text(data=circ2.counts, aes(x = reorder(term1,`-log10(adj.p-value)`), y=20,label = p),size=6/.pt) +
  scale_y_continuous(sec.axis = sec_axis(~.*1/f+1.2*median(circa2$`-log10(adj.p-value)`), name = "-log10(FDR)"))+
  ggtitle("PGCp VS PGCd")

bp2 <- bp2 + mytheme + theme( axis.title.y = element_blank(),#text(size=16),
                             #axis.text.y= element_blank(),
                             #strip.text.x = element_text(size = 5),
                             axis.text.y = element_blank(), legend.position = 'none',
                             #axis.title.x = element_blank(),
                             #axis.title.x.top = element_blank()
                             )

bp2


B <- ggarrange(bp + rremove("xlab"), bp2 + rremove("xlab"), ncol = 2,common.legend = F, widths = c(1,0.5), align = "h")

B <- annotate_figure(B,
                bottom = text_grob("Gene counts", color = "black",hjust = 0, x = 0.5, size = 7))
              
B

#----------------------------------------------------------------------------panel C GO-terms ---------------------------------------------------------------------------------------------



#plot#1
circaGO1 <- unique(GOplottingobject(terms=goEnrichment1[goEnrichment1$Term %in% terms_to_show2,], category="Category", ID="ID", term="Term",adj_pval="adj_pval",genes="genes", sep=", ",
                          FC=results1, identifier="ID", log2FC="log2FoldChange"))
circaGO1 <- circaGO1[order(circaGO1$ID),]
circaGO1$direction = NA
circaGO1$direction[circaGO1$logFC>0] = 'up'
circaGO1$direction[circaGO1$logFC<0] = 'down'
circaGO1$`-log10(adj.p-value)` <- (-log10(as.numeric(as.character(circaGO1$adj_pval))))
circaGO1 %>% arrange(desc(`-log10(adj.p-value)`)) -> circaGO1
circaGO1$group <- rep(1, dim(circaGO1)[1])

# prepare data for plotting
plotting_df <-   circaGO1 %>% group_by(term, direction, `-log10(adj.p-value)`) %>% 
  dplyr::summarise(Freq = n()) %>% mutate(Freq = if_else(direction == "down", -Freq, Freq)) %>% arrange(desc(`-log10(adj.p-value)`))
plotting_df$term
plotting_df$term1 <- c("- defense response to \n other organism",
                       "- defense response to \n other organism",
                       "- interferon-gamma-mediated \n signaling pathway",
                       "- interferon-gamma-mediated \n signaling pathway",
                       "- regulation of \n innate immune response",
                       "- regulation of \n innate immune response",
                       "- regulation of \n inflammatory response",
                       "- regulation of \n inflammatory response",
                       "- response to type I interferon",
                       "- cellular response to \n organic substance",
                       "- cellular response to \n organic substance"
                        )

circaGO1 <- merge(circaGO1, plotting_df[,c(1,5)], by.x = "term", by.y="term")
circaGO1 %>% group_by(term1, `-log10(adj.p-value)`, direction) %>% filter(direction == 'up') %>% dplyr::summarize(counts=mean(count),p=n()) -> circ2.counts

plotting_df$term1 <- factor(plotting_df$term1, levels=c("- cellular response to \n organic substance", "- regulation of \n inflammatory response", "- response to type I interferon",
                                                        "- regulation of \n innate immune response", "- interferon-gamma-mediated \n signaling pathway", "- defense response to \n other organism"))

# side-by-side plot
p <- plotting_df %>% 
  ggplot(aes(x = term1, y = Freq, group = direction, fill = direction)) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts")+
  scale_fill_manual(values= c("#B2182B","#2166AC"),
                    name="log2FC",
                    breaks=c("down", "up"),
                    labels=c("down", "up")) #+
  #geom_text(aes(label=abs(Freq),y=Freq/1.6, x=term1), vjust = 0.5, size=5/.pt)
p
f <- sum(abs(range(plotting_df$Freq)))/(diff(range(plotting_df$`-log10(adj.p-value)`))+1)*0.5 #factor for top y-axis (coordinates are flipped) scaling
circaGO1$y <- f*(circaGO1$`-log10(adj.p-value)`-1.2*median(circaGO1$`-log10(adj.p-value)`))
circaGO1$adj_pval <- round(circaGO1$adj_pval,3)
highlight_df <- circaGO1 %>% filter(adj_pval >= 0.05)

#add scatter plot on top of barplot
bp <-  p + 
  geom_line(data=circaGO1, aes(x=reorder(term1,`-log10(adj.p-value)`), y=y, group=group),color="#09e8b5",size = 0.5) +
  geom_point(data=circaGO1, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group),shape=21, color="#09e8b5", fill="#09e8b5", size=1) +
  geom_point(data=highlight_df, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group), color='#8c8f8e', size=1)+
  theme_bw()+ 
  scale_y_continuous(sec.axis = sec_axis(~.*1/f+1.2*median(circaGO1$`-log10(adj.p-value)`), name = "-log10(FDR)"))+
  ggtitle("PGCp VS GFDp")

bp <- bp + mytheme + theme(#axis.title.x.top = element_text(color='#4daf4a',size=10), 
                            # axis.text.x.top = element_text(color='#4daf4a',size=8), 
                             axis.title.y = element_blank(), legend.position = c(0.9, 0.3),legend.key.size = unit(0.2, 'cm'),
                             legend.key = element_rect(colour = "transparent", fill = "white"),
                             legend.background  = element_rect(colour = "transparent", fill = "white")#,text(size=16),
                             #strip.text.x = element_text(size = 6)
                             )
bp

#plot#2
results3 <- CDres[["PGCp VS PGCd"]]
goEnrichment3 <- CDres[["PGCp_VS_PGCd_ALL_GO_BP_terms"]]
circaGO3 <- unique(GOplottingobject(terms=goEnrichment3[goEnrichment3$Term %in% terms_to_show2,], category="Category", ID="GO.ID", term="Term",adj_pval="adj_pval",genes="genes", sep=", ",
                          FC=results3, identifier="ID", log2FC="log2FoldChange"))
circaGO3 <- circaGO3[order(circaGO3$ID),]
circaGO3$direction = NA
circaGO3$direction[circaGO3$logFC>0] = 'up'
circaGO3$direction[circaGO3$logFC<0] = 'down'
circaGO3$`-log10(adj.p-value)` <- (-log10(as.numeric(as.character(circaGO3$adj_pval))))
circaGO3 %>% arrange(desc(`-log10(adj.p-value)`)) -> circaGO3
circaGO3$group <- rep(1, dim(circaGO3)[1])

# prepare data for plotting
plotting_df1 <-   circaGO3 %>% group_by(term, direction, `-log10(adj.p-value)`) %>% 
  dplyr::summarise(Freq = n()) %>% mutate(Freq = if_else(direction == "down", -Freq, Freq)) %>% arrange(match(term, plotting_df$term))
plotting_df1$term1 <- c("- defense response to \n other organism",
                        "- defense response to \n other organism",
                        "- interferon-gamma-mediated \n signaling pathway",
                        "- interferon-gamma-mediated \n signaling pathway",
                        "- regulation of \n innate immune response",
                        "- regulation of \n innate immune response",
                        "- regulation of \n inflammatory response",
                        "- regulation of \n inflammatory response",
                        "- response to type I interferon",
                        "- response to type I interferon",
                        "- cellular response to \n organic substance",
                        "- cellular response to \n organic substance")


#plotting_df <- plotting_df[order(plotting_df$`-log10(adj.p-value)`, decreasing = F),]
plotting_df1$term1 <- factor(plotting_df1$term1, levels=levels(plotting_df$term1))

circaGO3 <- merge(circaGO3, plotting_df[,c(1,5)], by.x = "term", by.y="term")
circaGO3 %>% group_by(term1, `-log10(adj.p-value)`, direction) %>% filter(direction == 'up') %>% dplyr::summarize(counts=mean(count),p=n()) -> circ2.counts

# side-by-side plot
p <- plotting_df1 %>% 
  ggplot(aes(x = term1, y = Freq, group = direction, fill = direction)) +
  geom_bar(stat = "identity", width = 0.75) +  coord_flip()+
  ylab("Gene counts")+
  scale_fill_manual(values= c("#B2182B","#2166AC"),
                    name="log2FC",
                    breaks=c("down", "up"),
                    labels=c("down", "up")) +
  #geom_text(aes(label=abs(Freq),y=Freq/1.6, x=term1), vjust = 0.5, size=5/.pt)+
  scale_x_discrete(position = "top")

f <- sum(abs(range(plotting_df$Freq)))/(diff(range(plotting_df$`-log10(adj.p-value)`))+1)*0.11 #factor for top y-axis (coordinates are flipped) scaling
circaGO3$y <- f*(circaGO3$`-log10(adj.p-value)`-1.2*median(circaGO3$`-log10(adj.p-value)`))
highlight_df <- circaGO3 %>% filter(adj_pval > 0.049)

#add scatter plot on top of barplot
bp2 <-  p + 
  geom_line(data=circaGO3, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group),color="#09e8b5",size = 0.5) +
  geom_point(data=circaGO3, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group),shape=21, color="#09e8b5", fill="#09e8b5", size=1) +
  geom_point(data=highlight_df, aes(x=reorder(term1,`-log10(adj.p-value)`),y=y, group=group), color='#8c8f8e', size=1)+
  theme_bw()+ 
  scale_y_continuous(sec.axis = sec_axis(~.*1/f+1.2*median(circaGO3$`-log10(adj.p-value)`), name = "-log10(FDR)"))+
  ggtitle("PGCp VS PGCd")

bp2 <- bp2 + mytheme + theme(#axis.title.x.top = element_text(color='#4daf4a',size=10), 
                           #axis.text.x.top = element_text(color='#4daf4a',size=8), 
                           axis.title.y = element_blank(),#text(size=16),
                           #strip.text.x = element_text(size = 6),
                           axis.text.y = element_blank(), legend.position = "none")
bp2

Cpl <- ggarrange(bp + rremove("xlab"), bp2 + rremove("xlab"), ncol = 2,common.legend = F, widths = c(1,0.5), align = "h")

Cpl <- annotate_figure(Cpl,
                     bottom = text_grob("Gene counts", color = "black",hjust = 0, x = 0.5, size = 7))

Cpl



#---------------------------------------------------------------------------panel E--------------------------------------------------------------------------------------------------
SampleGSZ <- CDres[["SampleGSZ"]]
groupGSZ <- CDres[["groupGSZ"]]

gene.sets <- c("Transit amplifying cells", "Mature enterocytes" , "Immune cells", "Duodenal transporters")
SampleGSZ <- SampleGSZ[,c("Samples", gene.sets)]
groupGSZ <- groupGSZ[groupGSZ$gene.sets %in% gene.sets,]

groupGSZ %>% 
  mutate(p_full = ifelse(`pvalue` > 0.01, round(pvalue,2), `pvalue`))%>% 
  mutate(p_full = ifelse(`pvalue` < 0.01, round(pvalue,3), p_full)) %>%
  mutate(p_full = ifelse(`pvalue` < 0.001, formatC(pvalue, format = "e", digits = 2),p_full))->groupGSZ

GSZ_melt <- reshape2::melt(SampleGSZ)
colnames(GSZ_melt)[2] <- "gene.sets"
GSZ_melt <- GSZ_melt[GSZ_melt$gene.sets %in% gene.sets,]
GSZ_melt <- merge(GSZ_melt,as.data.frame(metadata), by.x="Samples", by.y = "ID")

groupGSZ <- groupGSZ[groupGSZ$gene.sets %in% gene.sets,]
groupGSZ <- na.omit(groupGSZ[groupGSZ$pvalue <= 0.05,])


GSZ_melt$gene.sets <- factor(as.character(GSZ_melt$gene.sets), levels=gene.sets)
groupGSZ$gene.sets <- factor(as.character(groupGSZ$gene.sets), levels=gene.sets)

#boxplot
df <-   as.data.frame(GSZ_melt) %>% group_by(gene.sets, label3) %>% 
  dplyr::summarise(max=max(value)*1.2)
df$y_pos <- c(13.8,14,13.5,
              17,16,16.5,
              16.5, 16.5, 15.8,
              4.9,4.26,4.6)
groupGSZ <- merge(groupGSZ,df, by.x=c("gene.sets", "group2"), by.y=c("gene.sets", "label3"))

bp<-ggboxplot(as.data.frame(GSZ_melt), x = "label3", y = "value", outlier.colour = NA, order = c("GFD","PGCd","PGCp"),
              ylab = "GSZ",palette = alpha(c("#00b347","#6500ff","#ff6500"),1),
              fill="label3",facet.by = 'gene.sets', ncol =  4,scales = "free_y")+
  scale_x_discrete(labels=c("GFDd+p","PGCd","PGCp"))+
  geom_jitter(size=0.5, alpha=0.5)+
  mytheme+
  geom_signif(data=groupGSZ,
              aes(xmin=group1, xmax=group2, annotations=p_new, y_position=y_pos),
              textsize = 6/.pt, 
              manual=TRUE, parse=T, size=0.3)

E <- bp +
  theme( legend.position = "none", axis.title.x = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

E

#----------------------------------------------------------------------------combine all together ver2---------------------------------------------------------------------------------
p1 <- ggarrange(B, Cpl, ncol = 2, widths = c(1,1), align = "h")
p1

p <- arrangeGrob(A, p1, E,
                  ncol = 1, nrow = 3, 
                  layout_matrix = rbind(c(1),c(2),c(3)))

p <- as_ggplot(p) +                        # transform to a ggplot
  draw_plot_label(label = c("a", "b", "c","d"), size = 11,
                  x = c(0, 0,0.5,0), 
                  y = c(1, 0.67, 0.67,0.35))
p

graph2ppt(file="./output_data/Figure_2.pptx", width=7.08, height=7.28)
