library(readxl)
library(reshape2)
library(biomaRt)
library(zellkonverter)
library(SummarizedExperiment)
library(MuSiC)
library(Biobase)
library(DESeq2)
library(dplyr)
library(purrr)
library(ggpubr)
library(ggbreak)
library(gridExtra)
library(cowplot)
library(colorspace)
library(export)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

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
                                                     strip.text.x = element_text(size = 5.5, colour = "black"))

palette_bio <- c("#00b347","#006729","#6500ff","#ff6500")

#----------------------------------------------------------------------panel A deconvolution
#Gut atlas survey data  https://www.gutcellatlas.org/

#sceHI <-  readH5AD("./Full_obj_log_counts_soupx_v2.h5ad", reader = "R") #Memory-heavy, performed once
#metadata <- colData(sceHI)
#metadataDUO <- metadata[metadata$Region.code == "DUO",]
#unique(metadataDUO$category) #[1] Endothelial  Epithelial   Plasma cells Mesenchymal  Neuronal     Myeloid      T cells   
#unique(metadataDUO$Integrated_05) 
#assays(sceHI)
#DUO.only <- sceHI[, sceHI$Region.code == "DUO"]
#colData(DUO.only)
#saveRDS(DUO.only, file = "./Full_obj_log_counts_soupx_v2_DUOdata.rds")


DUO.sce = readRDS("./input_data/Full_obj_log_counts_soupx_v2_DUOdata.rds")
#DUO.sce
names(assays(DUO.sce)) <- c("counts")
#colData(DUO.sce)
#unique(colData(DUO.sce)$Diagnosis)
#unique(colData(DUO.sce)$sample.name)

metaDuo <- colData(DUO.sce)


cell.types <- as.data.frame(metaDuo) %>% dplyr::select(category, Integrated_05) %>% group_by(Integrated_05) #%>% summarise(n = n())
cell.types <- unique(cell.types)
cell.numbers <- as.data.frame(metaDuo) %>% dplyr::select(category, Integrated_05) %>% group_by(Integrated_05) %>% summarise(n = n())

cells <- merge(cell.types, cell.numbers)
## Bulk RNA-seq data

CDeset <-readRDS("./input_data/Celiac biopsies ExpressionSet object.rds")
CDeset


# music deconvolution by category
bulk.mtx <- exprs(CDeset)

prop_music2=music_prop(bulk.mtx = bulk.mtx, sc.sce = DUO.sce, 
                       clusters = 'category', samples = 'sample.name',
                       select.ct = c("Endothelial",  "Epithelial",  "Plasma cells", "Mesenchymal",     "Myeloid",      "T cells"), verbose = F)$Est.prop.weighted

# plot estimated cell type proportions
prop_all = cbind('proportion'=c(prop_music2), 'sampleID'=rep(rownames(prop_music2),times=ncol(prop_music2)), 'celltype'=rep(colnames(prop_music2), each=nrow(prop_music2)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))
proportions1 <- reshape2::dcast(prop_all, sampleID  ~ celltype, value.var = "proportion")



metadata <- pData(CDeset)
prop_all <- merge(prop_all,metadata, by.x = "sampleID", by.y = "ID")

n = length(unique(prop_all$celltype))
h1 <- hcl.colors(n, palette = "Dynamic")

plotA <- ggplot(prop_all, aes(x=label2, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8, size=0.5)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = mean,
               geom = "crossbar", width = 0.6,linewidth=0.1,color='gray36')+
  facet_grid(.~celltype, scales = "fixed")+
  scale_color_manual(values=h1) + mytheme + theme(legend.position = "none", axis.text.x = element_text(angle = 45,hjust=1)) 

plotA <- plotA + scale_y_sqrt() + scale_y_break(c(0.24, 0.72)) + theme(axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank())
#----------------------------------------------------------------------panel B deconvolution
# music deconvolution by cell type
CellTypes <- c("Mature arterial EC", "Mature venous EC", "LEC5 (CLDN11+)", "LEC1 (ACKR4+)", "LEC6 (ADAMTS4+)", 
               "Paneth","Microfold cell","TA", "Enterocyte"," BEST4+ epithelial",  "Tuft", "Stem cells", " Goblet cell ",  "M/X cells (MLN/GHRL+)",
               "Contractile pericyte (PLN+)", "Transitional Stromal 3 (C3+)", "mLN Stroma (FMO2+)",
               "IgA plasma cell", "IgG plasma cell",
               "Macrophages", "gdT", "Activated CD4 T", "TRGV2 gdT")
prop_music1=music_prop(bulk.mtx = bulk.mtx, sc.sce = DUO.sce, 
                       clusters = 'Integrated_05', samples = 'sample.name',
                       select.ct = CellTypes, verbose = F)$Est.prop.weighted

# plot estimated cell type proportions
prop_all = cbind('proportion'=c(prop_music1), 'sampleID'=rep(rownames(prop_music1),times=ncol(prop_music1)), 'celltype'=rep(colnames(prop_music1), each=nrow(prop_music1)))
prop_all = as.data.frame(prop_all)
prop_all$proportion = as.numeric(as.character(prop_all$proportion))
proportions2 <- reshape2::dcast(prop_all, sampleID  ~ celltype, value.var = "proportion")


metadata <- pData(CDeset)
prop_all <- merge(prop_all,metadata, by.x = "sampleID", by.y = "ID")
unique(prop_all$celltype)
n = length(unique(prop_all$celltype))
h1 <- hcl.colors(palette = "TealGrn", n = n)


CellTypes_to_plot <- c( "Stem cells", "TA", "Enterocyte", "Paneth", "Tuft","Microfold cell")
prop_to_plot <- prop_all[prop_all$celltype %in% CellTypes_to_plot,]
prop_to_plot[prop_to_plot$celltype == "Microfold cell", "celltype"] <- "M-cell"

prop_to_plot$celltype <-  factor(prop_to_plot$celltype, levels=c(CellTypes_to_plot[1:5], "M-cell"))
plotB <- ggplot(prop_to_plot, aes(x=label2, y=proportion, color=celltype)) + xlab('')+
  geom_jitter(width=0.25,alpha=0.8, size=0.5)+ylab('Cell Type Proportions')+theme_bw()+
  stat_summary(fun = median,
               geom = "crossbar", width = 0.6,size=0.1,color='gray36')+
  facet_grid(.~celltype, scales = "fixed")+
    scale_color_manual(values=h1) + mytheme + theme(legend.position = "none", axis.text.x = element_text(angle = 45,hjust=1)) 
plotB <- plotB + scale_y_sqrt() + scale_y_break(c(0.24, 0.48)) + theme(axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank()) 

 


#----------------------------------------------------------------------------panelC ind genes--------------------------------------------------

CD <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")

df_markers <- data.frame(Genes = c("OLFM4",'SLC12A2','CDCA7','CDK8',"NOS2",
                                   'RBP2','ALDOB','ALPI','APOC3', "APOA4",
                                   "KLRD1", "GZMA","MICA","HLA-E","KLRK1",
                                   "SLC17A4", "ABCC2", "SLC46A1", "SLC39A5",'SLC7A1'),
                         Cell_type = c(rep("Stem cells",5), rep("Enterocytes",5), rep("cytotoxic IELs",5), rep("Duodenal transporters",5)))
                         
#searchAttributes(mart = human, pattern = "type")
annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
                     filters = "external_gene_name", values = df_markers$Genes, mart = human)
df_markers <- merge(df_markers, annot_human, by.x = "Genes", by.y = "external_gene_name")
#IDs <- c("ENSG00000213809", "ENSG00000180644", "ENSG00000145649", "ENSG00000100453",  "ENSG00000204592", "ENSG00000134545",  "ENSG00000134539","ENSG00000204520" )

ddsobj <- estimateSizeFactors(CD)
normalized_counts <- counts(ddsobj, normalized=TRUE) 
counts <- normalized_counts[rownames(normalized_counts) %in% annot_human$ensembl_gene_id,]
metadata <- colData(CD)


melt<-reshape2::melt(counts)
melt$log10_Expr <- log10(melt$value+1)
melt <-as.data.frame(merge(melt, metadata[,c("ID","label2")],by.x="Var2", by.y="ID"))
melt <-as.data.frame(merge(df_markers, melt, by.x="ensembl_gene_id", by.y="Var1"))

CT <- unique(df_markers$Cell_type)  
  CT %>%
    map(function(x){
      melt1 <- melt[melt$Cell_type == x,] #data.frame(Month = ts_data_1[,1], y = ts_data_1[[x]])
      pl <- ggviolin(melt1, x = "label2", y = "log10_Expr", outlier.colour = NA, order = c("GFDd","GFDp","PGCd","PGCp"), 
                     ylab = "log10(Expression, counts)",
                     palette = palette_bio, 
                     fill="label2", alpha = 0.8,
                     draw_quantiles = 0.5,
                     title = x)+
        facet_wrap("Genes", scales="free_y", ncol = 1)+
        geom_jitter(size=0.5, alpha=0.5)+ mytheme +
        theme(axis.title.x = element_blank(), legend.position = "none" #,axis.text.y=element_blank()
        ) + labs(fill = element_blank()) 
      return(pl)
      #legend <- g_legend(pl + theme(legend.position='bottom'))
    }) %>%
    arrangeGrob(grobs = ., ncol=4) %>% #add ncol/nrow argument here
    grid.arrange()  -> plotC


#----------------------------------------------------------------------------combine all together ver1---------------------------------------------------------------------------------

  
  A <- arrangeGrob(print(plotA),print(plotB),plotC,
                   ncol = 2, nrow = 2, 
                   layout_matrix = rbind(c(1,2),c(3,3)),
                   padding=unit(0.5, 'line'))    
p <-   as_ggplot(A) +                        # transform to a ggplot
    draw_plot_label(label = c("a", "b", "c"), size = 11,
                    x = c(0, 0.5,0), 
                    y = c(1, 1, 0.5))

p
graph2ppt(file="./output_data/Figure_S3.pptx", width=7.08, height=7.28)


MUSiC_res <- list(proportions2,proportions1)
names(MUSiC_res) <- c("Deconvolution_cell_type_results", "Deconvolution_cell_category_res")
openxlsx::write.xlsx(MUSiC_res, file = "./output_data/Deconvolution_results_CeD_biopsies.xlsx")

  