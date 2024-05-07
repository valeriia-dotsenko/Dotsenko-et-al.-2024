library(readxl)
library(openxlsx)
library(reshape2)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(rlist)
library(mGSZ)
library(numform)
library(ReactomePA)
library(topGO)
library(org.Hs.eg.db) #annotation package
library(GOxploreR)

#---------------------------------------------------------------------------CD data-----------------------------------------------------------------------------------------------
#prepare data
metadata <- read_xlsx(path = "./input_data/metadata-comb.xlsx")
#----------------------------------------------------------------------------CD patients characteristics--------------------------------------------------------------------------
metadata$VHCrD <- as.numeric(metadata$VHCrD)

PChar <- metadata %>% 
  group_by(Treatment, Timepoint) %>%
  dplyr::summarize(n = n(),mean_Age = mean(Age), sd_Age = sd(Age), sum(Gender== "F"), percent = sum(Gender== "F")/n()*100, sum(Genotype_PCR== "HLA-DQ2"), percent_DQ2 = sum(Genotype_PCR== "HLA-DQ2")/n()*100,
                   sum(Genotype_PCR== "HLA-DQ8"), percent_DQ8 = sum(Genotype_PCR== "HLA-DQ8")/n()*100,
                   sum(Genotype_PCR== "HLA-DQ2 + HLA-DQ8"), percent_DQ2DQ8 = sum(Genotype_PCR== "HLA-DQ2 + HLA-DQ8")/n()*100, 
                   median_TG2IGA = median(`TG2IGA, kU/L`,na.rm=T), `Q1-Q3_TG2IGA` = paste0(quantile(`TG2IGA, kU/L`, na.rm=T)[2], "-",quantile(`TG2IGA, kU/L`, na.rm=T)[4]),
                   mean_VHCrD = mean(VHCrD), sd_VHCrD = sd(VHCrD), 
                   mean_weight = mean(Weight_kg),sd_weight = sd(Weight_kg))%>% 
  mutate_if(is.numeric, round, 2)

metadata_wide <- reshape2::dcast(as.data.frame(metadata[,c("newID", "Timepoint", "VHCrD", "HLA_Genotype_Group", "Treatment" )]), newID + HLA_Genotype_Group +Treatment ~ Timepoint, value.var = "VHCrD")

metadata_wide$baseline <- as.numeric(metadata_wide$baseline)
metadata_wide$challenge <- as.numeric(metadata_wide$challenge)
metadata_wide <- na.omit(metadata_wide)

metadata_wide %>% 
  group_by(Treatment) %>%
  dplyr::summarize(n = n(),mean_baseline = mean(baseline), sd_baseline = round(sd(baseline),2), mean_challenge = mean(challenge), sd_challenge = round(sd(challenge),2),
                   p.val= round(t.test(challenge,baseline, paired=T)$p.value,3),
                   est = round(t.test(challenge,baseline, paired=T)$estimate,2), 
                   est_CI = paste(round(t.test(challenge,baseline, paired=T)$conf.int[1],2), "to", 
                                  round(t.test(challenge,baseline, paired=T)$conf.int[2],2), sep=" ")) %>% 
  mutate_if(is.numeric, round, 2)

writexl::write_xlsx(PChar, path = "./output_data/table_1_patient_characteristics.xlsx")

#---------------------------------------------------------------------------create DESeqDataSet-------------------------------------------------------
counts <- read.table("./input_data/Expression_raw_counts.txt", header=T)
row.names(counts) <- counts$Identifier
metadata <- metadata[metadata$ID!='P1_C88',]
counts <- as.matrix(counts[,metadata$ID]) #  file created in chunk1 

seobj <- SummarizedExperiment(assays=counts, colData=as.data.frame(metadata))
ddsobj <- DESeqDataSet(seobj,design = ~ label2+Pair+RIN)

#pre-filtering to keep only rows that have at least 575 reads total
keep <- rowSums(counts(ddsobj)) >= ncol(ddsobj)*5 #ncol(ddsobj)*5
ddsobj <- ddsobj[keep,]
vst_seobj <- vst(ddsobj)
dds <- DESeq(ddsobj)
#DESeq2 size factor normalization
ddsobj <- estimateSizeFactors(ddsobj)
#normalized_counts <- counts(ddsobj, normalized=TRUE)
saveRDS(ddsobj, "./input_data/Celiac biopsies DESeqDataSet object.rds")

#ExpressionSet object for deconvolution
#change id to names
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human =  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype", "entrezgene_id"),
                     filters = "ensembl_gene_id", values = rownames(counts), mart = human)

tmp <- annot_human[annot_human$ensembl_gene_id %in% rownames(counts), c("ensembl_gene_id","external_gene_name")]
tmp <- tmp[!(tmp$external_gene_name == ''),]
tmp <- tmp[!duplicated(tmp$external_gene_name), ]

counts_name <- counts[rownames(counts) %in% tmp$ensembl_gene_id,]
rownames(counts_name) <- tmp$external_gene_name[match(rownames(counts_name), tmp$ensembl_gene_id)]

rownames(metadata) <- metadata$ID
phenoData <- new("AnnotatedDataFrame", data=metadata)
CDeset <- ExpressionSet(assayData=counts_name, phenoData=phenoData)
saveRDS(CDeset, "./input_data/Celiac biopsies ExpressionSet object.rds")


group_comp <- data.frame(`Sample Group` = c("PGCp","PGCd","PGCp"),
                         `Control Group` = c("GFDp","GFDd","PGCd"),
                         `Group Comparisons` = c("PGCp VS GFDp","PGCd VS GFDd","PGCp VS PGCd"))
results_CD <- list()
for(i in 1:nrow(group_comp)){
  res <- results(dds, contrast=c("label2",group_comp$Sample.Group[i],group_comp$Control.Group[i]))
  resOrdered <- res[order(res$padj),]
  results1 <- as.data.frame(resOrdered)
  results1$ID <- rownames(results1)
  results1$Experiment<-group_comp$Group.Comparisons[i]
  results1$Gene.type = "non-DE"
  results1$Gene.type[results1$padj<=.05 & abs(results1$log2FoldChange)>=0.5] = 'DE'
  results1 <- merge(results1, unique(annot_human[, c("ensembl_gene_id","external_gene_name","gene_biotype")]), by.x="ID", by.y="ensembl_gene_id")
  nonannot <- data.frame(ID=c("ENSG00000237070","ENSG00000237643","ENSG00000260912","ENSG00000261012", "ENSG00000265206"),
                         name = c("LOC105375166","HSALNG0051132","HSALNG0070301","HSALNG0013478","HSALNG0117712"))
  results1[results1$ID %in% nonannot$ID,]$external_gene_name <- nonannot$name
  results_CD[[i]] <- results1
}
names(results_CD) <- group_comp$Group.Comparisons
#list.save(results_CD, './input_data/Celiac biopsies DEG results.rds')
write.xlsx(results_CD[1:3], file = "./output_data/DEGs_in_CeD_biopsies.xlsx")



#-------------------------------------------------------------------------------Reactome enrichment (done once and saved for later)----------------------------
#Reactome enrichment #1
results1 <- results_CD[["PGCp VS GFDp"]]
results1 <- merge(results1, annot_human[, c("ensembl_gene_id","entrezgene_id")], by.x = "ID", by.y = "ensembl_gene_id")
genesLFC <- results1[abs(results1$log2FoldChange)>=0.5 & as.numeric(results1$padj)<0.05 & results1$gene_biotype == "protein_coding" & is.na(results1$padj)==F,]$log2FoldChange
names(genesLFC) <- results1[abs(results1$log2FoldChange)>=0.5 & as.numeric(results1$padj)<0.05 & results1$gene_biotype == "protein_coding"& is.na(results1$padj)==F,]$entrezgene_id

Reactome1 <- as.data.frame(enrichPathway(gene=names(genesLFC), readable=TRUE, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "fdr", qvalueCutoff = 0.5))
#p-values are calculated using a hypergeometric distribution, the test is one-tailed by default. When using the hypergeometric distribution for enrichment analysis, or testing the overrepresentation of items 
#(e.g., genes, GO terms) in a subset, the focus is typically on detecting an excess of items in the subset compared to what would be expected by chance.
#the p-values obtained from the hypergeometric distribution represent the probability of observing the observed enrichment (or more extreme enrichment) by chance alone, 
#specifically in the direction of overrepresentation.

results_CD[["PGCp_VS_GFDp_ALL_Reactome_terms"]] <- Reactome1

#Reactome enrichment #2
results3 <- results_CD[["PGCp VS PGCd"]]
results3 <- merge(results3, annot_human[, c("ensembl_gene_id","entrezgene_id")], by.x = "ID", by.y = "ensembl_gene_id")

results3$log2FoldChange <- as.numeric(as.character(results3$log2FoldChange))
results3$padj <- as.numeric(as.character(results3$padj))
genesLFC2 <- results3[abs(results3$log2FoldChange)>=0.5 & as.numeric(results3$padj)<0.05 & results3$gene_biotype == "protein_coding"& is.na(results3$padj)==F,]$log2FoldChange
names(genesLFC2) <- results3[abs(results3$log2FoldChange)>=0.5 & as.numeric(results3$padj)<0.05 & results3$gene_biotype == "protein_coding"& is.na(results3$padj)==F,]$entrezgene_id

Reactome2 <- as.data.frame(enrichPathway(gene=names(genesLFC2), readable=TRUE, organism = "human", pvalueCutoff = 0.5, pAdjustMethod = "fdr", qvalueCutoff = 0.5))
results_CD[["PGCp_VS_PGCd_ALL_Reactome_terms"]] <- Reactome2
write.xlsx(results_CD[4:5], file = "./output_data/ALL_Reactome_terms_CeD_biopsies.xlsx")

#-------------------------------------------------------------------------------GOenrichment (done once and saved for later)----------------------------
#GO annotation #1

results1 <- results1 %>% mutate_at(c("log2FoldChange","padj"), as.numeric) 

logFC <- 0.27
genes <-results1[abs(results1$log2FoldChange)>=logFC & is.na(results1$padj)==F,]$padj
names(genes) <- results1[abs(results1$log2FoldChange)>=logFC & is.na(results1$padj)==F,]$ID

selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="ensembl") #c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=5)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 
#we perform a classical enrichment analysis by testing the over-representation of GO terms within the
#group of diferentially expressed genes. For the method classic each GO category is tested independently. 
#a one-tailed test is appropriate because the goal is typically to identify GO terms
# that are significantly enriched (overrepresented) in a gene set 
#compared to what would be expected by chance.

goEnrichment1 <- GenTable(GOdata,  classicFisher = resultFisher,topNodes = 1000)
goEnrichment1$adj_pval <- p.adjust(goEnrichment1$classicFisher, "fdr")
goEnrichment1$Category <- rep("BP", dim(goEnrichment1)[1])
#Extract list of genes, assosiated with each term
AnnotatedGenes = sapply(sapply(goEnrichment1$GO.ID, function(x) as.character(unlist(genesInTerm(object = GOdata, whichGO = x)))), paste, collapse=', ') 
goEnrichment1$genes <- AnnotatedGenes
colnames(goEnrichment1)[1] <- c("ID")
results_CD[["PGCp_VS_GFDp_ALL_GO_BP_terms"]] <- goEnrichment1


#GO annotation #2
results3 <- results3 %>% mutate_at(c("log2FoldChange","padj"), as.numeric) 


genes <-results3[abs(results3$log2FoldChange)>=logFC & is.na(results3$padj)==F,]$padj
names(genes) <- results3[abs(results3$log2FoldChange)>=logFC & is.na(results3$padj)==F,]$ID

selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="ensembl") #c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene")
GOdata <- new("topGOdata",
              ontology="BP",
              allGenes=genes,
              annot=annFUN.GO2genes,
              GO2genes=allGO2genes,
              geneSel=selection,
              nodeSize=5)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

goEnrichment3 <- GenTable(GOdata,  classicFisher = resultFisher,topNodes = 1000)
goEnrichment3$adj_pval <- p.adjust(goEnrichment3$classicFisher, "fdr")
goEnrichment3$Category <- rep("BP", dim(goEnrichment3)[1])

#Extract list of genes, associated with each term
AnnotatedGenes = sapply(sapply(goEnrichment3$GO.ID, 
                               function(x) as.character(unlist(genesInTerm(object = GOdata, whichGO = x)))), 
                        paste, collapse=', ') 
goEnrichment3$genes <- AnnotatedGenes

results_CD[["PGCp_VS_PGCd_ALL_GO_BP_terms"]] <- goEnrichment3

write.xlsx(results_CD[6:7], file = "./output_data/ALL_GO_BP_terms_CeD_biopsies.xlsx")

list.save(results_CD, './input_data/Celiac biopsies DEG results.rds')
