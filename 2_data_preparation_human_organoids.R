library(readxl)
library(openxlsx)
library(reshape2)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(rlist)

#----------------------------------------------------------------------------Org data--------------------------------------------------------------------------------
metadata <- as.data.frame(readxl::read_xlsx("./input_data/metadata_hDUO.xlsx"))

counts <- as.data.frame(read_xlsx("./input_data/gene_count_hDUO.xlsx"))
counts %>% mutate_at(c(metadata$`Sample_Name`), as.numeric) -> counts
row.names(counts) <- counts$gene_id
counts2 <- as.matrix(counts[,c(metadata$`Sample_Name`)]) #  file created in chunk1 


seobj <- SummarizedExperiment(assays=counts2, colData=as.data.frame(metadata))
ddsobj <- DESeqDataSet(seobj,design = ~ Treatment+RIN+passage+RNA_isolation_batch+Pair)

#pre-filtering to keep only rows that have at least 575 reads total
keep <- rowSums(counts(ddsobj)) >= ncol(ddsobj)*5 #ncol(ddsobj)*5
ddsobj <- ddsobj[keep,]
ddsobj <- estimateSizeFactors(ddsobj)
#normalized_counts <- counts(ddsobj, normalized=TRUE)
saveRDS(ddsobj, "./input_data/hDuo ZED1277 DESeqDataSet object.rds")

group_comp <- data.frame(#`No`= c(1:length(combn(unique(metadata$Treatment),2)[2,])),
  `Sample Group` = combn(unique(metadata$Treatment),2)[2,],
  `Control Group` = combn(unique(metadata$Treatment),2)[1,],
  `Group Comparisons` = paste(combn(unique(metadata$label),2)[2,], "VS", combn(unique(metadata$label),2)[1,], sep="_"))
group_comp <- group_comp[-c(4,5),]
group_comp$No <- c(1:nrow(group_comp))
group_comp <- merge(group_comp, unique(metadata[,c("Treatment","label")]), by.x="Sample.Group",by.y="Treatment",all.x = TRUE)
group_comp <- merge(group_comp, unique(metadata[,c("Treatment","label")]), by.x="Control.Group",by.y="Treatment",all.x = TRUE)
group_comp <- group_comp[order(group_comp$No),]

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#searchAttributes(mart = human, pattern = "type")
annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
                     filters = "ensembl_gene_id", values = rownames(counts), mart = human)
annot_human <- rbind(annot_human, data.frame(ensembl_gene_id ="ENSG00000112096", external_gene_name = "SOD2", gene_biotype = "protein_coding"))

dds <- DESeq(ddsobj)

results_ORG <- list()
for(i in 1:nrow(group_comp)){
  res <- results(dds, contrast=c("Treatment",group_comp$Sample.Group[i],group_comp$Control.Group[i]))
  resOrdered <- res[order(res$padj),]
  results1 <- as.data.frame(resOrdered)
  results1$ID <- rownames(results1)
  results1$Experiment<-group_comp$Group.Comparisons[i]
  results1$Gene.type = "non-DE"
  results1$Gene.type[results1$padj<=.05 & abs(results1$log2FoldChange)>=0.5] = 'DE'
  results1 <- merge(results1, annot_human, by.x="ID", by.y="ensembl_gene_id")
  results_ORG[[i]] <- results1
}
names(results_ORG) <- group_comp$Group.Comparisons

list.save(results_ORG, './input_data/hDuo ZED1277 DEG results.rds')
write.xlsx(results_ORG[1:4], file = "./output_data/DEGs_in_hDUO_organoids.xlsx")
