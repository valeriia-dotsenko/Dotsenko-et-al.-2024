library(readxl)
library(rlist)
library(biomaRt)


#------------------------------------------------------------------------------1. HUMAN DUODENAL MARKERs-------------------------------------------------------------------------------------------


# retrieve human duodenal marker genes from https://doi.org/10.1016/j.celrep.2021.108819

#1. open excel file in R
#download the supp file from paper using the link 

URL <- 'https://ars.els-cdn.com/content/image/1-s2.0-S2211124721001339-mmc6.xlsx'
pf <- tempfile()
download.file(URL, pf, mode="wb")
#open file

read_excel_allsheets <- function(filename, tibble = FALSE, skip) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X, skip = skip))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
} #function author Jeromy Anglim https://stackoverflow.com/questions/12945687/read-all-worksheets-in-an-excel-workbook-into-an-r-list-with-data-frames 
mysheets <- read_excel_allsheets(pf, skip=6)

mglist <- vector("list")
mgtable <- data.frame()
for(i in names(mysheets)){
  tmp <- mysheets[[i]]
  genes <- t(as.data.frame(strsplit(tmp[tmp$padj<=.05 & tmp$log2FoldChange>=1.0,]$`Gene name`, "__")))
  if(length(genes)==0){
    genes <- data.frame(V1="NA", V2="NA")
  }
  mglist[[i]] <- genes[,1]
  gene.sets <- rep(i,length(genes))
  tmp2 <- cbind(gene.sets,genes[,1])
  mgtable <- rbind(mgtable,tmp2)
}

mgtable <- unique(mgtable)

#list.save(mglist, 'Human duodenum marker genes list.rds')
#write.table(mgtable, file="Human duodenum marker genes list.txt",row.names = F)


#2. annotate with ENSEMBL ID
#human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
human =  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
#searchAttributes(mart = human, pattern = "type")
annot_human <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                     filters = "external_gene_name", values = mgtable$V2, mart = human)
mgtable <- merge(mgtable, annot_human, by.x="V2", by.y="external_gene_name")

mglist2 <- vector("list")
for(i in unique(mgtable$gene.sets)){
  tmp <- mgtable[mgtable$gene.sets==i,]
  genes <- tmp$ensembl_gene_id
  mglist2[[i]] <- genes
  }
#list.save(mglist2, 'Human duodenum marker genes list ENSEMBL.rds')


#------------------------------------------------------------------------------2. HUMAN DUODENAL TRANSPORTERS-------------------------------------------------------------------------------------------
# retrieve from the same paper https://doi.org/10.1016/j.celrep.2021.108819
#1. open excel file in R
#download the supp file from paper using the link 

URL <- 'https://ars.els-cdn.com/content/image/1-s2.0-S2211124721001339-mmc7.xlsx'
pf <- tempfile()
download.file(URL, pf, mode="wb")
#open file
transporters <- read_excel(pf, skip=3)
annot_humanTR <- getBM(attributes = c("ensembl_gene_id","external_gene_name"),
                       filters = "external_gene_name", values = transporters$GeneID, mart = human)
transporters <- merge(transporters, annot_humanTR, by.x="GeneID", by.y="external_gene_name", all.x = T)
nonannotTR <- data.frame(ID=c("ENSG00000119686", "ENSG00000168394", "ENSG00000140519"),
                         name = c("SLC49A2","ABCB2","SLC42A3"))
transporters[transporters$GeneID %in% nonannotTR$name,]$ensembl_gene_id <- nonannotTR$ID
transporters_duo <- transporters[is.na(transporters$Duodenum)==F,]

transporters_list <- list("Duodenal transporters" = transporters[is.na(transporters$Duodenum)==F,]$ensembl_gene_id,
                          "Esophagus transporters" = transporters[is.na(transporters$Esophagus)==F,]$ensembl_gene_id,
                          "Stomach transporters" = transporters[is.na(transporters$Stomach)==F,]$ensembl_gene_id)
#list.save(transporters_list, 'Human duodenum transporters list ENSEMBL.rds')

#------------------------------------------------------------------------------3. REACTOME PATHWAYS-------------------------------------------------------------------------------------------
reactome <- read.table("./input_data/gene_sets/c2.cp.reactome.v7.4.entrez.txt", header=F, sep = '\t') #downloaded from GSEA
reactome_melt <- na.omit(reshape2::melt(reactome))

annot_human <-getBM(attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                    filters = "entrezgene_id", values = unique(reactome_melt$value), mart = human)
reactome_melt <- merge(reactome_melt, annot_human, by.x="value", by.y="entrezgene_id")

#transform list to HUMAN ENSEMBL
reactome_list <- vector("list")
for(i in unique(reactome_melt$V1)){
  genes <- reactome_melt[reactome_melt$V1==i,]$ensembl_gene_id
  reactome_list[[i]] <- genes
}

#list.save(reactome_list, 'REACTOME sets list human ENSEMBL.rds')

#------------------------------------------------------------------------------ KEGG PATHWAYS-------------------------------------------------------------------------------------------
kegg <- read.table("./input_data/gene_sets/c2.cp.kegg.v7.4.entrez.txt", header=F, sep = '\t') #downloaded from GSEA
kegg_melt <- na.omit(reshape2::melt(kegg))

annot_human <-getBM(attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                    filters = "entrezgene_id", values = unique(kegg_melt$value), mart = human)
kegg_melt <- merge(kegg_melt, annot_human, by.x="value", by.y="entrezgene_id")

#transform list to HUMAN ENSEMBL
kegg_list <- vector("list")
for(i in unique(kegg_melt$V1)){
  genes <- kegg_melt[kegg_melt$V1==i,]$ensembl_gene_id
  kegg_list[[i]] <- genes
}

#list.save(kegg_list, 'kegg sets list human ENSEMBL.rds')

#------------------------------------------------------------------------------ biocarta PATHWAYS-------------------------------------------------------------------------------------------
biocarta <- read.table("./input_data/gene_sets/c2.cp.biocarta.v7.4.entrez.txt", header=F, sep = '\t') #downloaded from GSEA
biocarta_melt <- na.omit(reshape2::melt(biocarta))

annot_human <-getBM(attributes = c("ensembl_gene_id","entrezgene_id", "external_gene_name"),
                    filters = "entrezgene_id", values = unique(biocarta_melt$value), mart = human)
biocarta_melt <- merge(biocarta_melt, annot_human, by.x="value", by.y="entrezgene_id")

#transform list to HUMAN ENSEMBL
biocarta_list <- vector("list")
for(i in unique(biocarta_melt$V1)){
  genes <- biocarta_melt[biocarta_melt$V1==i,]$ensembl_gene_id
  biocarta_list[[i]] <- genes
}

#list.save(biocarta_list, 'biocarta sets list human ENSEMBL.rds')

#------------------------------------------------------------------------------4. COMBINE and SELECT GENE SETS---------------------------------------------------------------------------------
gene.sets <- c("REACTOME_INTERFERON_GAMMA_SIGNALING", "REACTOME_REGULATION_OF_LIPID_METABOLISM_BY_PPARALPHA","REACTOME_LINOLEIC_ACID_LA_METABOLISM", "REACTOME_METABOLISM_OF_STEROIDS", "REACTOME_METABOLISM_OF_STEROID_HORMONES",
               "REACTOME_PHOSPHOLIPID_METABOLISM", "REACTOME_GLYCOSPHINGOLIPID_METABOLISM","REACTOME_FATTY_ACID_METABOLISM",
               "REACTOME_METABOLISM_OF_LIPIDS", "REACTOME_ABC_TRANSPORTERS_IN_LIPID_HOMEOSTASIS", "REACTOME_INTERLEUKIN_21_SIGNALING","Immune cells", "Mature enterocytes", 
               "Transit amplifying cells", "Duodenal transporters",
               "KEGG_PPAR_SIGNALING_PATHWAY", "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY", "KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION",
               "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION", "BIOCARTA_NO2IL12_PATHWAY")

list_of_interest <- c(mglist2, transporters_list, reactome_list, kegg_list, biocarta_list)


list_of_interest <- list_of_interest[gene.sets]
list.save(list_of_interest, './input_data/gene_sets/list_of_interest_human_ENSEMBL.rds')



#------------------------------------------------------------------------------5. OUR CUSTOM GENE SETS-----------------------------------------------------------------------------------------
CD <- readRDS('./input_data/Celiac biopsies DEG results.rds')
ORG <- readRDS('./input_data/hDuo ZED1277 DEG results.rds')
a <- unique(CD[[1]][CD[[1]]$Gene.type == "DE",]$ID)
b <- unique(ORG[[2]][ORG[[2]]$Gene.type == "DE",]$ID)
epithelial_IFNg_responce <- b
epithelial_IFNg_responce_CeD <- intersect(a, b)
reactome_IFNg <- reactome_list[["REACTOME_INTERFERON_GAMMA_SIGNALING"]]
IFNg_pathways <- list("epithelial_IFNg_responce" = epithelial_IFNg_responce, "epithelial_IFNg_responce_CeD" = epithelial_IFNg_responce_CeD, "REACTOME_INTERFERON_GAMMA_SIGNALING" = reactome_IFNg)
list.save(IFNg_pathways, './input_data/gene_sets/IFNg_list_of_interest_human_ENSEMBL.rds')
