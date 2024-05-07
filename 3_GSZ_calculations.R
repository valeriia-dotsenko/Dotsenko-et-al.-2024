library(mGSZ) #package version 1.0 Not available from cran anymore
source("./mGSZ_function.R") # from http://ekhidna.biocenter.helsinki.fi/downloads/pashupati/pre_print.pdf
source("./sampleGSZ_function.R")

library(rlist)
library(dplyr)
library(numform)

CD <-readRDS("./input_data/Celiac biopsies DESeqDataSet object.rds")
CDres <- readRDS('./input_data/Celiac biopsies DEG results.rds')
metadata <-colData(CD)

#open gene sets lists
mglist <- readRDS("./input_data/gene_sets/list_of_interest_human_ENSEMBL.rds")
mglist_2 <- readRDS("./input_data/gene_sets/IFNg_list_of_interest_human_ENSEMBL.rds")
mglist_2 <- mglist_2[-3]
mglist <- c(mglist, mglist_2)

#transform counts
normalised_counts <- counts(CD, normalized=TRUE) 
tran_data<- log2(normalised_counts+1)


#calculate Sample GSZ
SampleGSZ <- Sample.GSZ(gene.set=mglist,ex.data=tran_data)


#Gene set analysis based on Gene Set Z-scoring function and asymptotic p-value
#calculate GSZs for groups
groups <- data.frame(group1=c("GFD","PGCp","PGCd"),group2=c("PGCd","GFD","PGCp"))

#DO NOT RE-RUN
#groupGSZ <- data.frame() #method contains permutations, so each time performed gives slightly different p-values 
#for(i in 1:nrow(groups)){
#  groups2 <- metadata[metadata$label3 %in% c(groups[i,1],groups[i,2]),]$label3
#  tran_data2 <- tran_data[,metadata[metadata$label3 %in% c(groups[i,1],groups[i,2]),]$ID]
#  set.seed(50) 
#  mGSZ.obj2 <- mGSZ(tran_data2, mglist, groups2, perm.number= 100, min.cl.sz = 5) #perm.number = 100, min.cl.sz = 5)
  
#  group1 <- c(rep(groups[i,1],length(mglist)))
#  group2 <- c(rep(groups[i,2],length(mglist)))
#  comparison <- paste(group1,"-", group2, sep="")
#  tmp <- cbind(mGSZ.obj2$mGSZ, group1, group2, comparison)
#  groupGSZ <- rbind(groupGSZ, tmp)
#}

#groupGSZ %>% 
#  mutate(p_new = ifelse(`pvalue` > 0.01, c(paste("italic('P')~`=", f_num(`pvalue`,2), "`")), `pvalue`))%>% 
#  mutate(p_new = ifelse(`pvalue` < 0.01, c(paste("italic('P')~`=", f_num(`pvalue`,3), "`")), p_new)) %>%
#  mutate(p_new = ifelse(`pvalue` < 0.001, c(paste("italic('P')~`", "<.001", "`")),p_new))->groupGSZ

#Save data
CDres[["SampleGSZ"]] <- SampleGSZ
CDres[["groupGSZ"]] <- groupGSZ

list.save(CDres, './input_data/Celiac biopsies DEG results.rds')
openxlsx::write.xlsx(CDres[8:9], file = "./output_data/GSZ_results_CeD_biopsies.xlsx")

