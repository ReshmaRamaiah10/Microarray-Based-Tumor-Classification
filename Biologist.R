# set the directory 
setwd("/projectnb/bf528/users/tinman_2022/project_1/")

# Installing the required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db"),force=TRUE)
BiocManager::install("GSEABase")

# Loading libraries
library(hgu133plus2.db)
library(tidyverse)
library(GSEABase)
library(magrittr)

# reading the file and creating gene symbols
t_summary <- read.csv('/projectnb2/bf528/users/tinman_2022/project_1/5.6_biologist_welch_comp_trimmed_matrix.csv')
gene_symbol <- AnnotationDbi::select(hgu133plus2.db, keys = t_summary$Probe_ID, columns = ("SYMBOL"), keytype = "PROBEID")

# All symbols in dataframe:
t_summary$symbol <- gene_symbol[!duplicated(gene_symbol$PROBEID), "SYMBOL"]
head(t_summary)
cat("Total gene symbols: ",length(t_summary$symbol))

# Removing duplicate symbols
symbols_in_df <- unique(t_summary$symbol)
cat("Unique length: ",length(symbols_in_df), "\n")

# some probeset IDs map to the same gene symbol
# choosing the probset IDs with higher significance
# if adj_p_values are same then chose the ID with max value
t_summary_no_dup <- t_summary %>% 
  group_by(symbol) %>%
  filter(adj_p_value == min(adj_p_value),symbol!='NA') %>%
  filter(Probe_ID == max(Probe_ID))

# top 1000 and 10 up/down-regulated genes
top_bottom_1000 <- rbind(t_summary_no_dup[order(t_summary_no_dup$statistic, decreasing = T),][1:1000,],t_summary_no_dup[order(t_summary_no_dup$statistic, decreasing = F),][1:1000,])
top_bottom_10 <- rbind(t_summary_no_dup[order(t_summary_no_dup$statistic, decreasing = T),][1:10,],t_summary_no_dup[order(t_summary_no_dup$statistic, decreasing = F),][1:10,])
knitr::kable(top_bottom_10, digits=100)

# top 10 upregulated and downregulated genes stored in a csv file
write.csv(top_bottom_10, file="top_and_bottom_10.csv", row.names=FALSE)

# importing KEGG, GO, and Hallmark gene sets from MSigDB
kegg_set <- getGmt(con="c2.cp.kegg.v7.5.1.symbols.gmt.txt")
go_set <- getGmt(con="c5.go.v7.5.1.symbols.gmt.txt")
hm_set <- getGmt(con="h.all.v7.5.1.symbols.gmt.txt")

# number of gene sets in each datbases
cat('Number of genesets in KEGG:',length(kegg_set))
cat('Number of genesets in GO:',length(go_set))
cat('Number of genesets in Halmark:',length(hm_set))

# creating the contigency table and conducting a fisher test
contigency_table <- function(geneset,geneset_collection){
  deg <- geneset[(order(abs(geneset$statistic),decreasing = TRUE)[1:1000]),'symbol']
  non_deg <- geneset[(-order(abs(geneset$statistic),decreasing = TRUE)[1:1000]),'symbol']
  ingene_deg <- sum(deg %in% geneset_collection)
  ingene_non_deg <- sum(non_deg %in% geneset_collection)
  not_ingene_deg <- sum(! deg %in% geneset_collection)
  not_ingene_non_deg <- sum(! non_deg %in% geneset_collection)
  
  res <- fisher.test(matrix(c(ingene_deg,not_ingene_deg,ingene_non_deg,not_ingene_non_deg),nrow=2))
  return(c(res$estimate,res$p.value))
}

# adjusting the p-values for multiple hypotheses using the Benjamini-Hochberg (FDR) 
adjust_p <- function(query, genesetcollection){
  res <- sapply(seq_len(length(genesetcollection)), function(x){
    p <- contigency_table(query, geneIds(genesetcollection[[x]]))
    tmp_res <- c(setName(genesetcollection[[x]]), p)
    names(tmp_res) <- c("Term", "Odds_ratio", "p_value")
    return(tmp_res)
  })
  res <- t(res) %>% data.frame()
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res <- res[order(res$p_adj, decreasing = F),]
  return(res)
}

# upregulated and down regulated genesets in KEGG, GO and Halmark databases
kegg_up_res <- adjust_p(t_summary_no_dup[t_summary_no_dup$statistic > 0,],kegg_set@.Data)
kegg_down_res <- adjust_p(t_summary_no_dup[t_summary_no_dup$statistic < 0,],kegg_set@.Data)
go_up_res <- adjust_p(t_summary_no_dup[t_summary_no_dup$statistic > 0,],go_set@.Data)
go_down_res <- adjust_p(t_summary_no_dup[t_summary_no_dup$statistic < 0,],go_set@.Data)
halmark_up_res <- adjust_p(t_summary_no_dup[t_summary_no_dup$statistic > 0,],hm_set@.Data)
halmark_down_res <- adjust_p(t_summary_no_dup[t_summary_no_dup$statistic < 0,],hm_set@.Data)

# number of significantly enriched gene sets at adjusted p<0.005
kegg_significant <- length((filter(rbind(kegg_up_res,kegg_down_res),p_adj<0.05))$p_adj)
go_significant <- ((filter(rbind(go_up_res,go_down_res),p_adj<0.05))$p_adj)
hm_significant <- ((filter(rbind(halmark_up_res,hallmark_down_res),p_adj<0.05))$p_adj)

# Top 3 genesets for each datbase
kegg_top3 <- head(kegg_up_res, order_by=p_adj, n=3)
go_top3 <- head(go_up_res, order_by=p_adj, n=3)
hm_top3 <- head(hallmark_up_res, order_by=p_adj, n=3)
top3s <- rbind(kegg_top3, go_top3, hm_top3)
write.csv(top3s, file="geneset_top3_results.csv")