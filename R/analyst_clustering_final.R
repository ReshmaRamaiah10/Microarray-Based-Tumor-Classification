#load required libraries
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(stats)
library(genefilter)
library(matrixTests)

#read filtered gene list
filtered_genes <- read.csv("filtered_genes.csv")

###5.1###


fg_trans <- as.data.frame(t(as.matrix(filtered_genes)))

#hclust for filtered_genes
clusters <- hclust(dist(fg_trans[2:35,]), method = 'complete')
plot(clusters)


###5.2###

#cut cluster in half
clusterCut <- cutree(clusters, k=2) %>%
    #turn named vector into a tibble
    enframe() %>%
    #rename columns
    rename(Sample = name, Cluster = value)

head(clusterCut)


###5.3###

#extract headers of filtered_genes and store into df


fgm <- as.matrix(filtered_genes)

# heatmap() requires a R matrix, and cannot accept a tibble or a dataframe
marker_matrix <- as.matrix(
  dplyr::select(filtered_genes,c('GSM972389':'GSM972521'))
)

# rownames of the matrix become y labels
rownames(marker_matrix) <- filtered_genes$Genes

heatmap(marker_matrix)






#above was the base heatmap, I need a heatmap with colsidecolors
#color class on trans data then trans back? tidy? 4.4?

#preparing column colors
#pull columns and make matrix
filtered_genes_columns <- colnames(filtered_genes)
filtered_genes_columns_M <- as.matrix(filtered_genes_columns)
fgc <- filtered_genes_columns_M[-1,]
fgcm <- as.matrix(fgc)

#add column with class for each sample
fgcm_class <- cbind(fgcm, Class = c('C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C4'
                                    ,'C4'
                                    ,'C4'
                                    ,'C4'
                                    ,'C4'
                                    ,'C3'
                                    ,'C4'
                                    ,'C4'
                                    ,'C4'
                                    ,'C3'
                                    ,'C4'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C4'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C4'
                                    ,'C3'
                                    ,'C4'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C3'
                                    ,'C4')
                    )
##old stuff I tried, can disregard
#fg_cols <- filtered_genes_columns(c(-1))

#read in classifications
#class <- read.csv("proj_metadata.csv")

#class_2 <- select(class, c(SixSubtypesClassification, geo_accession))
#fgt_samplesonly <- select(fg_trans, c(Genes))

#this but with the output of fg_trans with class
condition_colors <-
  transmute(
    data.frame(fgcm_class),
    color=if_else(Class == "C3","red","blue")
  )

#reassociate Sample w color
cond_color_w_col<- cbind(condition_colors, Sample = fgcm[,1])
                    
heatmap(
  marker_matrix,
  ColSideColors=condition_colors$color
)

### 5.4 ###

#make df for cluster 1
cluster_2_samples = c('GSM972401',
                     'GSM972409',
                     'GSM972413',
                      'GSM972422',
                      'GSM972433',
                      'GSM972438',
                      'GSM972440',
                      'GSM972444',
                      'GSM972467',
                      'GSM972476',
                      'GSM972479')
fg_cluster_1_only<- dplyr::select(filtered_genes,-c(all_of(cluster_2_samples)))

#make df for cluster 2
fg_cluster_2_only <- dplyr::select(filtered_genes,c(all_of(cluster_2_samples)))

#Welch's t-test
t_test_output <- row_t_welch(as.matrix(fg_cluster_1_only[-c(1)]), as.matrix(fg_cluster_2_only[-c(1)]), alternative = "two.sided", mu = 0, conf.level = 0.95)

#add back probes
t_test_output_probes <- as_tibble(t_test_output) %>%
                            mutate(Gene = filtered_genes[c(1)])

#p.adjusted
pvalue <- t_test_output_probes$pvalue
welch_comp <- t_test_output_probes %>%
                  mutate(adj_p_value = p.adjust(pvalue, 
                  method = p.adjust.methods, n = length(pvalue)))
#trim welch comp df
welch_comp_trimmed <- select(welch_comp,c('Gene', 'statistic', 'pvalue', 'adj_p_value' ))

#check how many pass filter for pvalue <0.05
sub_p<- filter(welch_comp_trimmed, pvalue < 0.05)
sub_p.adj<- filter(welch_comp_trimmed, adj_p_value < 0.05)

#export welch_comp
welch_comp_trimmed_matrix <- as.matrix(welch_comp_trimmed)
typeof(welch_comp_trimmed_matrix)
write.csv(welch_comp_trimmed_matrix,"C:\\Users\\rkafrawi\\Desktop\\BF528-Project1-analyst\\analyst_4_5\\welch_comp_trimmed_matrix.csv", 
          row.names = TRUE)

### 5.6 ###
#for biologist

fgb_cluster_1_only<- dplyr::select(data_filt_2_only_clean,-c(all_of(cluster_2_samples)))

#make df for cluster 2
fgb_cluster_2_only <- dplyr::select(data_filt_2_only_clean,c(all_of(cluster_2_samples)))

#Welch's t-test
b_t_test_output <- row_t_welch(as.matrix(fgb_cluster_1_only[-c(1)]), as.matrix(fgb_cluster_2_only[-c(1)]), alternative = "two.sided", mu = 0, conf.level = 0.95)

#add back probes
b_t_test_output_probes <- as_tibble(b_t_test_output) %>%
  mutate(Gene = data_filt_2_only_clean[c(1)])

#p.adjusted
pvalue2 <- b_t_test_output_probes$pvalue
biologist_welch_comp <- b_t_test_output_probes %>%
  mutate(adj_p_value = p.adjust(pvalue2, 
                                method = p.adjust.methods, n = length(pvalue2)))

#trim welch comp df
biologist_welch_comp_trimmed <- select(biologist_welch_comp,c('Gene', 'statistic', 'pvalue', 'adj_p_value' ))

#export welch_comp
biologist_welch_comp_trimmed_matrix <- as.matrix(biologist_welch_comp_trimmed)
typeof(biologist_welch_comp_trimmed_matrix)
write.csv(biologist_welch_comp_trimmed_matrix,"C:\\Users\\rkafrawi\\Desktop\\BF528-Project1-analyst\\analyst_4_5\\biologist_welch_comp_trimmed_matrix.csv", 
          row.names = TRUE)