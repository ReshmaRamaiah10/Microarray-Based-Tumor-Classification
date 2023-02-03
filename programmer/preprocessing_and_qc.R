#Author- Aneeq Husain
#Description- This code imports the data, normalises it and performs quality control step. It also plots median histograms and a PCA plot. 

#load libraries
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
library(readr)
library(magrittr)

#Set file path
filepath <- '/projectnb/bf528/users/tinman_2022/project_1/samples/'

#import the data
affy_batch <- affy :: ReadAffy(celfile.path = filepath)

#Normalize data
norm_data <- rma(affy_batch)

#Computing Relative Log Expression (RLE) and NUSE scores-
rle_data <- fitPLM(affy_batch, normalize = TRUE, background = TRUE)

#Find medians
med_NUSE <- NUSE(rle_data, type= 'stats')  #finds the median NUSE for samples 
med_RLE <-RLE(rle_data, type = "stats") #find median RLE for samples

#plotting histograms
hist(med_RLE['median',], xlab = 'Median RLE Scores', ylab = 'Number of samples', main = 'Histogram of median NUSE Scores' ,border ='black', col = 'coral2')
hist(med_NUSE['median',], xlab = 'Median NUSE Scores', ylab = 'Number of samples', main = 'Histogram of median RLE Scores' ,border ='black', col = 'cadetblue2')
#Correcting for batch effects
anno_path = "/project/bf528/project_1/doc/proj_metadata.csv"

#import the required batch and model data
anno <- read_csv(anno_path, show_col_types = FALSE) %>% 
  transmute(combatbatch = normalizationcombatbatch, combatmod = normalizationcombatmod)
batch = anno$combatbatch
model = model.matrix(~combatmod, data = anno)

#perform batch correction
result <- ComBat(dat=exprs(affy_batch), batch = batch, mod = model)

#exporting the csv file
write.csv(result, file = "expression_data.csv")

#Performing PCA Analysis-
pca <- result %>%
  t() %>%
    scale() %>%
      t() %>%
        prcomp(scale = FALSE, center = FALSE)

#extracting plotting data from class prcomp
PC1 <- pca$rotation[,1]
PC2 <- pca$rotation[,2]
summ <- summary(pca)
per1 <- paste('PC1-', round(summ$importance[2,1] *100,2), '%')
per2 <- paste('PC2-', round(summ$importance[2,2] *100,2), '%')

#plotting PCA
plot(PC1, PC2, main = 'PCA plot of normalized data', xlab= per1, ylab =per2)
