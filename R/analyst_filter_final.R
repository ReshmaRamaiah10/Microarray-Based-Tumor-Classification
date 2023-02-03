#load libs
library(tidyverse)
library(tidyr)
library(dplyr)
library(matrixStats)


#read file
data <- read.csv("example_intensity_data_RK_36_inc_no_parenth.csv")

###4.1: filter 1###

##Next 3 lines are things I tried that didn't work:

#data_mean<- transform(data, Mean=apply(data[,-c(1)],1, mean, na.rm = TRUE)) 
#sum(data_mean$Mean >3.9)
#mutate(data, Count = count(data[,-c(1)]) >3.9)

n <- log2(15)

data_count <- data %>%
    mutate(Count = rowSums(select(data, 'GSM972389':'GSM972521') > n)) %>%
    mutate(Frequency = (Count/35)*100)

data_filt_1 <- filter(data_count, Frequency>20)



###4.2: filter 2###
#(Var-medianVariance)^2/medianVariance
#((n???1)*xvar(P)/varmed)
#two sided

#variance by row i.e. probeset
data_filt_1_var <- data_filt_1 %>% 
  mutate(row_wise_var = rowVars(as.matrix(data_filt_1[,c(2:36)])))

#store median variance
m <- median(data_filt_1_var$row_wise_var)

#number of columns? i.e. length of data
c <- ncol(data_filt_1[,c(2:36)])

#example for filtering by desired rows in dataframe
df <- data_filt_1[,c(2:36)]

#calculate test-statistic for data_filt_1
data_filt_1_t <- data_filt_1_var %>%
  mutate(t = (c-1)*(row_wise_var/m))

#p value?
p <-0.01

#Confidence level?
CL <- 1-p
#calculate qchisq upper/lower
#chiupper = qchisq(1 - 0.01/2, c-1)
#chilower = qchisq((0.01)/2, c-1)
chi_L = qchisq((1 -CL)/2, c-1)
chi_U = qchisq((CL)/2, c-1, lower.tail = FALSE)

#filter probes passing filter 1 by qchisq lower|upper and t values
data_filt_1_2 <- filter(data_filt_1_t, t > chi_U | t<chi_L) 

###4.3: filter3###
#output of filter 2 is data_filt_1_2
data_filt_1_2_cv <- data_filt_1_2 %>%
  mutate(row_wise_sd = rowSds(as.matrix(data_filt_1_2[,c(2:36)]))) %>%
  mutate(row_wise_mean = rowMeans(as.matrix(data_filt_1_2[,c(2:36)]))) %>%
  mutate(CV = (row_wise_sd/row_wise_mean))
  

data_filt_1_2_3 <- filter(data_filt_1_2_cv, CV>0.186)
filtered_genes <- data_filt_1_2_3[,c(1:36)]

####4.4###
#write out csv of all passing genes, without analysis columns
write.csv(filtered_genes,"C:\\Users\\rkafrawi\\Desktop\\BF528-Project1-analyst\\analyst_4_5\\filtered_genes.csv")

### 4.5 ###

#variance by row i.e. probeset
data_var <- data %>% 
  mutate(row_wise_var = rowVars(as.matrix(data[,c(2:36)])))

#store median variance
m2 <- median(data_var$row_wise_var)

#number of columns? i.e. length of data
c2 <- ncol(data[,c(2:36)])

#calculate test-statistic for data_filt_1
data_t <- data_var %>%
  mutate(t = (c-1)*(row_wise_var/m))

#p value?
p <-0.01

#Confidence level?
CL <- 1-p
#calculate qchisq upper/lower
#chiupper = qchisq(1 - 0.01/2, c-1)
#chilower = qchisq((0.01)/2, c-1)
chi_L = qchisq((1 -CL)/2, c-1)
chi_U = qchisq((CL)/2, c-1, lower.tail = FALSE)

#filter probes passing filter 1 by qchisq lower|upper and t values
data_filt_2_only <- filter(data_t, t > chi_U | t<chi_L) 
data_filt_2_only_clean <- data_filt_2_only[, -c(37,38)]

#generate csv
write.csv(data_filt_2_only_clean,"C:\\Users\\rkafrawi\\Desktop\\BF528-Project1-analyst\\analyst_4_5\\data_filt_2_only_clean.csv")

