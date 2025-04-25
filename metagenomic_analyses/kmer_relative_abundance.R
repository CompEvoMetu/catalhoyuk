# Visualizing k-mer counts by PCA -----------------------------------------
# Loading Libraries -------------------------------------------------------
library(tidyverse)
library(readxl)
#library(vegan)
#library(plotly)
#library(ggpubr)
#library(svMisc)
#library(plyr)
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

size <- args[1]

`%nin%` = Negate(`%in%`)
groups <- read.table("catalhoyuk_kmer.txt", header = T,sep = "\t")

kmer_sizes <- c("1K", "10K", "50K", "100K", "500K", "1M")
kmer <- as.character(size)


# Loading absolute abundance table ----------------------------------------
header <- read.table("./header", header = T)
data <- read.table(file = paste0("./kmer_", kmer), header = F, sep = " ", col.names = colnames(header))
header <- as.character(read.table("./header", header = F)[1,])
colnames(data) <- header

# Remove cth (capture)
data <- data[,-c(which(colnames(data) %in% 
                         colnames(data)[grepl("cth|scp1|cp1", x = colnames(data), ignore.case = T)]))]

header <- colnames(data)

substring(header, 7, 7 ) <- "."

header_short <- c()
for(i in 1:length(header)){
  
  header_short[i] <- tolower(str_split(header[i], pattern = "_")[[1]][[1]])
  
}

header_short <- unique(sort(header_short))
header_short <- header_short[header_short != "kmer"]
header_short <- header_short[which(nchar(header_short) %in% 15)]
header_short <- substr(header_short, start = 1, 13)

temp <- matrix(ncol = length(header_short)+1, nrow = nrow(data))
temp[,1] <- data$kmer

colnames(temp) <- c("kmer", header_short)

i = 2
for(i in 1:length(header_short)){
  
  temp[,header_short[i]] <- 
    rowSums(data[,c(which(colnames(data) 
      %in% colnames(data)[grepl(header_short[i], 
       colnames(data), ignore.case = T)])), drop =F], )
}

data <- temp
highprop <- groups$short[groups$age %in% c("adult_highprop", "subadult_highprop", "adult_lowprop", "subadult_lowprop")]

cols <- which(colnames(data) %in% colnames(data)[grepl(paste0(highprop, collapse = "|")
                                                       , colnames(data))])
data <- data[,c(1, cols)]

kmers <- as.character(data[,1])

# Calculating relative abundance ------------------------------------------
data <- apply(data,2, as.numeric)
relative_adundance <- data

relative_adundance[,2:ncol(relative_adundance)] <- as.numeric(relative_adundance[,2:ncol(relative_adundance)])
str(relative_adundance)

for(j in 2:ncol(data)){
  
  sum <- sum(data[,j])
  print(j)
  
  for(i in 1:nrow(data)){
    relative_adundance[,j] <- data[,j] / sum
  }
  
}

relative_adundance[,1] <- kmers

temp <- t(relative_adundance)
colnames(temp) <- temp[1,]
temp <- temp[-1,]

temp_num <- matrix(as.numeric(temp), ncol = ncol(temp))
colnames(temp_num) <- colnames(temp)
row.names(temp_num) <- row.names(temp)


write.table(temp_num, file = paste0("./rel_abundance/relabundance_", kmer), row.names = T, 
            col.names = T, quote = F)

