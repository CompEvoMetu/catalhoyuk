#@author: goztag

#!/usr/local/bin/Rscript

arg = commandArgs(trailingOnly = TRUE)

vcf = arg[1] #Name of the vcf file that contains imputed genomes
outfile = arg[2] #Name of the output file
windowsize = as.numeric(arg[3]) #Window size

#Read vcf file
vcf = read.table(vcf, header = F, sep = '\t', comment.char = "#")

#Convert vcf file to genotype matrix
geno <- t(apply(vcf[, 10:ncol(vcf)], 1, function(x){
        as.numeric(c(substr(x, 1, 1), substr(x, 3, 3)))}))

#Get position list
pos <- as.numeric(vcf[, c(2)])

#Create function to calculate mean difference between haplotypes
find_summary_window = function(gtfile, snppos, windowsize){
summary_window <- t(sapply(seq(1, (length(pos) - windowsize), by = windowsize), function(i){

                  #Get the window
                  WINDOW = i:(i+windowsize-1)

                  #Create table of genotypes within the window
                  haplotypes <- geno[WINDOW,]

                  #Calculate haplotype mismatch
                  distances <- dist(t(haplotypes), method = "manhattan") 

                  #Calculate the proportion of 0's, ie same haplotype
                  prop0 <- mean(distances == 0)

                  #Calculate the mean difference between haplotypes
                  diff <- mean(distances[distances != 0])

                  #Calculate size of window in bp
                  pos_in_window = pos[WINDOW]
                  windowsize_bp <- c(max(pos_in_window) - min(pos_in_window))
                  c(prop0, diff, windowsize_bp)}))

                  #Save results to a table
                  colnames(summary_window) <- c('Proportion_same', 'Mean_difference', 'Size_in_bp')
                  summary_window$Proportion_difference = 1 - summary_window$Proportion_same
                  write.table(summary_window, outfile, quote = F, row.names = F, col.names = T)
}

#Start function
find_summary_window(gtfile, snppos, windowsize)
