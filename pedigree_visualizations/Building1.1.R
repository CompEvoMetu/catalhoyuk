library(kinship2)
bul1 <- pedigree(
  id = c(1,2,5,6,7,8,9,10,11,12,14,15,16,17,18,19),
  dadid = c(0,0,1,0,0,1,1,0,1,0,5,10,12,7,7,7), 
  momid = c(0,0,2,0,0,2,2,0,2,0,6,9,11,8,8,8),
  sex = c(1,2,1,2,1,2,2,1,2,1,1,1,1,1,2,1),
  affected = c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1))

names=c("","","","","","","","","","","cch288","cch155","cch152","cch131","cch331","cch3208")
mt=c("","","","","","","","","","","-","U3b2","U3b2","U3b2","U3b2","U3b")
y=c("","","","","","","","","","","-","C1a2b","J2a","T1a3","-","T1a3")
building=c("","","","","","","","","","","-","Building 1","Building 1","Building 1","Building 1","Building 1")
age=c("","","","","","","","","","","child","infant","infant","child","child","child")

pdf("building1.pdf",width=8,height=8)
plot(bul1,col=c(rep("black",13),"darkorchid4","darkorchid4","darkorchid4"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch131-cch152 Third Degree, First Degree on X with 269 SNPs
#cch131-cch155 Third Degree
#cch131-cch288 Third Degree
#cch131-cch3208 First Degree, Identical on X with 1555 SNPs
#cch131-cch331 First Degree

#cch152-cch155 Third Degree
#cch152-cch288 Third Degree
#cch152-cch3208 Third Degree, First Degree on X with 2225 SNPs
#cch152-cch331 Third Degree

#cch155-cch288 Third Degree
#cch155-cch3208 Third Degree, Identical on X with 506 SNPs
#cch155-cch331 Third Degree on Read2 with 2273 SNPs

#cch288-cch3208 Third Degree, Unrelated on X with 840 SNPs
#cch288-cch331 Third Degree

#cch3208-cch331 First Degree, First Degree on X with 858 SNPs

#cch170-cch288 Third Degree on X with 562 SNPs
#cch170-cch331 Third Degree on X with 500 SNPs
#cch155-cch170 Second Degree on X with 304 SNPs
