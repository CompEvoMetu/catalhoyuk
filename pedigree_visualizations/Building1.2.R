library(kinship2)
bul1.2 <- pedigree(
  id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21),
  dadid = c(0,0,1,1,0,0,4,0,4,0,4,17,6,8,10,12,0,0,17,0,20), 
  momid = c(0,0,2,2,0,0,5,0,5,0,5,18,7,7,9,11,0,0,18,0,19),
  sex = c(1,2,1,1,2,1,2,1,2,1,2,1,1,1,1,2,1,2,2,1,2),
  affected = c(0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,1))

names=c("","","cch3418","","","","","","","","","cch107","cch251","cch117","cch339","cch170","","","","","cch080")
mt=c("","","T2e","","","","","","","","","T2","X2B","X2B","X2B","X2B","","","","","T")
y=c("","","H2a","","","","","","","","","-","","-","J2a1a","","","","","","")
building=c("","","Building 1","","","","","","","","","Building 1","Building 1","Building 1","Building 1","Building 1","","","","","Building 1")
age=c("","","child_infant","","","","","","","","","young adult","prenatal","child","infant","neonate","","","","","child")

pdf("building1.2.pdf",width=8,height=8)
plot(bul1.2,id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch080-cch107 Second Degree with 1183 SNPs 
#cch080-cch117 Third Degree or Unrelated wtih 2321 SNPs
#cch080-cch170 Third Degree, Third Degree on X with 778 SNPs
#cch080-cch251 Unrelated, Unrelated on X with 1649 SNPs
#cch080-cch339 Unrelated
#cch080-cch3418 Unrelated, Third Degree on X with 959 SNPs

#cch107-cch117 not enough SNPs for anything
#cch107-cch170 First Degree
#cch107-cch251 Unrelated, Unrelated on X with 255 SNPs
#cch107-cch339 Unrelated
#cch107-cch3418 Unrelated 

#cch117-cch170 Third Degree, Unrelated on X with 206 SNPs
#cch117-cch251 Second Degree 510 SNPs Unrelated on X with
#cch117-cch339 Third Degree
#cch117-cch3418 Third Degree, Second Degree on X with 285 SNPs

#cch170-cch251 Third Degree, Second Degree on X with 111267 SNPs
#cch170-cch339 Third Degree, Second Degree on X with 1006 SNPs
#cch170-cch3418 Third Degree, Unrelated on X with 6514 SNPs

#cch251-cch339 Third Degree or Unrelated, Unrelated on X with 2159 SNPs
#cch251-cch3418 Third Degree, Third Degree on X with 14200 SNPs

#cch339-cch3418 Third Degree, Unrelated on X with 1240 SNPs
