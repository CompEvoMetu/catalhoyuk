library(kinship2)
bul129 <- pedigree(
  id = c(1,2,3,4,5,6,7,8),
  dadid = c(0,0,2,0,3,3,0,7), 
  momid = c(0,0,1,0,4,4,0,6),
  sex = c(2,1,1,2,2,2,1,2),
  affected = c(1,0,0,0,1,0,0,1))
mt = c("N1b1a1","","","","T2","","","T2")
y =c("-","","","","-","","","-")
names=c("cch242","","","","cch463","","","cch346")
building= c("Building 129","","","","Building 129","","","Building 129")
age= c("old adult","","","","adult","","","infant")

pdf("building129.s1.pdf",width=10,height=10)
plot(bul129,col=c("pink4","black","black","black","pink3","black","black","black"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

bul129.2 <- pedigree(
  id = c(1,2,3,4,5,6,7,8),
  dadid = c(0,0,0,2,2,0,3,5), 
  momid = c(0,0,0,1,1,0,4,6),
  sex = c(2,1,1,2,1,2,2,2),
  affected = c(1,0,0,0,0,0,1,1))
mt = c("T2","","","","","","T2","N1b1a1")
y = c("-","","","","","","-","-")
names=c("cch463","","","","","","cch346","cch242")
building= c("Building 129","","","","","","Building 129","Building 129")
age= c("adult","","","","","","infant","old adult")

pdf("building129.s2.pdf",width=10,height=10)
plot(bul129.2,col=c("pink3","black","black","black","black","black","black","pink4"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch162-cch242 Third Degree, Third Degree on X with 969 SNPs
#cch162-cch346 Unrelated, Unrelated on X with 314 SNPs
#cch162-cch463 Third Degree, Second Degree on X with 1021 SNPs

#cch242-cch346 Third Degree, Third Degree on X with 376 SNPs
#cch242-cch463 Second Degree, First Degree on X with 1301 SNPs

#cch346-cch463 Second Degree, Second Degree on X with 433 SNPs

bul129 <- pedigree(
  id = c(3,4,5,6,7,10,11,12,13),
  dadid = c(0,0,0,3,3,6,6,0,11), 
  momid = c(0,0,0,4,4,5,5,0,12),
  sex = c(1,2,2,1,2,2,1,2,2),
  affected = c(0,0,1,0,1,1,0,0,1))
mt = c("","","T2","","K1a12a","T2","","","N1b1a2")
y = c("","","-","","-","-","","","-")
names=c("","","cch346","","cch162","cch463","","","cch242")
building= c("","","Building 129","","Building 57","Building 129","","","Building 129")
age= c("","","infant","","neonate","adult","","","old adult")

pdf("building129.pdf",width=10,height=10)
plot(bul129,col=c("black","black","black","black","black","pink3","black","black","pink4"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()
