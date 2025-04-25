library(kinship2)
bul17 <- pedigree(
  id = c(1,2,3,4,5),
  dadid = c(0,0,1,1,1), 
  momid = c(0,0,2,2,2),
  sex = c(1,2,1,2,2),
  affected = c(0,0,1,1,1))
mt = c("","","H","H","H")
y = c("","","","-","-")
names = c("","","cch172","cch226","cch165")
building = c("","","Building 17","Building 17","Building 17")
age = c("","","neonate","child","infant")
pdf("building17.pdf",width=8,height=8)
plot(bul17,col=c("black","black","rosybrown3","rosybrown3","rosybrown3"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch165-cch172 First Degree
#cch165-cch226 First Degree, First Degree on X with 368 SNPs
#cch172-cch226 First Degree

bul17 <- pedigree(
  id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
  dadid = c(0,0,0,0,0,1,1,3,3,0,3,0,5,10,11,7,7,7), 
  momid = c(0,0,0,0,0,2,2,4,4,0,4,0,6,9,12,8,8,8),
  sex = c(1,2,1,2,1,2,1,2,2,1,1,2,2,2,2,1,2,2),
  affected = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1))
mt = c("","","","","","","","","","","","","K1a4","H","K1a+195T!","H","H","H")
y = c("","","","","","","","","","","","","T1a3","-","-","","-","-")
names = c("","","","","","","","","","","","","cch115","cch126","cth842","cch172","cch226","cch165")
building = c("","","","","","","","","","","","","Building 161","Building 128","Building 50","Building 17","Building 17","Building 17")
age = c("","","","","","","","","","","","","young adult","young adult","infant","neonate","child","infant")
pdf("building17.pdf",width=8,height=8)
plot(bul17,id=paste(names,age,mt,y,building,sep="\n"))
dev.off()
