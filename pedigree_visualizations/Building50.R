library(kinship2)
bul50 <- pedigree(
  id = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
  dadid = c(0,0,0,1,0,1,0,1,3,5,5,5,7), 
  momid = c(0,0,0,2,0,2,0,2,4,6,6,6,8),
  sex = c(1,2,1,2,1,2,1,2,1,1,2,2,1),
  affected = c(0,0,0,0,0,0,0,0,1,1,1,1,1))
mt = c("","","","","","","","","K1a4","-","K1a4","K1a+195T! ","-")
y = c("","","","","","","","","G2a2a1","C1a2","","","")
names= c("","","","","","","","","cch289","cch124","cth728","cth842","cch288")
building= c("","","","","","","","","building 50","building 50","building 50","building 50","-")
age= c("","","","","","","","","child","child","infant","infant","child")

pdf("building50.1.pdf",width=8,height=8)
plot(bul50,col= c("black","black","black","black","black","black","black","black","lightpink2","lightpink2","lightpink2","lightpink2","black"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch124-cch288 Third Degree on Read2
#cch124-cch289 Second Degree
#cch124-cth728 First Degree
#cch124-cth842 First Degree

#cch288-cch289 Third Degree on Read2, Third Degree on X with 387 snps
#cch288-cth728 
#cch288-cth842 Third Degree on Read2,

#cch289-cth728 Third Degree, First Degree on X with 212 snps
#cch289-cth842 Third Degree, Unrelated on X with 859 snps

#cth728-cth842 First Degree

bul50.2 <- pedigree(
  id = c(1,2,3,4,5,6),
  dadid = c(0,0,2,2,0,5), 
  momid = c(0,0,1,1,0,4),
  sex = c(2,1,1,2,1,1),
  affected = c(0,0,1,0,0,1))
mt = c("","","-","","","HV")
y = c("","","G2a2a1","","","C1a2b")
age= c("","","middle adult","","","neonate")
names= c("","","cch223","","","cch294")
building= c("","","building 50","","","building 50")

pdf("building50.2.pdf",width=8,height=8)
plot(bul50.2,col=c("black","black","darkseagreen1","black","black","darkseagreen1"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch223-cch294 Second Degree, First Degree on X with 920 snps
