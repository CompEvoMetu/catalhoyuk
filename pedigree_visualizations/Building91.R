library(kinship2)
bul91 <- pedigree(
  id = c(1,2,5,6,7,9,10,11,12,15,16,17,19),
  dadid = c(0,0,1,1,0,0,6,6,0,0,9,11,15), 
  momid = c(0,0,2,2,0,0,7,7,0,0,10,12,16),
  sex = c(1,2,1,1,2,1,2,1,2,1,2,1,2),
  affected = c(0,0,1,0,0,0,0,0,0,0,1,1,1))
names = c("","","cch254","","","","","","","","cch531","cch288","cth747")
mt =  c("","","T2+16189!","","","","","","","","T2c","H","T2c1+146T!")
building =  c("","","Building 43","","","","","","","","Building 91","-","Building 91")
y =  c("","","J2a","","","","","","","","-","-","-")
age = c("","","neonate","","","","","","","","middle adult","-","infant")
pdf("building91.pdf",width=8,height=8)
plot(bul91,id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch254-cch288 Third Degree on Read2
#cch254-cch531 Third Degree
#cch254-cth747 Third Degree on Read2
#cch288-cch531 Third Degree on Read2
#cch288-cth747 Third Degree on Read2
#cch531-cth747 First Degree

#cch129-cch254 Third Degree on Read2, Unrelated on X with 259 SNPs
#cch129-cch288 Third Degree on Read2, Unrelated on X with 176 SNPs
#cch129-cch339 Unrelated, Third Degree on X with 350 SNPs
#cch129-cch531 Third Degree on Read2, 
#cch129-cch747 Unrelated, First Degree on X with 418 SNPs

#cch254-cch339 Unrelated
#cch288-cch339 Third Degree on Read2
#cch339-cch531 Third Degree on Read2
#cch339-cth747 Unrelated. Third Degree on X with 242 SNPs
