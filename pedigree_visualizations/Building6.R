library(kinship2)
bul6 <- pedigree(
  id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
  dadid = c(0,0,0,1,1,1,0,1,0,3,7,9,9,0,14,14), 
  momid = c(0,0,0,2,2,2,0,2,0,4,6,8,8,0,13,13),
  sex = c(1,2,1,2,2,2,1,2,1,2,1,1,2,1,1,2),
  affected = c(0,0,0,0,1,0,0,0,0,1,1,1,0,0,1,1))
mt = c("","","","","N1a1a1b","","","","","N1a1a1b","N1a1a1b","N1a1a1b","","","N1a1a1b","N1a1")
y = c("","","","","-","","","","","-","J2a1a2b","C1a2b","","","C1a2b","-")
names = c("","","","","cch140","","","","","cch161","cch139","cch130","","","cch158","cch176")
building = c("","","","","Building 6","","","","","Building 6","Building 6","Building 6","","","Building 6","Building 6")
age = c("","","","","infant","","","","","neonate","neonate","prenatal","","","prenatal","infant")

pdf("building6.pdf",width=10,height=10)
plot(bul6,col=c("black","black","black","black","darkseagreen1","black","black","black","black","darkseagreen1","darkseagreen1","darkseagreen1","black","black","darkseagreen1","darkseagreen1"),id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch130-cch139 Third Degree
#cch130-cch140 Third Degree, Second Degree on X with 328 SNPs
#cch130-cch158 Second Degree
#cch130-cch161 Third Degree
#cch130-cch176 Second Degree, Third Degree on X

#cch139-cch140 Second Degree
#cch139-cch158 Unrelated
#cch139-cch161 Third Degree 
#cch139-cch176 Third Degree

#cch140-cch158 Third Degree
#cch140-cch161 Second Degree
#cch140-cch176 Third Degree, Third Degree on X with 618 SNPs

#cch158-cch161 Unrelated
#cch158-cch176 First Degree

#cch161-cch176 Third Degree
