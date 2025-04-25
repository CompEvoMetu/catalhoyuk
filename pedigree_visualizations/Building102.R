library(kinship2)
bul102 <- pedigree(
  id = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
  dadid = c(0,0,0,2,2,0,2,0,0,3,3,6,8,0,9,14,3,0,18), 
  momid = c(0,0,0,1,1,0,1,0,0,4,4,5,7,0,10,15,4,0,7),
  sex = c(2,1,1,2,2,1,2,1,1,2,2,1,2,1,2,1,1,1,2),
  affected = c(0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1))

mt=c("","","","","","","","","","H","H","H","H","","","H27d","-","","H")
y=c("","","","","","","","","","-","-","-","-","","","J2a1a2b","-","-")
names=c("","","","","","","","","","cch103","cch102","cch471","cch126","","","cch109","cch073","","cch119")
building=c("","","","","","","","","","Building 102","Building 102","Building 128","Building 128","","","Building 102","Building 102","","Building 128")
age=c("","","","","","","","","","old adult","young adult","adult","young adult","","","child","young adult","","child")

pdf("building102.pdf",width=10,height=10)
plot(bul102,id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch073-cch205 Third Degree
#cch102-cch205 Unrelated, Unrelated on X with 1100 SNPs
#cch103-cch205 Unrelated, Third Degree on X with 1804 SNPs
#cch106-cch205 Third Degree, Second Degree on X with 156 SNPs
#cch109-cch205 Unrelated, Second Degree on X with 458 SNPs
#cch119-cch205 Unrelated, Unrelated on X with 2321 SNPs
#cch126-cch205 Unrelated, Unrelated on X with 403 SNPs
#cch205-cch471 Unrelated, Second Degree on X with 931 SNPs

#cch073-cch106 
#cch102-cch106 
#cch103-cch106 Third Degree on Read2
#cch106-cch109 
#cch106-cch119 Third Degree on Read2, Second Degree on X with 165 SNPs
#cch106-cch126 
#cch106-cch205 Third Degree, Second Degree on X with 156 SNPs
#cch106-cch471 Unrelated

#cch073-cch119 Third Degree on NGSRelate
#cch102-cch119 Third Degree on NGSRelate, Third Degree on X with 1109 SNPs
#cch103-cch119 Third Degree, Third Degree on X with 1793 SNPs
#cch109-cch119 Unrelated, Second Degree on X with 461 SNPs
#cch119-cch126 Second Degree, First Degree on X with 406 SNPs
#cch119-cch471 Third Degree, Second Degree on X with 926 SNPs

#cch073-cch102 First or Second Degree
#cch073-cch103 First or Second Degree
#cch073-cch109 674 common SNPs
#cch073-cch126 316 common SNPs
#cch073-cch471 Third Degree or Unrelated

#cch102-cch103 First Degree, First Degree on X
#cch102-cch109 Third Degree, First Degree on X
#cch102-cch126 Third Degree, Third Degree on X with 190 SNPs
#cch102-cch471 Third Degree, Third  Degree on X

#cch103-cch109 Second Degree, Second Degree on X
#cch103-cch126 Third Degree, Third Degree on X
#cch103-cch471 Second Degree, First Degree on X

#cch109-cch126 Unrelated 
#cch109-cch471 Third Degree on Read2 only, Third Degree on X

#cch126-cch471 Third Degree, Identical on X with 162 SNPs
