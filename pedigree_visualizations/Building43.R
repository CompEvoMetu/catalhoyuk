library(kinship2)
bul43 <- pedigree(
  id = c(1,2,3,4,5,6,7,8,9,10),
  dadid = c(0,0,1,0,4,0,1,0,7,7), 
  momid = c(0,0,2,0,3,0,2,0,6,8),
  sex = c(1,2,2,1,1,2,1,2,1,1),
  affected = c(0,0,0,0,1,0,0,0,1,1))

names = c("","","","","cch180","","","","cch153","cch254")
mt = c("","","","","T2e+152C!","","","","T2e ","T2+16189!")
y = c("","","","","G2a2a1","","","","J2a1a2b2a2a","J2a")
building=c("","","","","Building 43","","","","Building 43","Building 43")
age=c("","","","","prenatal","","","","prenatal","prenatal")

pdf("building43.pdf",width=10,height=10)
plot(bul43,id=paste(names,age,mt,y,building,sep="\n"))
dev.off()

#cch153-cch180 Third Degree, Unrelated on X with 397 SNPs
#cch153-cch254 Second Degree, Unrelated on X with 921 SNPs
#cch180-cch254 Third Degree

#cch153        cch180        cch254
#T2e           T2e+152C!     T2+16189!
#B43           B43           B43
#J2a1a2b2a2a   G2a2a1        J2a
#XY            XY            XY
#prenatal      prenatal      prenatal
