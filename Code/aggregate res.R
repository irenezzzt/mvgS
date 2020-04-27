files = list.files(path="/Users/tingzhang/Desktop/compute_canada/temp_res_sg/",pattern="*.txt")

res<-NULL
setwd("/Users/tingzhang/Desktop/compute_canada/temp_res_sg/")
read.delim("/Users/tingzhang/Desktop/compute_canada/temp_res_rg/2.txt", sep=",")[,2]
for (f in files){
  res<-rbind(res, read.delim(f, sep=",")[,2])
}
colnames(res)<-c("gS+abs+F", "gS+sq+F","gS+abs+SKAT", "gS+sq+SKAT","iSKAT","GESAT")
res.sg<-colMeans(res<0.05, na.rm=T)
xtable(t(as.matrix(res.sg)),digits=4)
write.csv(res.rg, "~/Desktop/compute_canada/res.rg.csv")
xtable(t(as.matrix(read.csv("~/Desktop/compute_canada/res.rg.csv"))))
library(xtable)




files = list.files(path="/Users/tingzhang/Desktop/compute_canada/temp_res_rg/",pattern="*.txt")

res<-NULL
getwd()
setwd("/Users/tingzhang/Desktop/compute_canada/temp_res_rg/")
read.delim("/Users/tingzhang/Desktop/compute_canada/temp_res_rg/2.txt", sep=",")[,2]
for (f in files){
  res<-rbind(res, read.delim(f, sep=",")[,2])
}
colnames(res)<-c("gS+abs+F", "gS+sq+F","gS+abs+SKAT", "gS+sq+SKAT","iSKAT","GESAT")
res.rg<-colMeans(res<0.05, na.rm=T)
xtable(as.data.frame(res.rg),digits=4)
write.csv(res.rg, "~/Desktop/compute_canada/res.rg.csv")






files = list.files(path="/Users/tingzhang/Desktop/compute_canada/temp_res_be/",pattern="*.txt")

res<-NULL
getwd()
setwd("/Users/tingzhang/Desktop/compute_canada/temp_res_be/")
read.delim("/Users/tingzhang/Desktop/compute_canada/temp_res_be/2.txt", sep=",")[,2]
for (f in files){
  res<-rbind(res, read.delim(f, sep=",")[,2])
}
colnames(res)<-c("gS+abs+F", "gS+sq+F","gS+abs+SKAT", "gS+sq+SKAT","iSKAT","GESAT")
res.be<-colMeans(res<0.05, na.rm=T)
xtable(as.data.frame(res.be),digits=4)
write.csv(res.rg, "~/Desktop/compute_canada/res.be.csv")




files = list.files(path="/Users/tingzhang/Desktop/compute_canada/temp_res_gt/",pattern="*.txt")

res<-NULL
getwd()
setwd("/Users/tingzhang/Desktop/compute_canada/temp_res_gt/")
read.delim("/Users/tingzhang/Desktop/compute_canada/temp_res_gt/2.txt", sep=",")[,2]
for (f in files){
  res<-rbind(res, read.delim(f, sep=",")[,2])
}
colnames(res)<-c("gS+abs+F", "gS+sq+F","gS+abs+SKAT", "gS+sq+SKAT","iSKAT","GESAT")
res.gt<-colMeans(res<0.05, na.rm=T)
xtable(as.data.frame(res.gt),digits=4)
write.csv(res.rg, "~/Desktop/compute_canada/res.gt.csv")

res4<-rbind(res.gt,res.be,res.rg, res.sg)
rownames(res4)<-c("genotypic", "beta2","ridge", "sigma" )



files = list.files(path="/Users/tingzhang/Desktop/e0res/",pattern="*.txt")

res<-NULL
setwd("/Users/tingzhang/Desktop/e0res")
read.delim("/Users/tingzhang/Desktop/e0res/2.txt", sep=",")[,2]
for (f in files){
  res<-rbind(res, read.delim(f, sep=",")[,2])
}
colnames(res)<-c("gS+abs+F", "gS+sq+F","gS+abs+SKAT", "gS+sq+SKAT","iSKAT","GESAT")

res.sg<-colMeans(res<0.05, na.rm=T)
xtable(t(as.matrix(res.sg)),digits=4)
write.csv(res.rg, "~/Desktop/compute_canada/res.rg.csv")
xtable(t(as.matrix(read.csv("~/Desktop/compute_canada/res.rg.csv"))))
library(xtable)



#-------------



colnames(res)<-c("gS_f", "gS_SKAT", "gS_SKATO","iSKAT","SKAT","GESAT")
#setting.c<-colMeans(res<0.05, na.rm=T)
#setting.d<-colMeans(res<0.05, na.rm=T)
#setting.a<-colMeans(res<0.05, na.rm=T)
#setting.e<-colMeans(res<0.05, na.rm=T)
#setting.g<-colMeans(res<0.05, na.rm=T)
#setting.g[1]/setting.g[6]
#setting.h<-colMeans(res<0.05, na.rm=T)
#setting.h[1]/setting.h[6]
setting.h1<-colMeans(res<0.05, na.rm=T)
#change sigma epsilon can get result better
 
sum(is.na(res[,]))
files[which(is.na(res[,4]))]
dim(res)
# iSKAT & GESAT returned 4 NA, # 175.txt, 465.txt, 469.txt, 545.txt 

sum(is.na(res[,6])) 
files[which(is.na(res[,6]))]

0.691/0.772

fin.res<-rbind(setting.a, setting.b, setting.c, setting.d,
               setting.e, setting.f, setting.g, setting.h)
fin.res

#random coefficient
beta1<-matrix(rnorm(p,mean=0.01,sd=0.02)*rbinom(p,size=1,p=0.6),nrow=p)
beta3<-matrix(rnorm(p,mean=0.1,sd=0.2)*rbinom(p,size=1,p=0.4),nrow=p)

temp<-data.frame(beta1, beta3)
write.csv()

xtable(res4,digits=3)





