# what happen to iSKAT under perfect null

library(devtools)
#install_github("dsoave/gJLS")
library(gJLS)
#install.packages("car")
library(car)
#install.packages("SKAT")
library(SKAT)
#install_github("lin-lab/iSKAT-GESAT")
library(iSKAT)
#install.packages("quantreg")
library(quantreg)
#install.packages("L1pack")
library(L1pack)
library(assertthat)
library(matrixcalc)
library(xtable)

setwd("~/Desktop/MgS")
options(warn = -1)
dum<-function(i){
  temp<-cbind(ifelse(Z[,i]==1,1,0),ifelse(Z[,i]==2,1,0) )
  colnames(temp)<-c(file.path(colnames(Z)[i],paste0(1)), file.path(colnames(Z)[i],paste0(2)))
  return(temp)
}
gS_test_gt <-function(Z,Y,resid="absolute", type = "f"){
  #first stage L1 regression to obtain residual
  
  #Y_norm<-qnorm( (rank(Y)-0.5)/(length(Y)) )
  #summary(Y_norm)
  #Z_cent<-apply(FUN=scale, center=T, X=Z, MARGIN=2) 
  
  
  #if 0 G left, quit 
  if(dim(as.matrix(Z))[2]<2){stop(return(rep(NA,6)))}
  
  #expand Z to dummy code 
  cat.Z<-dum(1:dim(as.matrix(Z))[2]) #summary(cat.Z)
  cat.Z<-cat.Z[,colMeans(cat.Z)!=0]
  
  # if singular Z, quit
  if(is.singular.matrix(t(cat.Z)%*%cat.Z)){stop(return(rep(NA,6)) )}
  
  data<-as.data.frame(cat.Z,Y)
  lm1 <- rq( Y ~ cat.Z, data=data, tau=0.5)
  
  # residual choice
  if (resid=="absolute"){ 
    data$d1 <- abs(lm1$residuals)
  }else if (resid=="square"){
    data$d1 <- lm1$residuals^2
  }else{ stop("residual type not recognized")}
  
  # second stage regress test variance heterogeneity
  if (type=="f"){
    fit0<-lm(d1~1, data = data)
    fit <-lm(d1~cat.Z, data = data)
    aovfit<-anova(fit0,fit)
    #gS_F<-aovfit[2,5];numDF<-aovfit[2,3];denDF<-aovfit[2,1];
    gS_p<-aovfit[2,6]
  } else if (type=="skat"){
    obj<-SKAT_Null_Model(data$d1~1, out_type="C")
    gS_p<-SKAT(as.matrix(cat.Z), obj)$p.value
  } else if(type=="skato"){
    obj<-SKAT_Null_Model(data$d1~1, out_type="C")
    gS_p<-SKAT(as.matrix(cat.Z), obj, method="SKATO")$p.value
  }else{
    stop("stage II method not recognized")
  }
  
  # return the p-value
  return(c(gS_p))
} 
gS_test_ad <-function(Z,Y,resid="absolute", type = "f"){
  #first stage L1 regression to obtain residual
  
  #Y_norm<-qnorm( (rank(Y)-0.5)/(length(Y)) )
  #summary(Y_norm)
  #Z_cent<-apply(FUN=scale, center=T, X=Z, MARGIN=2) 
  
  
  #if 0 G left, quit 
  #if(dim(as.matrix(Z))[2]<2){stop(return(rep(NA,6)))}
  #expand Z to dummy code 
  #cat.Z<-dum(1:dim(as.matrix(Z))[2]) #summary(cat.Z)
  #cat.Z<-cat.Z[,colMeans(cat.Z)!=0]
  
  
  # if singular Z, quit
  if(is.singular.matrix(t(Z)%*%Z)){stop(return(rep(NA,6)) )}
  
  data<-as.data.frame(Z,Y)
  lm1 <- rq( Y ~ Z, data=data, tau=0.5)
  
  # residual choice
  if (resid=="absolute"){ 
    data$d1 <- abs(lm1$residuals)
  }else if (resid=="square"){
    data$d1 <- lm1$residuals^2
  }else{ stop("residual type not recognized")}
  
  # second stage regress test variance heterogeneity
  if (type=="f"){
    fit0<-lm(d1~1, data = data)
    fit <-lm(d1~Z, data = data)
    aovfit<-anova(fit0,fit)
    #gS_F<-aovfit[2,5];numDF<-aovfit[2,3];denDF<-aovfit[2,1];
    gS_p<-aovfit[2,6]
  } else if (type=="skat"){
    obj<-SKAT_Null_Model(data$d1~1, out_type="C")
    gS_p<-SKAT(as.matrix(Z), obj)$p.value
  } else if(type=="skato"){
    obj<-SKAT_Null_Model(data$d1~1, out_type="C")
    gS_p<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
  }else{
    stop("stage II method not recognized")
  }
  
  # return the p-value
  return(c(gS_p))
} 

gene<-read.csv("~/Desktop/MgS/gene.csv",header=T)[,-1]
coef<-read.table("~/Desktop/MgS/coef.txt")
beta1<-coef$beta1; beta3<-coef$beta3
interact.ind<-which(beta3!=0)

N.in<-10^5
N.out<-1000
n=dim(gene)[2]
p=11
alpha=0.05
res.isk<-NULL
#check the null dist by p val

for(xx in 1618:10000){
  print(xx)
  set.seed(xx*17+3)
  E <- matrix(rbinom(n,size=1,prob=0.3)); E.obs=E
  select.null.snp<-sample(seq(1,dim(gene)[1],1), size=p)
  geno.new<-as.data.frame(t(gene[c(select.null.snp),]))
  colnames(geno.new)<-c("G1", "G2", "G3","G4","G5","G6","G7","G8","G9","G10","G11")
  Z<-as.matrix(geno.new, ncol=p)
  #Y<-Z%*%beta1+0.015*E+E*(Z%*%beta3)+matrix(rnorm(n,0,2*0.17^2), ncol=1)
  Y0<-rnorm(n, mean= 1, sd= 1)
  Y1<-rchisq(n, df=5, ncp=1)
  # delete all indiv=0 SNPs
  Z<-Z[,colMeans(Z)!=0]
 Z<-matrix(Z, nrow=n)
  if(is.singular.matrix(t(Z)%*%Z)){res.isk<-res.isk}else{
    perm.dat<-data.frame(Y=Y,E=E)
    perm.dat<-perm.dat[sample(nrow(perm.dat)),]
    
    Y_perm<-perm.dat$Y
    E_perm<-perm.dat$E
    data<-as.data.frame(Z,Y_perm, E_perm)
    p1<-GESAT(Z, matrix(Y_perm,ncol=1), matrix(E_perm, ncol=1), missing_cutoff = 0.9)$pvalue
    p2<-iSKAT(Z=Z,Y=matrix(Y_perm,ncol=1),E=matrix(E_perm,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    p3<-GESAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
    p4<-iSKAT(Z=Z,Y=matrix(Y0,ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    p5<-GESAT(Z, matrix(Y1,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
    p6<-iSKAT(Z=Z,Y=matrix(Y1,ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    
    res.isk<-rbind(res.isk, c(p1,p2,p3,p4,p5,p6) )
  }
}

#res<-matrix(unlist(lapply(X=c(1:10000), FUN=sim.null)), ncol=9, byrow=T)
colnames(res.isk)<-c("perm gesat", "perm isk", "norm gesat", "norm isk", "chisq gesat", "chisq isk")
colMeans(res.isk<0.05, na.rm=T)    

# 1. both sensitive to non-normality (score is derived from likelihood), perm = mixture of chi sq
# 2. iSKAT combine burden and skat, burden more affected by non-normality the approximation is not working(?)


# iSKAT under folded normal Y0

set.seed(1127)
E <- matrix(rbinom(n,size=1,prob=0.3)); E.obs=E
select.null.snp<-sample(seq(1,dim(gene)[1],1), size=p)
geno.new<-as.data.frame(t(gene[c(select.null.snp),]))
colnames(geno.new)<-c("G1", "G2", "G3","G4","G5","G6","G7","G8","G9","G10","G11")
Z<-as.matrix(geno.new, ncol=p)
Y<-Z%*%beta1+0.015*E+E*(Z%*%beta3)+matrix(rnorm(n,0,0.17^2), ncol=1)
Z<-Z[,colMeans(Z)!=0]
Z<-matrix(Z, nrow=n)
res.isk.2<-NULL
xx=10001
for(xx in 10001:1000){
  print(xx)
  set.seed(xx*17+3)
  Y0<-rnorm(n, mean= 0, sd= 4)
  Y1<-rchisq(n, df=5, ncp=1)
  Y2<-abs(Y0)

    perm.dat<-data.frame(Y=Y,E=E)
    perm.dat<-perm.dat[sample(nrow(perm.dat)),]
    Y_perm<-perm.dat$Y
    E_perm<-perm.dat$E
    
    data<-as.data.frame(Z,Y_perm, E_perm)
    p1<-GESAT(Z, matrix(Y_perm,ncol=1), matrix(E_perm, ncol=1), missing_cutoff = 0.9)$pvalue
    p2<-GESAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
    p3<-GESAT(Z, matrix(Y1,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
    p4<-GESAT(Z, matrix(Y2,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
   
    p5<-iSKAT(Z=Z,Y=matrix(Y_perm,ncol=1),E=matrix(E_perm,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    p6<-iSKAT(Z=Z,Y=matrix(Y0,ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    p7<-iSKAT(Z=Z,Y=matrix(Y1,ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    p8<-iSKAT(Z=Z,Y=matrix(Y2,ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
    
    res.isk.2<-rbind(res.isk.2, c(p1,p2,p3,p4,p5,p6,p7,p8) )

}
res.isk.2<-res.isk.2[1:10000,]
dim(res.isk.2)[]
write.csv(res.isk.2,"~/Desktop/res.isk3.csv")
colMeans(res.isk.2<0.05,na.rm=T)

hist(Y_perm, main="Histogram of permuted Y")
hist(Y0, main="Histogram of normal Y0")
hist(Y1, main="Histogram of chisq Y0")
hist(Y2,main="Histogram of folded normal Y0")

# iskat center and scaled Y inside, the pvalues are exactly the same.
iSKAT(Z=Z,Y=matrix(Y1,ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
iSKAT(Z=Z,Y=matrix(scale(Y1,center=T),ncol=1),E=matrix(E,ncol=1),missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue



