# F, SKAT(GESAT), SKATO(iSKAT) compare
# Y0 ~ normal, chisq, folded normal

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

setwd("~/Desktop/mvgS")
options(warn = -1)
gene<-read.csv("gene.csv",header=T)[,-1]
coef<-read.table("coef.txt")
beta1<-coef$beta1; beta3<-coef$beta3

N.in<-5000
N.out<-1000
n=dim(gene)[2]
p=11
alpha=0.05

# real Gs -----------

set.seed(1127)
E <- matrix(rbinom(n,size=1,prob=0.3))
select.null.snp<-sample(seq(1,dim(gene)[1],1), size=p)
Z<- as.matrix(t(gene[c(select.null.snp),]))
summary(Z)
GE<- apply(X=Z, MARGIN=2, FUN=function(x){x*E})


sim1<-function(ind){
  set.seed(ind)
  Y0 <-rnorm(n=n)
  
  fit0 <- lm(Y0 ~ Z + E )
  fit1 <- lm(Y0 ~ Z + E + Z*E )
  p1<-anova(fit0, fit1)$`Pr(>F)`[2]

  obj<-SKAT_Null_Model(Y0 ~ Z + E, out_type="C")

  p2<-SKAT(GE, obj)$p.value
  p3<-SKAT(GE, obj, method="SKATO")$p.value
  
  p4<-GESAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
  p5<-iSKAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
  return(c(p1, p2, p3, p4, p5))
  }
temp1<-t(sapply(X=c(1:N.in), FUN=sim1))
colMeans(temp1<0.05)


sim2<-function(ind){
  set.seed(ind)
  Y0 <-rchisq(n=n, df=3)/2.4

  fit0 <- lm(Y0 ~ Z + E )
  fit1 <- lm(Y0 ~ Z + E + Z*E )
  p1<-anova(fit0, fit1)$`Pr(>F)`[2]
  
  obj<-SKAT_Null_Model(Y0 ~ Z + E, out_type="C")
  p2<-SKAT(GE, obj)$p.value
  p3<-SKAT(GE, obj, method="SKATO")$p.value
  
  p4<-GESAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
  p5<-iSKAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
  return(c(p1, p2, p3, p4, p5))
}
temp2<-t(sapply(X=c(1:N.in), FUN=sim2))
colMeans(temp2<0.05)


sim3<-function(ind){
  set.seed(ind)
  Y0 <-abs(rnorm(n=n))
  
  fit0 <- lm(Y0 ~ Z + E )
  fit1 <- lm(Y0 ~ Z + E + Z*E )
  p1<-anova(fit0, fit1)$`Pr(>F)`[2]
  
  obj<-SKAT_Null_Model(Y0 ~ Z + E, out_type="C")
  p2<-SKAT(GE, obj)$p.value
  p3<-SKAT(GE, obj, method="SKATO")$p.value
  
  p4<-GESAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9)$pvalue
  p5<-iSKAT(Z, matrix(Y0,ncol=1), matrix(E, ncol=1), missing_cutoff = 0.9, MAF_cutoff = 0.5)$pvalue
  return(c(p1, p2, p3, p4, p5))
}
temp3<-t(sapply(X=c(1:N.in), FUN=sim3))
colMeans(temp3<0.05)

save(temp1, temp2, temp2, file="realG.RData")
load("realG.RData")
matrix(colMeans(cbind(temp1,temp2,temp3)<0.01),byrow=F, ncol=3)
1.96*sqrt(0.01*0.99/5000)+0.01
# rare  maf 0.01-0.05 -----------

set.seed(2002)
E <- matrix(rbinom(n,size=1,prob=0.3))

maf<-runif(n=11,min=0.01,max=0.04)
print(maf)
Z<-matrix(unlist(lapply(X=maf,FUN=rbinom,n=n, size=2)), byrow=F, ncol=11)

GE<- apply(X=Z, MARGIN=2, FUN=function(x){x*E})

temp1<-t(sapply(X=c(1:N.in), FUN=sim1))
colMeans(temp1<0.05)
temp2<-t(sapply(X=c(1:N.in), FUN=sim2))
colMeans(temp2<0.05)
temp3<-t(sapply(X=c(1:N.in), FUN=sim3))
colMeans(temp3<0.05)
save(temp1, temp2, temp2, file="rareG.RData")
load("rareG.RData")
matrix(colMeans(cbind(temp1,temp2,temp3)<0.01),byrow=F, ncol=3)
# common  maf 0.1-0.2 ----------


set.seed(2020)
E <- matrix(rbinom(n,size=1,prob=0.3))
maf<-runif(n=11,min=0.05,max=0.2)
print(maf)
Z<-matrix(unlist(lapply(X=maf,FUN=rbinom,n=n, size=2)), byrow=F, ncol=11)
GE<- apply(X=Z, MARGIN=2, FUN=function(x){x*E})


temp1<-t(sapply(X=c(1:N.in), FUN=sim1))
colMeans(temp1<0.05)
temp2<-t(sapply(X=c(1:N.in), FUN=sim2))
colMeans(temp2<0.05)
temp3<-t(sapply(X=c(1:N.in), FUN=sim3))
colMeans(temp3<0.05)
save(temp1, temp2, temp2, file="commonG.RData")

load("commonG.RData")

matrix(colMeans(cbind(temp1,temp2,temp3)<0.05),byrow=F, ncol=3)



plot the observed quantile in the test set
