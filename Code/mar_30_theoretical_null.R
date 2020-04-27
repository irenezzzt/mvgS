# perfect theoretical null make sure iSKAT is fine


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

options(warn = -1)

dum<-function(i){
  temp<-cbind(ifelse(Z[,i]==1,1,0),ifelse(Z[,i]==2,1,0) )
  colnames(temp)<-c(file.path(colnames(Z)[i],paste0(1)), file.path(colnames(Z)[i],paste0(2)))
  return(temp)
}
gS_test <-function(Z,Y,resid="absolute", type = "f"){
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

gene<-read.csv("~/Desktop/MgS/gene.csv",header=T)[,-1]
coef<-read.table("~/Desktop/MgS/coef.txt")
beta1<-coef$beta1; beta3<-coef$beta3

N.in<-10^5
n=dim(gene)[2]
p=11
alpha=0.05


# additive coding
t1e.sim<-function(xx){
  set.seed(xx*17+13)
  E <- matrix(rbinom(n,size=1,prob=0.3))
  select.null.snp<-sample(seq(1,dim(gene)[1],1), size=p)
  geno.new<-as.data.frame(t(gene[c(select.null.snp),]))
  #colnames(geno.new)<-c("G1", "G2", "G3","G4","G5","G6","G7","G8","G9","G10","G11")
  Z<-as.matrix(geno.new, ncol=p)
  
  #Y<-Z%*%beta1+0.015*E+E*(Z%*%beta3)+matrix(rnorm(n,0,0.27^2), ncol=1)
  #mean(Y);sd(Y)
  Y<-rnorm(n, 0.002, 0.073)
  # delete all indiv=0 SNPs
  #Z<-Z[,colMeans(Z)!=0]
  
  #perm.dat<-data.frame(Y_perm=Y, E_perm=E.obs)
  #perm.dat<-perm.dat[sample(nrow(perm.dat)), ]
  
  p1<-gS_test_ad(Z=Z, Y=Y,resid="absolute", type = "f")
  p2<-gS_test_ad(Z=Z, Y=Y,resid="absolute", type = "skat")
  p3<-gS_test_ad(Z=Z, Y=Y,resid="absolute", type = "skato")
  
  p4<-gS_test_ad(Z=Z, Y=Y, resid="square", type="f")
  p5<-gS_test_ad(Z=Z, Y=Y, resid="square", type="skat")
  p6<-gS_test_ad(Z=Z, Y=Y, resid="square", type="skato")
  
  return(c(p1,p2,p3,p4,p5,p6))
}

s<-Sys.time()
t1e.res<-matrix(unlist(lapply(X=1:10^4, FUN=t1e.sim)),ncol=6,byrow=T)
time.cost<-Sys.time()-s


sum(complete.cases(t1e.res))#99315/100000 non-singular
colnames(t1e.res)<-c("abs+F", "abs+SKAT","abs+SKATO", "sq+F", "sq+SKAT", "sq+SKATO")
colMeans(t1e.res<0.05,na.rm=T)

