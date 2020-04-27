# what happened to SKAT under perfect null?

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
res<-NULL

#check the null dist by p val
for(xx in 2496:10000){
  print(xx)
  set.seed(xx*17+3)
  E <- matrix(rbinom(n,size=1,prob=0.3)); E.obs=E
  select.null.snp<-sample(seq(1,dim(gene)[1],1), size=p)
  geno.new<-as.data.frame(t(gene[c(select.null.snp),]))
  colnames(geno.new)<-c("G1", "G2", "G3","G4","G5","G6","G7","G8","G9","G10","G11")
  Z<-as.matrix(geno.new, ncol=p)
  #Y<-Z%*%beta1+0.015*E+E*(Z%*%beta3)+matrix(rnorm(n,0,2*0.17^2), ncol=1)
  Y0<-rnorm(n, mean= 1, sd= 1)
  # delete all indiv=0 SNPs
  Z<-Z[,colMeans(Z)!=0]
  if(is.singular.matrix(t(Z)%*%Z)){res<-res}else{
  perm.dat<-data.frame(Y=Y0,E=E)
  perm.dat<-perm.dat[sample(nrow(perm.dat)),]
  
  Y_perm<-perm.dat$Y
  data<-as.data.frame(Z,Y_perm)
  
  lm1 <- rq( Y_perm ~ Z, data=data, tau=0.5)
  lm2 <- rq( Y0 ~ Z, tau=0.5 )
  data$d1 <- abs(lm1$residuals) # folded normal under the null
  data$d2 <- (lm1$residuals)^2 # chi sq under the null
  data$d3 <- rnorm(n, mean=mean(data$d1), sd=sd(data$d1))
  data$d4 <-abs(lm2$residuals)
  data$d5 <- (lm2$residuals)^2
  ##### folded normal #####
  # F 
  fit0<-lm(d1~1, data = data)
  fit <-lm(d1~Z, data = data)
  p1<-anova(fit0,fit)[2,6]
  # SKAT
  obj<-SKAT_Null_Model(data$d1~1, out_type="C")
  p2<-SKAT(as.matrix(Z), obj)$p.value
  # SKATO
  p3<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value

  #### square #####
  # F 
  fit0<-lm(d2~1, data = data)
  fit <-lm(d2~Z, data = data)
  p4<-anova(fit0,fit)[2,6]
  # SKAT
  obj<-SKAT_Null_Model(data$d2~1, out_type="C")
  p5<-SKAT(as.matrix(Z), obj)$p.value
  # SKATO
  p6<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
  
  #### normal #####
  # F 
  fit0<-lm(d3~1, data = data)
  fit <-lm(d3~Z, data = data)
  p7<-anova(fit0,fit)[2,6]
  # SKAT
  obj<-SKAT_Null_Model(data$d3~1, out_type="C")
  p8<-SKAT(as.matrix(Z), obj)$p.value
  # SKATO
  p9<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
  
  res<-rbind(res, c(p1,p2,p3,p4,p5,p6,p7,p8,p9) )
  }
}

#res<-matrix(unlist(lapply(X=c(1:10000), FUN=sim.null)), ncol=9, byrow=T)
colMeans(res<0.05)    
colnames(res)<-c("perm F", "perm SKAT", "perm SKATO", 
                 "sq F", "sq SKAT", "sq SKATO",
                 "norm F", "norm SKAT", "norm SKATO")

res.bur<-NULL

# burden in second stage?
for(xx in 9869:10000){
  print(xx)
  set.seed(xx*17+3)
  E <- matrix(rbinom(n,size=1,prob=0.3)); E.obs=E
  select.null.snp<-sample(seq(1,dim(gene)[1],1), size=p)
  geno.new<-as.data.frame(t(gene[c(select.null.snp),]))
  colnames(geno.new)<-c("G1", "G2", "G3","G4","G5","G6","G7","G8","G9","G10","G11")
  Z<-as.matrix(geno.new, ncol=p)
  Y<-Z%*%beta1+0.015*E+E*(Z%*%beta3)+matrix(rnorm(n,0,2*0.17^2), ncol=1)
  Y0<-rnorm(n, mean= 1, sd= 1)
  # delete all indiv=0 SNPs
  Z<-Z[,colMeans(Z)!=0]
  
  if(is.singular.matrix(t(Z)%*%Z)){res.bur<-res.bur}else{
    perm.dat<-data.frame(Y=Y,E=E)
    perm.dat<-perm.dat[sample(nrow(perm.dat)),]
    
    Y_perm<-perm.dat$Y
    data<-as.data.frame(Z,Y_perm)
    lm1 <- rq( Y_perm ~ Z, data=data, tau=0.5)
    lm2 <- rq( Y0 ~ Z, tau=0.5 )
    # all kinds of di ------
    data$d1 <- abs(lm1$residuals) # folded normal under the empirical null
    data$d2 <- (lm1$residuals)^2 # chi sq under the empirical null
    data$d3 <- rnorm(n, mean=mean(data$d1), sd=sd(data$d1)) # normal di
    data$d4 <- abs(data$d3) # folded normal di
    data$d5 <-abs(lm2$residuals) #folded normal under theoretical null
    data$d6 <- (lm2$residuals)^2 #chi sq under theoretical null

    
    # F  --------
    fit0<-lm(d1~1, data = data)
    fit <-lm(d1~Z, data = data)
    p11<-anova(fit0,fit)[2,6]
    fit0<-lm(d2~1, data = data)
    fit <-lm(d2~Z, data = data)
    p12<-anova(fit0,fit)[2,6]
    fit0<-lm(d3~1, data = data)
    fit <-lm(d3~Z, data = data)
    p13<-anova(fit0,fit)[2,6]
    fit0<-lm(d4~1, data = data)
    fit <-lm(d4~Z, data = data)
    p14<-anova(fit0,fit)[2,6]
    fit0<-lm(d5~1, data = data)
    fit <-lm(d5~Z, data = data)
    p15<-anova(fit0,fit)[2,6]
    fit0<-lm(d6~1, data = data)
    fit <-lm(d6~Z, data = data)
    p16<-anova(fit0,fit)[2,6]
   
    # SKAT n SKATO-------

    obj<-SKAT_Null_Model(data$d1~1, out_type="C")
    p21<-SKAT(as.matrix(Z), obj)$p.value
    p31<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
    
    obj<-SKAT_Null_Model(data$d2~1, out_type="C")
    p22<-SKAT(as.matrix(Z), obj)$p.value
    p32<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
    
    obj<-SKAT_Null_Model(data$d3~1, out_type="C")
    p23<-SKAT(as.matrix(Z), obj)$p.value
    p33<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
    
    obj<-SKAT_Null_Model(data$d4~1, out_type="C")
    p24<-SKAT(as.matrix(Z), obj)$p.value
    p34<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
    
    obj<-SKAT_Null_Model(data$d5~1, out_type="C")
    p25<-SKAT(as.matrix(Z), obj)$p.value
    p35<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
    
    obj<-SKAT_Null_Model(data$d6~1, out_type="C")
    p26<-SKAT(as.matrix(Z), obj)$p.value
    p36<-SKAT(as.matrix(Z), obj, method="SKATO")$p.value
    # burden ---------
 
   if( ncol(as.matrix(Z, nrow=n))== 1 ){ data$supG= Z}else{rowSums(Z)}
    p41<-summary(lm(d1~supG, data=data))$coefficient[2,4]
    p42<-summary(lm(d2~supG, data=data))$coefficient[2,4]
    p43<-summary(lm(d3~supG, data=data))$coefficient[2,4]
    p44<-summary(lm(d4~supG, data=data))$coefficient[2,4]
    p45<-summary(lm(d5~supG, data=data))$coefficient[2,4]
    p46<-summary(lm(d6~supG, data=data))$coefficient[2,4]
    
    # res ------
    res.bur<-rbind(res.bur, c(p11,p12,p13,p14,p15,p16,
                              p21,p22,p23,p24,p25,p26,
                              p31,p32,p33,p34,p35,p36,
                              p41,p42,p43,p44,p45,p46) )
  }
}
write.csv(res.bur, "t1e.burden.csv")
xtable(matrix(colMeans(res.bur<0.05, na.rm=T),ncol=6,byrow=T)[,-4],digits=4)
Fres<-colMeans(res.bur[,1:6]<0.0, na.rm = T)
SKATres<-colMeans(res.bur[,7:12]<0.05, na.rm = T)
SKATOres<-colMeans(res.bur[,13:18]<0.05, na.rm = T)
Bres<-colMeans(res.bur[,19:24]<0.05, na.rm = T)

getwd()

