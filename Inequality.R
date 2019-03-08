setwd("C:/Users/jianleih/Dropbox/GMM R code")
#install.packages("sandwich");                ## install these packages at the first time run
#install.packages("optimx");
#install.packages("expm");
#install.packages("numDeriv");
#install.packages("readr");
#install.packages("MASS");
library(expm)
library(sandwich)
library(optimx)
library(numDeriv)
library(MASS)
testdata=read.csv("MIM.csv")
restrictions=3; #number of restrictions
###function to get inversed vcov matrix
Iomega=function(G1){                          
  lmobj=lm(G1~1);                          ##construct a lm object
  HACm=kernHAC(lmobj,kernel="Quadratic Spectral",approx="AR(1)")  ##Get HAC matrix with kernel in Andrews (1991)
  W=solve(HACm)/nrow(G1);                              ##inversed Omega
  return(W)
}
###objective function
OBJ=function(theta,theta_hat){
  a=t(theta_hat-theta)%*%IOmega_hat%*%(theta_hat-theta)
  return(as.numeric(a))
}

theta_hat=colMeans(testdata) ##get mean values
IOmega_hat=Iomega(as.matrix(testdata))
inipa=rep(0.1,restrictions);  #set initial parameters
solution=optimx(inipa,OBJ,lower=rep(0,restrictions),method="nlminb",theta_hat=theta_hat)
theta_bar=as.numeric(solution[1:restrictions])
W=nrow(testdata)*OBJ(theta_bar,theta_hat)
##simulation for probablity matrix
set.seed(5256);
rep=1000                                     ##number of repetitions
vcov=solve(IOmega_hat);
sample=mvrnorm(n=rep,rep(0,restrictions),vcov);
count=rep(0,restrictions+1);
for (i in 1:rep){
  theta_star=as.matrix(sample[i,])
  solution=optimx(inipa,OBJ,lower=rep(0,restrictions),method="nlminb",theta_hat=theta_star)
  number=1;
  for (k in 1:restrictions) {
    number=number+(solution[k]>0)
  }
  count[number]=count[number]+1
}
pr=count/rep; #probablity vector
pvalue=0;
for (i in 2:length(pr)){
  pvalue=pvalue+(1-pchisq(W,(restrictions+1-i)))*pr[i]
}
cat("W-stat:",W,"\n","P-value:",pvalue)
write.csv(cbind(W,pvalue),"ineqoutput.csv",row.names=FALSE)
