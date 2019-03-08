setwd("C:/Users/jianleih/Dropbox/GMM R code") ##set work directory
#install.packages("sandwich");                ## install these packages at the first time run
#install.packages("optimx");
#install.packages("expm");
#install.packages("numDeriv");
#install.packages("readr");
library(expm);
library(sandwich);
library(optimx);
library(numDeriv);
############Beginning of User Defined Settings#########
NOBS=924;                                     ###Number of observations in the input data set
N=20;                                         ###Number of parameters to be estimated (unrestricted), if there's no unrestricted model, define it same to restrited model
NR=10;                                        ###Number of parameters to be estimated (restricted)
M=20;                                         ###Number of Moment conditions 
NZ=0;                                         ###Number of observations that are needed to form the initial values of any lagged variables in the data set
COVEST="Quadratic Spectral";                  ###select kernel used to estimate the covariance matrix: Quadratic Spectral, Parzens, Bartlett, Truncated
AR="AR(1)";                                    ###bandwidth selection: AR(1) or ARMA(1)
DIFSTAT=1;                                  
#         Switch that determines whether the GMM difference 
#               statistic is calculated.
#
#               If DIFSTAT=0, GMM difference statistic not calculated, only one restricted model estimated
#               If DIFSTAT=1, it is calculated
#
#     Note: Setting DIFSTAT=1 insures that the same weighting matrix
#     is used for both the unrestricted and restricted estimations.
#     Holding the weighting matrix fixed is necessary in order for
#     the test to be valid. The moment conditions for the unrestricted
#     estimation go in the GMM and FUNC1 subroutines. The moment
#     conditions for the restricted estimation go in FUNC2.
if (DIFSTAT==0) N=NR;  
iniPU=rep(1,N);      ##set up initial parameter guess for unrestricted model. If no unrestricted model, same to restricted
iniPR=rep(1,NR);                       ##set up initial parameter guess for restricted model
#read input

library(readr);
data1=read_csv("example1.csv");
###example1: value-we return, equal-we return and 10 deciles return
###variables: vwretd: value-weighted return; ewretd: equal-weighted return; d1-d10: 10 portfolio return.
data2=read_csv("example2.csv"); ## risk-free rate
data1$rf=data2$bid/1200;
data1$ratiom=(1+data1$vwretd)/(1+data2$bid/1200);
########End of User Defined Settings###########################
###objective function2: restricted model with NR parameters
FUN2=function(PR){                                  ###PR:Parameters for restricted model, for example, 10 betas
  Z=matrix(0,nrow(data1),M);
  for ( i in 1:10){                                 ### moment conditions for restricted model
    Z[,i]=as.matrix(data1[,3+i]-data1$rf)-PR[i]*(data1$vwretd-data1$rf);}
  for ( i in 11:20){
    Z[,i]=Z[,i-10]*(data1$vwretd-data1$rf);}
  g=colMeans(Z);
  Q=0;
  for (i in 1:M){
    Q=Q+mean(Z[,i])^2*WEGR[i]
  }
  result=list("residual"=Z,"moments"=g,"obj"=Q);
  return(result);
}
###objective function for function2
OBJ2=function(PR){
  return(FUN2(PR)$obj)
}
###function1: unrestricted model with N parameters
FUN1=function(PU){                                  ###PU: Parameters for unrestricted model, for example, 10 alphas and 10 betas
  Z=matrix(0,nrow(data1),M);
  for ( i in 1:10){                                    ### moment conditions for unrestricted model
  Z[,i]=as.matrix(data1[,3+i]-data1$rf)-PU[i]-PU[i+10]*(data1$vwretd-data1$rf);} 
  for ( i in 11:20){
  Z[,i]=Z[,i-10]*(data1$vwretd-data1$rf);}
  g=colMeans(Z);
  Q=0;
  for (i in 1:M){
    Q=Q+mean(Z[,i])^2*WEG[i]
  }
  result=list("residual"=Z,"moments"=g,"obj"=Q);
  return(result);
}
###objective function for function1
OBJ1=function(PU){
  return(FUN1(PU)$obj)
}

if (DIFSTAT==0) {FUN1=FUN2; OBJ2=OBJ1;N=NR;}     ### if there is no unrestricted model to be estimate, set unrestricted model same to restricted model


###function to get weight matrix W
WEM=function(G1){                          
  #y=rep(1,nrow(G1));
  lmobj=lm(G1~1);                          ##construct a lm object
  HACm=kernHAC(lmobj,kernel=COVEST,approx=AR)*NOBS  ##Get HAC matrix with kernel in Andrews (1991)
  W=solve(HACm);                              ##weight matrix as inverse(S)
  return(W)
}

###function to get D matrix
Dmatrix=function(P,restrict){
  
  if (!restrict){
    D=matrix(0,M,N);
    for (i in 1:M){
        dfun=function(pu){
          return(FUN1(pu)$moments[i])
        }
        D[i,]=grad(dfun,P)
    }
  }
  else{
    D=matrix(0,M,NR);
    for (i in 1:M){
       dfun=function(pu){
        return(FUN2(pu)$moments[i])
      }
      D[i,]=grad(dfun,P)
    }
  }
  return(D)
}
###step 1 with weight matrix of I(M)
WEG=WEGR=rep(1,M);

optpu=optimx(iniPU,OBJ1,method="nlminb"); ##optimizer "nlminb
PU1=as.matrix(optpu[1:N]);
G1=FUN1(PU1)$residual;                            ##store all residuals
###step 2 calculate GMM covariance matrix 
WE=WEM(G1)
for (i in 1:M){
  WEG[i]=WE[i,i];
}
WEGR=WEG;
optpu2=optimx(iniPU,OBJ1,method="nlminb"); ## second step optimize for unrestricted model

optpr2=optimx(iniPR,OBJ2,method="nlminb"); ## second step optimize for restricted model
## store estimated coefficients and residuals
Coef1=as.matrix(optpu2[1:N]);
Z1=FUN1(Coef1)$residual;
Coef2=as.matrix(optpr2[1:NR]);
Z2=FUN2(Coef2)$residual;
## calculate new weight matrix
WE1=WEM(Z1);
WE2=WEM(Z2);
##Call Dmatrix to calculate D
DU=Dmatrix(Coef1,restrict=FALSE);
DR=Dmatrix(Coef2,restrict=TRUE);
##Unrestricted model standard errors
STDE=solve(t(DU) %*% WE1 %*% DU)/NOBS;
if (nrow(STDE)>1) STDEU=sqrtm(STDE);
if (nrow(STDE)==1) STDEU=sqrt(STDE);
STDEEU=diag(STDEU);
tu=Coef1/STDEEU;
row.names(tu)="t-stat"
##Restricted model standard errors
STDE=solve(t(DR) %*% WE2 %*% DR)/NOBS;
if (nrow(STDE)>1) STDER=sqrtm(STDE);
if (nrow(STDE)==1) STDER=sqrt(STDE);
STDEER=diag(STDER);
tr=Coef2/STDEER;
row.names(tr)="t-stat"
row.names(Coef1)=row.names(Coef2)="Coefficients";
outputu=rbind(Coef1,tu);
outputr=rbind(Coef2,tr);
###J stats
gu=colMeans(Z1);
gr=colMeans(Z2);
Ju=NOBS*t(gu)%*%WE1%*%gu;
Jr=NOBS*t(gr)%*%WE2%*%gr;
###writeoutput
if (DIFSTAT==1){
  tp=1-pt(tu,NOBS-M);
  outdata=cbind(t(Coef1),t(tu),t(tp));
  J=matrix("",nrow(outdata),1);
  J[1,1]=Ju;
  outdata=cbind(outdata,J);
  outdata=data.frame(outdata);
  names(outdata)[3:4]=c("P-value","J-stat");
  tp=1-pt(tr,NOBS-NR);
  tr2=c(tr,rep("",N-NR));
  tp=c(tp,rep("",N-NR));
  co2=c(Coef2,rep("",N-NR));
  J[1,1]=Jr;
  outdata=cbind(outdata,co2,tr2,tp,J);
  J[1,1]=1-pchisq(Jr,M-NR);
  outdata=cbind(outdata,J);
  names(outdata)[5:9]=c("Coefficients_restricted","t.stat","P-value","J-stat","P-value")
  write.csv(outdata,"GMMoutput.csv",row.names=FALSE)
}else
  {
  tp=1-pt(tr,NOBS-NR);
  outdata=cbind(t(Coef2),t(tr),t(tp));
  J=matrix("",nrow(outdata),1);
  J[1,1]=Jr;
  outdata=cbind(outdata,J);
  J[1,1]=1-pchisq(Jr,M-NR);
  outdata=cbind(outdata,J);
  outdata=data.frame(outdata);
  names(outdata)[3:5]=c("P-value","J-stat","P-value")
  write.csv(outdata,"GMMoutput.csv",row.names=FALSE)
}
###Outputs
if (DIFSTAT==1){
  cat("Unrestricted model fit Coefficents:\n",Coef1,"\n","Unrestricted model fit Coefficents t-stats:\n",tu,"\n",
      "Unrestricted model J-stats:\n",Ju,"\n","Restricted model fit Coefficents:\n",Coef2,"\n",
      "Restricted model fit Coefficents t-stats:\n",tr,"\n","Restricted model J-stats:\n",Jr,"\n");
} else
{
  cat("Restricted model fit Coefficents:\n",Coef2,"\n",
      "Restricted model fit Coefficents t-stats:\n",tr,"\n","Restricted model J-stats:\n",Jr,"\n");
}


