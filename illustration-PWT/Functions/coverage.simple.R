##################################################################################
## The following codes will report the results for estimating the aggregate MPI
##
##   x1 is n times np matrix for inputs in the first period
##   y1 is n times nq matrix for outputs in the first period
##   x2 is n times np matrix for inputs in the second period
##   y2 is n times nq matrix for oututs in the second period
##
## see the following papers for more details:
## Inference for DEA Estimators of Malmquist Productivity Indices: 
## An Overview, Further Improvements and a Guide for Practitioners
## 
## The programming codes used in this paper involve 
## some earlier codes from Paul Wilson
## All rights reserved. 
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors.
##################################################################################

coverage.simple <- function(x1,y1,x2,y2) {
  L=100
  np=ncol(x1)
  nq=ncol(y1)
  n=nrow(x1)
  na=floor(n/2)
  nb=n-na
  kappa=2/(np+nq+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  # evaluate cononical efficiency, Output Orientation
  d11= 1/rDEA::dea(X=x1,Y=y1,XREF=x1,YREF=y1, 
                   model="output", RTS="constant")$thetaOpt
  d12= 1/rDEA::dea(X=x1,Y=y1,XREF=x2,YREF=y2, 
                   model="output", RTS="constant")$thetaOpt
  d21= 1/rDEA::dea(X=x2,Y=y2,XREF=x1,YREF=y1, 
                   model="output", RTS="constant")$thetaOpt
  d22= 1/rDEA::dea(X=x2,Y=y2,XREF=x2,YREF=y2, 
                   model="output", RTS="constant")$thetaOpt
  #
  Mi=0.5*(log(d11)+log(d12)-log(d21)-log(d22))
  M=mean(Mi)
  sig=sd(Mi)
  ########################
  # data sharpening method
  tau0=n^(-kappa)
  #
  epison11=rep(1,length(d11))
  ii=which(1/d11>=1-tau0)
  epison11[ii]=runif(length(ii),1-tau0,1)
  d11.nsz=d11/epison11
  #
  y1.nsz=y1*epison11
  #
  d12.nsz= 1/rDEA::dea(X=x1,Y=y1.nsz,XREF=x2,YREF=y2, 
                   model="output", RTS="constant")$thetaOpt
  epison12=d12/d12.nsz
  #
  epison22=rep(1,length(d22))
  ii=which(1/d22>=1-tau0)
  epison22[ii]=runif(length(ii),1-tau0,1)
  d22.nsz=d22/epison22
  #
  y2.nsz=y2*epison22
  #
  #
  d21.nsz= 1/rDEA::dea(X=x2,Y=y2.nsz,XREF=x1,YREF=y1, 
                       model="output", RTS="constant")$thetaOpt
  epison21=d21/d21.nsz
  #
  Mi.nsz=0.5*(log(d11.nsz)+log(d12.nsz)-log(d21.nsz)-log(d22.nsz))
  M.nsz=mean(Mi.nsz)
  sig.nsz=sd(Mi.nsz)
  ######################
  # compute bias corrections via generalized jackknife:
  # involving some earlier codes from Paul Wilson
  tbar=0
  tbar.nsz=0
  tbar11=rep(0,n)
  tbar12=rep(0,n)
  tbar21=rep(0,n)
  tbar22=rep(0,n)
  ind=c(1:n)
  for (j in 1:L) {
    if (j==1) {
      ind1=c(1:n)
      x1.b=x1
      y1.b=y1
      x2.b=x2
      y2.b=y2
      epison11.b=epison11
      epison12.b=epison12
      epison21.b=epison21
      epison22.b=epison22
    } else {
      ind1=sample(ind,size=n)
      x1.b[1:n,]=x1[ind1,]
      y1.b[1:n,]=y1[ind1,]
      x2.b[1:n,]=x2[ind1,]
      y2.b[1:n,]=y2[ind1,]
      epison11.b=epison11[ind1]
      epison12.b=epison12[ind1]
      epison21.b=epison21[ind1]
      epison22.b=epison22[ind1]
    }
    #
    x1a=matrix(x1.b[1:na,],ncol=np)
    y1a=matrix(y1.b[1:na,],ncol=nq)
    x1b=matrix(x1.b[(na+1):n,],ncol=np)
    y1b=matrix(y1.b[(na+1):n,],ncol=nq)
    #
    x2a=matrix(x2.b[1:na,],ncol=np)
    y2a=matrix(y2.b[1:na,],ncol=nq)
    x2b=matrix(x2.b[(na+1):n,],ncol=np)
    y2b=matrix(y2.b[(na+1):n,],ncol=nq)
    #
    d11a= 1/rDEA::dea(X=x1a,Y=y1a,XREF=x1a,YREF=y1a, 
                    model="output", RTS="constant")$thetaOpt
    d12a= 1/rDEA::dea(X=x1a,Y=y1a,XREF=x2a,YREF=y2a, 
                      model="output", RTS="constant")$thetaOpt
    d21a= 1/rDEA::dea(X=x2a,Y=y2a,XREF=x1a,YREF=y1a, 
                      model="output", RTS="constant")$thetaOpt
    d22a= 1/rDEA::dea(X=x2a,Y=y2a,XREF=x2a,YREF=y2a, 
                      model="output", RTS="constant")$thetaOpt
    #
    Mia=0.5*(log(d11a)+log(d12a)-log(d21a)-log(d22a))
    Ma=mean(Mia)
    #
    d11a.nsz=d11a/epison11.b[1:na]
    d12a.nsz=d12a/epison12.b[1:na]
    d21a.nsz=d21a/epison21.b[1:na]
    d22a.nsz=d22a/epison22.b[1:na]
    Mia.nsz=0.5*(log(d11a.nsz)+log(d12a.nsz)-log(d21a.nsz)-log(d22a.nsz))
    Ma.nsz=mean(Mia.nsz)
    #
    d11b= 1/rDEA::dea(X=x1b,Y=y1b,XREF=x1b,YREF=y1b, 
                      model="output", RTS="constant")$thetaOpt
    d12b= 1/rDEA::dea(X=x1b,Y=y1b,XREF=x2b,YREF=y2b, 
                      model="output", RTS="constant")$thetaOpt
    d21b= 1/rDEA::dea(X=x2b,Y=y2b,XREF=x1b,YREF=y1b, 
                      model="output", RTS="constant")$thetaOpt
    d22b= 1/rDEA::dea(X=x2b,Y=y2b,XREF=x2b,YREF=y2b, 
                      model="output", RTS="constant")$thetaOpt
    #
    Mib=0.5*(log(d11b)+log(d12b)-log(d21b)-log(d22b))
    Mb=mean(Mib)
    #
    d11b.nsz=d11b/epison11.b[(na+1):n]
    d12b.nsz=d12b/epison12.b[(na+1):n]
    d21b.nsz=d21b/epison21.b[(na+1):n]
    d22b.nsz=d22b/epison22.b[(na+1):n]
    Mib.nsz=0.5*(log(d11b.nsz)+log(d12b.nsz)-log(d21b.nsz)-log(d22b.nsz))
    Mb.nsz=mean(Mib.nsz)
    #
    tbar=tbar+(Ma+Mb)/2-M
    tbar.nsz=tbar.nsz+(Ma.nsz+Mb.nsz)/2-M.nsz
    #
    tbar11[ind1[1:na]]=tbar11[ind1[1:na]] +
      d11a - d11[ind1[1:na]]
    tbar11[ind1[(na+1):n]]=tbar11[ind1[(na+1):n]] +
      d11b - d11[ind1[(na+1):n]]
    
    tbar12[ind1[1:na]]=tbar12[ind1[1:na]] +
      d12a - d12[ind1[1:na]]
    tbar12[ind1[(na+1):n]]=tbar12[ind1[(na+1):n]] +
      d12b - d12[ind1[(na+1):n]]
    
    tbar21[ind1[1:na]]=tbar21[ind1[1:na]] +
      d21a - d21[ind1[1:na]]
    tbar21[ind1[(na+1):n]]=tbar21[ind1[(na+1):n]] +
      d21b - d21[ind1[(na+1):n]]
    
    tbar22[ind1[1:na]]=tbar22[ind1[1:na]] +
      d22a - d22[ind1[1:na]]
    tbar22[ind1[(na+1):n]]=tbar22[ind1[(na+1):n]] +
      d22b - d22[ind1[(na+1):n]]
  }
  # compute the bias
  M.bias=(1/L)*bc.fac*tbar
  tbar11=(1/L)*bc.fac*tbar11
  tbar12=(1/L)*bc.fac*tbar12
  tbar21=(1/L)*bc.fac*tbar21
  tbar22=(1/L)*bc.fac*tbar22
  #
  M.bias.nsz=(1/L)*bc.fac*tbar.nsz
  tbar11.nsz=tbar11/epison11
  tbar12.nsz=tbar12/epison12
  tbar21.nsz=tbar21/epison21
  tbar22.nsz=tbar22/epison22
  ##############################
  # using szz (2023,JPA) method
  Mi.szz=0.5*(log(d11-tbar11)+log(d12-tbar12)-log(d21-tbar21)-log(d22-tbar22))
  sig.szz=sd(Mi.szz)
  #
  Mi.nsz.szz=0.5*(log(d11.nsz-tbar11.nsz)+log(d12.nsz-tbar12.nsz)-log(d21.nsz-tbar21.nsz)-log(d22.nsz-tbar22.nsz))
  sig.nsz.szz=sd(Mi.nsz.szz)
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  ts=sig/sqrt(n)
  ts.szz=sig.szz/sqrt(n)
  ts.nsz=sig.nsz/sqrt(n)
  ts.nsz.szz=sig.nsz.szz/sqrt(n)
  bounds.st=matrix((M-ts*crit),nrow=3,ncol=2)
  if (np+nq<4) {
    bounds.ksw=matrix((M-M.bias-ts*crit),nrow=3,ncol=2)
    bounds.szz=matrix((M-M.bias-ts.szz*crit),nrow=3,ncol=2)
    bounds.nsz=matrix((M.nsz-M.bias.nsz-ts.nsz*crit),nrow=3,ncol=2)
    bounds.nsz.szz=matrix((M.nsz-M.bias.nsz-ts.nsz.szz*crit),nrow=3,ncol=2)
    # make a list of results to return to calling routine and then quit:
    res=list(bounds.st=bounds.st,
             bounds.ksw=bounds.ksw,
             bounds.szz=bounds.szz,
             bounds.nsz=bounds.nsz,
             bounds.nsz.szz=bounds.nsz.szz,
             sig=c(sig,sig.szz,sig.nsz,sig.nsz.szz),
             estimate=c(M,M.nsz,M.bias,M.bias.nsz,M-M.bias,M.nsz-M.bias.nsz)
             )
  } else {
    M.nk=mean(Mi[1:nk])
    M.nsz.nk=mean(Mi.nsz[1:nk])
    ts.nk=sig/sqrt(nk)
    ts.szz.nk=sig.szz/sqrt(nk)
    ts.nsz.nk=sig.nsz/sqrt(nk)
    ts.nsz.szz.nk=sig.nsz.szz/sqrt(nk)
    bounds.ksw=matrix((M.nk-M.bias-ts.nk*crit),nrow=3,ncol=2)
    bounds.szz=matrix((M.nk-M.bias-ts.szz.nk*crit),nrow=3,ncol=2)
    bounds.nsz=matrix((M.nsz.nk-M.bias.nsz-ts.nsz.nk*crit),nrow=3,ncol=2)
    bounds.nsz.szz=matrix((M.nsz.nk-M.bias.nsz-ts.nsz.szz.nk*crit),nrow=3,ncol=2)
    res=list(bounds.st=bounds.st,
             bounds.ksw=bounds.ksw,
             bounds.szz=bounds.szz,
             bounds.nsz=bounds.nsz,
             bounds.nsz.szz=bounds.nsz.szz,
             sig=c(sig,sig.szz,sig.nsz,sig.nsz.szz),
             estimate=c(M.nk,M.nsz.nk,M.bias,M.bias.nsz,M.nk-M.bias,M.nsz.nk-M.bias.nsz)
    )
  }
  return(res)
}