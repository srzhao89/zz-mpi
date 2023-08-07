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


coverage.agg <- function(x1,y1,x2,y2) {
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
  U1=d21*as.vector(y2)
  U2=d22*as.vector(y2)
  U3=d11*as.vector(y1)
  U4=d12*as.vector(y1)
  U5=as.vector(y2)
  U6=as.vector(y1)
  mu1=mean(U1)
  mu2=mean(U2)
  mu3=mean(U3)
  mu4=mean(U4)
  mu5=mean(U5)
  mu6=mean(U6)
  xi=0.5*(log(mu3)+log(mu4)-log(mu1)-log(mu2))+log(mu5)-log(mu6)
  #
  Ti=cbind(U1,U2,U3,U4,U5,U6)
  Sigma=cov(Ti)
  g1=-1/(2*mu1)
  g2=-1/(2*mu2)
  g3=1/(2*mu3)
  g4=1/(2*mu4)
  g5=1/mu5
  g6=-1/mu6
  g=c(g1,g2,g3,g4,g5,g6)
  V.Sigma=t(g)%*%Sigma%*%g
  sig=sqrt(V.Sigma)
  sig=as.vector(sig)
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
  U1.nsz=d21.nsz*as.vector(y2)
  U2.nsz=d22.nsz*as.vector(y2)
  U3.nsz=d11.nsz*as.vector(y1)
  U4.nsz=d12.nsz*as.vector(y1)
  mu1.nsz=mean(U1.nsz)
  mu2.nsz=mean(U2.nsz)
  mu3.nsz=mean(U3.nsz)
  mu4.nsz=mean(U4.nsz)
  xi.nsz=0.5*(log(mu3.nsz)+log(mu4.nsz)-log(mu1.nsz)-log(mu2.nsz))+log(mu5)-log(mu6)
  #
  Ti.nsz=cbind(U1.nsz,U2.nsz,U3.nsz,U4.nsz,U5,U6)
  Sigma.nsz=cov(Ti.nsz)
  g1.nsz=-1/(2*mu1.nsz)
  g2.nsz=-1/(2*mu2.nsz)
  g3.nsz=1/(2*mu3.nsz)
  g4.nsz=1/(2*mu4.nsz)
  g5=1/mu5
  g6=-1/mu6
  g.nsz=c(g1.nsz,g2.nsz,g3.nsz,g4.nsz,g5,g6)
  V.Sigma.nsz=t(g.nsz)%*%Sigma.nsz%*%g.nsz
  sig.nsz=sqrt(V.Sigma.nsz)
  sig.nsz=as.vector(sig.nsz)
  ######################
  # compute bias corrections via generalized jackknife:
  # involving some earlier codes from Paul Wilson
  tbar=0
  tbar.nsz=0
  tbard21=rep(0,n)
  tbard22=rep(0,n)
  tbard11=rep(0,n)
  tbard12=rep(0,n)
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
    U1a=d21a*as.vector(y2a)
    U2a=d22a*as.vector(y2a)
    U3a=d11a*as.vector(y1a)
    U4a=d12a*as.vector(y1a)
    mu1a=mean(U1a)
    mu2a=mean(U2a)
    mu3a=mean(U3a)
    mu4a=mean(U4a)
    xia=0.5*(log(mu3a)+log(mu4a)-log(mu1a)-log(mu2a))+log(mu5)-log(mu6)
    ###
    d11a.nsz=d11a/epison11.b[1:na]
    d12a.nsz=d12a/epison12.b[1:na]
    d21a.nsz=d21a/epison21.b[1:na]
    d22a.nsz=d22a/epison22.b[1:na]
    U1a.nsz=d21a.nsz*as.vector(y2a)
    U2a.nsz=d22a.nsz*as.vector(y2a)
    U3a.nsz=d11a.nsz*as.vector(y1a)
    U4a.nsz=d12a.nsz*as.vector(y1a)
    mu1a.nsz=mean(U1a.nsz)
    mu2a.nsz=mean(U2a.nsz)
    mu3a.nsz=mean(U3a.nsz)
    mu4a.nsz=mean(U4a.nsz)
    xia.nsz=0.5*(log(mu3a.nsz)+log(mu4a.nsz)-log(mu1a.nsz)-log(mu2a.nsz))+log(mu5)-log(mu6)
    ################
    d11b= 1/rDEA::dea(X=x1b,Y=y1b,XREF=x1b,YREF=y1b, 
                      model="output", RTS="constant")$thetaOpt
    d12b= 1/rDEA::dea(X=x1b,Y=y1b,XREF=x2b,YREF=y2b, 
                      model="output", RTS="constant")$thetaOpt
    d21b= 1/rDEA::dea(X=x2b,Y=y2b,XREF=x1b,YREF=y1b, 
                      model="output", RTS="constant")$thetaOpt
    d22b= 1/rDEA::dea(X=x2b,Y=y2b,XREF=x2b,YREF=y2b, 
                      model="output", RTS="constant")$thetaOpt
    U1b=d21b*as.vector(y2b)
    U2b=d22b*as.vector(y2b)
    U3b=d11b*as.vector(y1b)
    U4b=d12b*as.vector(y1b)
    mu1b=mean(U1b)
    mu2b=mean(U2b)
    mu3b=mean(U3b)
    mu4b=mean(U4b)
    xib=0.5*(log(mu3b)+log(mu4b)-log(mu1b)-log(mu2b))+log(mu5)-log(mu6)
    #
    d11b.nsz=d11b/epison11.b[(na+1):n]
    d12b.nsz=d12b/epison12.b[(na+1):n]
    d21b.nsz=d21b/epison21.b[(na+1):n]
    d22b.nsz=d22b/epison22.b[(na+1):n]
    U1b.nsz=d21b.nsz*as.vector(y2b)
    U2b.nsz=d22b.nsz*as.vector(y2b)
    U3b.nsz=d11b.nsz*as.vector(y1b)
    U4b.nsz=d12b.nsz*as.vector(y1b)
    mu1b.nsz=mean(U1b.nsz)
    mu2b.nsz=mean(U2b.nsz)
    mu3b.nsz=mean(U3b.nsz)
    mu4b.nsz=mean(U4b.nsz)
    xib.nsz=0.5*(log(mu3b.nsz)+log(mu4b.nsz)-log(mu1b.nsz)-log(mu2b.nsz))+log(mu5)-log(mu6)
    #############
    tbar=tbar+(xia+xib)/2-xi
    tbar.nsz=tbar.nsz+(xia.nsz+xib.nsz)/2-xi.nsz
    #
    tbard21[ind1[1:na]]=tbard21[ind1[1:na]] +
      d21a - d21[ind1[1:na]]
    tbard21[ind1[(na+1):n]]=tbard21[ind1[(na+1):n]] +
      d21b - d21[ind1[(na+1):n]]
    tbard22[ind1[1:na]]=tbard22[ind1[1:na]] +
      d22a - d22[ind1[1:na]]
    tbard22[ind1[(na+1):n]]=tbard22[ind1[(na+1):n]] +
      d22b - d22[ind1[(na+1):n]]
    tbard11[ind1[1:na]]=tbard11[ind1[1:na]] +
      d11a - d11[ind1[1:na]]
    tbard11[ind1[(na+1):n]]=tbard11[ind1[(na+1):n]] +
      d11b - d11[ind1[(na+1):n]]
    tbard12[ind1[1:na]]=tbard12[ind1[1:na]] +
      d12a - d12[ind1[1:na]]
    tbard12[ind1[(na+1):n]]=tbard12[ind1[(na+1):n]] +
      d12b - d12[ind1[(na+1):n]]
    #
  }
  # compute the bias
  xi.bias=(1/L)*bc.fac*tbar
  tbard21=(1/L)*bc.fac*tbard21
  tbard22=(1/L)*bc.fac*tbard22
  tbard11=(1/L)*bc.fac*tbard11
  tbard12=(1/L)*bc.fac*tbard12
  #
  xi.bias.nsz=(1/L)*bc.fac*tbar.nsz
  tbard11.nsz=tbard11/epison11
  tbard12.nsz=tbard12/epison12
  tbard21.nsz=tbard21/epison21
  tbard22.nsz=tbard22/epison22
  ##############################
  # using szz (2023,JPA) method
  U1.szz=(d21-tbard21)*as.vector(y2)
  U2.szz=(d22-tbard22)*as.vector(y2)
  U3.szz=(d11-tbard11)*as.vector(y1)
  U4.szz=(d12-tbard12)*as.vector(y1)
  mu1.szz=mean(U1.szz)
  mu2.szz=mean(U2.szz)
  mu3.szz=mean(U3.szz)
  mu4.szz=mean(U4.szz)
  Ti.szz=cbind(U1.szz,U2.szz,U3.szz,U4.szz,U5,U6)
  Sigma.szz=cov(Ti.szz)
  g1.szz=-1/(2*mu1.szz)
  g2.szz=-1/(2*mu2.szz)
  g3.szz=1/(2*mu3.szz)
  g4.szz=1/(2*mu4.szz)
  g5=1/mu5
  g6=-1/mu6
  g.szz=c(g1.szz,g2.szz,g3.szz,g4.szz,g5,g6)
  V.Sigma.szz=t(g.szz)%*%Sigma.szz%*%g.szz
  sig.szz=sqrt(V.Sigma.szz)
  sig.szz=as.vector(sig.szz)
  #
  U1.nsz.szz=(d21.nsz-tbard21.nsz)*as.vector(y2)
  U2.nsz.szz=(d22.nsz-tbard22.nsz)*as.vector(y2)
  U3.nsz.szz=(d11.nsz-tbard11.nsz)*as.vector(y1)
  U4.nsz.szz=(d12.nsz-tbard12.nsz)*as.vector(y1)
  mu1.nsz.szz=mean(U1.nsz.szz)
  mu2.nsz.szz=mean(U2.nsz.szz)
  mu3.nsz.szz=mean(U3.nsz.szz)
  mu4.nsz.szz=mean(U4.nsz.szz)
  Ti.nsz.szz=cbind(U1.nsz.szz,U2.nsz.szz,U3.nsz.szz,U4.nsz.szz,U5,U6)
  Sigma.nsz.szz=cov(Ti.nsz.szz)
  g1.nsz.szz=-1/(2*mu1.nsz.szz)
  g2.nsz.szz=-1/(2*mu2.nsz.szz)
  g3.nsz.szz=1/(2*mu3.nsz.szz)
  g4.nsz.szz=1/(2*mu4.nsz.szz)
  g5=1/mu5
  g6=-1/mu6
  g.nsz.szz=c(g1.nsz.szz,g2.nsz.szz,g3.nsz.szz,g4.nsz.szz,g5,g6)
  V.Sigma.nsz.szz=t(g.nsz.szz)%*%Sigma.nsz.szz%*%g.nsz.szz
  sig.nsz.szz=sqrt(V.Sigma.nsz.szz)
  sig.nsz.szz=as.vector(sig.nsz.szz)
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  ts=sig/sqrt(n)
  ts.szz=sig.szz/sqrt(n)
  ts.nsz=sig.nsz/sqrt(n)
  ts.nsz.szz=sig.nsz.szz/sqrt(n)
  bounds.st=matrix((xi-ts*crit),nrow=3,ncol=2)
  if (np+nq<4) {
    bounds.psz=matrix((xi-xi.bias-ts*crit),nrow=3,ncol=2)
    bounds.szz=matrix((xi-xi.bias-ts.szz*crit),nrow=3,ncol=2)
    bounds.nsz=matrix((xi.nsz-xi.bias.nsz-ts.nsz*crit),nrow=3,ncol=2)
    bounds.nsz.szz=matrix((xi.nsz-xi.bias.nsz-ts.nsz.szz*crit),nrow=3,ncol=2)
    # make a list of results to return to calling routine and then quit:
    res=list(bounds.st=bounds.st,
             bounds.psz=bounds.psz,
             bounds.szz=bounds.szz,
             bounds.nsz=bounds.nsz,
             bounds.nsz.szz=bounds.nsz.szz,
             sig=c(sig,sig.szz,sig.nsz,sig.nsz.szz),
             estimate=c(xi,xi.nsz,xi.bias,xi.bias.nsz,xi-xi.bias,xi.nsz-xi.bias.nsz)
    )
  } else {
    U1.nk=d21[1:nk]*as.vector(y2)[1:nk]
    U2.nk=d22[1:nk]*as.vector(y2)[1:nk]
    U3.nk=d11[1:nk]*as.vector(y1)[1:nk]
    U4.nk=d12[1:nk]*as.vector(y1)[1:nk]
    U5.nk=as.vector(y2)[1:nk]
    U6.nk=as.vector(y1)[1:nk]
    mu1.nk=mean(U1.nk)
    mu2.nk=mean(U2.nk)
    mu3.nk=mean(U3.nk)
    mu4.nk=mean(U4.nk)
    mu5.nk=mean(U5.nk)
    mu6.nk=mean(U6.nk)
    xi.nk=0.5*(log(mu3.nk)+log(mu4.nk)-log(mu1.nk)-log(mu2.nk))+log(mu5.nk)-log(mu6.nk)
    ###
    U1.nsz.nk=d21.nsz[1:nk]*as.vector(y2)[1:nk]
    U2.nsz.nk=d22.nsz[1:nk]*as.vector(y2)[1:nk]
    U3.nsz.nk=d11.nsz[1:nk]*as.vector(y1)[1:nk]
    U4.nsz.nk=d12.nsz[1:nk]*as.vector(y1)[1:nk]
    mu1.nsz.nk=mean(U1.nsz.nk)
    mu2.nsz.nk=mean(U2.nsz.nk)
    mu3.nsz.nk=mean(U3.nsz.nk)
    mu4.nsz.nk=mean(U4.nsz.nk)
    xi.nsz.nk=0.5*(log(mu3.nsz.nk)+log(mu4.nsz.nk)-log(mu1.nsz.nk)-log(mu2.nsz.nk))+log(mu5.nk)-log(mu6.nk)
    ts.nk=sig/sqrt(nk)
    ts.szz.nk=sig.szz/sqrt(nk)
    ts.nsz.nk=sig.nsz/sqrt(nk)
    ts.nsz.szz.nk=sig.nsz.szz/sqrt(nk)
    bounds.psz=matrix((xi.nk-xi.bias-ts.nk*crit),nrow=3,ncol=2)
    bounds.szz=matrix((xi.nk-xi.bias-ts.szz.nk*crit),nrow=3,ncol=2)
    bounds.nsz=matrix((xi.nsz.nk-xi.bias.nsz-ts.nsz.nk*crit),nrow=3,ncol=2)
    bounds.nsz.szz=matrix((xi.nsz.nk-xi.bias.nsz-ts.nsz.szz.nk*crit),nrow=3,ncol=2)
    res=list(bounds.st=bounds.st,
             bounds.psz=bounds.psz,
             bounds.szz=bounds.szz,
             bounds.nsz=bounds.nsz,
             bounds.nsz.szz=bounds.nsz.szz,
             sig=c(sig,sig.szz,sig.nsz,sig.nsz.szz),
             estimate=c(xi.nk,xi.nsz.nk,xi.bias,xi.bias.nsz,xi.nk-xi.bias,xi.nsz.nk-xi.bias.nsz)
    )
  }
  return(res)
}