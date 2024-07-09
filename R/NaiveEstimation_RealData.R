#' Estimating parameters along with corresponding variances based on Naive model with real data.
#'
#' @param ITER The number of simulation runs
#' @param MHiteration The number of iterations in Metropolis–Hastings algorithm
#' @param eps Stopping value for MCECM algorithm
#' @param mm Number of areas.
#' @param time Maximum time.
#' @param MuxStar0 Mean vector of unobserved areal level covariates.
#' @param MuxInd0 Mean vector of unobserved individual level covariates.
#' @param SigmaxStar0 Variance of unobserved areal level covariates.
#' @param SigmaxInd0 Variance of unobserved individual level covariates.
#' @param Sigmav0 Variance of areal level measurement error variable.
#' @param Sigmaw0 Variance of individual level measurement error variable.
#' @param lambda0 Spatial dependency parameter.
#' @param sigma0 Over dispersion parameter.
#' @param delta0 The spatial parameter.
#' @param alpha0 Initial value for intercept.
#' @param beta10 Initial value for coefficient of observed individual level covariates.
#' @param beta20 Initial value for coefficient of observed areal level covariates.
#' @param beta30 Initial value for coefficient of unobserved individual level covariates.
#' @param beta40 Initial value for coefficient of unobserved areal level covariates.
#' @param InfPeriod The infectious period length.
#' @param Di Euclidean distance between individuals
#' @param D Neibourhood structure
#' @param Nlabel Label for each sample from the area
#' @param n Total number of individuals
#' @param cov1 observed individual level covariates
#' @param cov2 observed areal level covariates
#' @param ww Unobserved individual level covariates
#' @param vv unobserved areal level covariates
#' @param tau tau
#'
#' @return The results of the function
#'
#' @export
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom MASS mvrnorm
#' @importFrom numDeriv hessian
#' @importFrom stats optim
#' @importFrom stats rnorm
#' @importFrom stats dnorm
#' @importFrom stats runif
#' @importFrom ngspatial adjacency.matrix
#' @importFrom corpcor make.positive.definite
#' @importFrom psych tr
#'
#'
#' @examples NaiveEstimation_RealData(1,5,0.05,4,20,0.1,0.1,0.15,0.8,0.6,0.6,
#' 0.85,1.1,2.7,0,1,0,1,1,3,
#' matrix(runif(900,min = 4,max = 20),nrow=30, byrow = TRUE),
#' matrix(c(2,-1,-1,0,-1,2,0,-1,-1,0,2,-1,0,-1,-1,2),nrow=4,byrow=TRUE),
#' rep(1:4,c(7,6,8,9)),30,runif(30, 0, 1),
#' runif(4,0,1),runif(30,-2,2),runif(4,0,1),
#' sample(c(0,1),replace = TRUE, size = 30))
#'
#'
NaiveEstimation_RealData=function(ITER,MHiteration,eps,mm,time,MuxStar0,MuxInd0,SigmaxStar0,SigmaxInd0,Sigmav0,Sigmaw0,lambda0,sigma0,delta0,alpha0,beta10,beta20,beta30,beta40,InfPeriod,Di,D,Nlabel,n,cov1,cov2,ww,vv,tau){


  ############################a######Generate Data############################
  nn=MHiteration
  NaivSimPAR=matrix(0,ITER,8)
  NaivSimVarPar=matrix(0,ITER,6)
  NaivCSimVarPar=matrix(0,ITER,6)
  NaivSimInfPar=matrix(0,ITER,6)
  NaivSimVarSpatial=matrix(0,ITER,2)
  NaivSimPos=list()


  lambda=rep(InfPeriod,n) #########infected period##############
  I=diag(mm)
  mu=rep(0,mm)
  Sigma0=sigma0^2*solve((lambda0*D+(1-lambda0)*I))
  Sigma0 = make.positive.definite(Sigma0, tol = 1e-3)
  phi = mvrnorm(1,mu,Sigma0,tol=1e-6)
  ###############################Covariates without measurement error#######################

  Newlabel=matrix(0,n,mm)
  for(mmm in 1:mm){
    for (i in 1:n) {
      if(D[Nlabel[i],mmm]!=0){
        Newlabel[i,mmm] = -1
      }
    }
  }

  CovXStar1=rnorm(mm, MuxStar0, sqrt(SigmaxStar0))
  CovX=rnorm(n, MuxInd0, sqrt(SigmaxInd0))

  ############################## f_y|x,xStar,u###############################
  Fy1=function(IndLe,AreLe,RaE,alpha1,beta1,beta2,beta3,beta4,delta,time,tau,mm,lambda){
    fy1=array(0,c(n,time,mm))
    for(i in 1:n){
      for(t in 1:time){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){


              dx1=rep(0,n)
              for(j in 1:n){
                if (Newlabel[j,mmm]!=0){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx1[j]=Di[i,j]^(-delta)
                  }
                }
              }
              dx1=replace(dx1,is.infinite(dx1),0)
              dx=sum(dx1)

              prob1=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*IndLe[i]+beta4*AreLe[mmm]+RaE[mmm])*dx)

              fy1[i,t,mmm]=(1-prob1)
            }

            if (tau[i]==(t+1)){

              dx1=rep(0,n)
              for(j in 1:n){
                if (Newlabel[j,mmm]!=0){
                  if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
                    dx1[j]=Di[i,j]^(-delta)
                  }
                }
              }
              dx1=replace(dx1,is.infinite(dx1),0)
              dx=sum(dx1)

              prob1=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*IndLe[i]+beta4*AreLe[mmm]+RaE[mmm])*dx)

              fy1[i,t,mmm]=prob1
            }
          }
        }
      }
    }

    PP=c()
    H=list()
    fy2=c()
    for(mmm in 1:mm){
      for(i in 1:n){
        PP[i]=round(prod(fy1[i,,mmm][fy1[i,,mmm]>0]),10)
      }
      H[[mmm]]=which(PP!=1)
      fy2[mmm]=prod(PP[H[[mmm]]][PP[H[[mmm]]]>0])
    }

    return(fy2)

  }




  ##############################################################################################
  fRanE=function(RaE,mmm,lambda1,sigma1,D){
    Nh=-D[,mmm]*RaE
    dnorm(RaE[mmm],lambda1/(1-lambda1+D[mmm,mmm]*lambda1)*sum(Nh[-mmm]),sqrt(1/(sigma1^2*(1-lambda1+D[mmm,mmm]*lambda1))), log = FALSE)
  }

  #################Sum h(z-z_{I}) #######################
  SumH=function(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda){

    SumH1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH1[j]=Di[i,j]^(-delta)
        }
      }
    }
    SumH1=replace(SumH1,is.infinite(SumH1),0)
    SumH2=sum(SumH1)

    SumH3=exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm])*SumH2

    return(SumH3)
  }


  #################Sum h(z-z_{I})log(d_{ij}) #######################
  SumHnew1=function(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda){

    SumH4=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH4[j]=Di[i,j]^(-delta)*(log(Di[i,j]))
        }
      }
    }
    SumH4=replace(SumH4,is.infinite(SumH4),0)
    SumH5=sum(SumH4)

    SumH6=exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm])*SumH5

    return(SumH6)
  }



  #################Sum h(z-z_{I})(log(d_{ij}))^2 #######################
  SumHnew2=function(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda){


    SumH7=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          SumH7[j]=Di[i,j]^(-delta)*(log(Di[i,j]))^2
        }
      }
    }
    SumH7=replace(SumH7,is.infinite(SumH7),0)
    SumH8=sum(SumH7)

    SumH9=exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm])*SumH8

    return(SumH9)
  }





  E.sigma.first=function(Pos,D,lambda1){

    d1=c()
    for(L in 1:nn){
      d1[L]=as.numeric(t(Pos[L,(mm+1):(2*mm)])%*%(lambda1*D+(1-lambda1)*I)%*%Pos[L,(mm+1):(2*mm)])
    }
    SMeansigma1=mean(d1)
    return(SMeansigma1)
  }


  E.lambda1.first=function(Pos,D){

    d2=c()
    for(L in 1:nn){
      d2[L]= as.numeric(t(Pos[L,(mm+1):(2*mm)])%*%(D-I)%*%Pos[L,(mm+1):(2*mm)])
    }
    SMeanlambda11=mean(d2)
    return(SMeanlambda11)

  }


  ######################################### Expectation #######################
  mean1=function(Pos,beta3,beta4,mmm,i){
    d1=c()
    for(L in 1:nn){
      d1[L]=exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }
    SMean1=mean(d1)
    return(SMean1)
  }




  meanXstar1=function(Pos,beta3,beta4,mmm,i){
    d1=c()
    for(L in 1:nn){
      d1[L]=Pos[L,mmm]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }
    SMeanXstar1=mean(d1)
    return(SMeanXstar1)
  }



  meanX1=function(Pos,beta3,beta4,mmm,i){
    d1=c()
    for(L in 1:nn){
      d1[L]=Pos[L,(2*mm+i)]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }
    SMeanX1=mean(d1)
    return(SMeanX1)
  }




  meanXstar3=function(Pos,beta3,beta4,mmm,i){
    d1=c()
    for(L in 1:nn){
      d1[L]=Pos[L,mmm]^2*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }
    SMeanXstar3=mean(d1)
    return(SMeanXstar3)
  }




  meanX3=function(Pos,beta3,beta4,mmm,i){
    d1=c()
    for(L in 1:nn){
      d1[L]=Pos[L,(2*mm+i)]^2*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }
    SMeanX3=mean(d1)
    return(SMeanX3)
  }

  ###############################################################

  mean2=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,(mmm+mm)])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMean2=0}

    if(prob11!=0){
      SMean2=(1-prob11)/prob11*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }

    return(SMean2)
  }




  meanXstar2=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMeanXstar2=0}

    if(prob11!=0){
      SMeanXstar2=Pos[L,mmm]*(1-prob11)/prob11*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }

    return(SMeanXstar2)
  }





  meanX2=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMeanX2=0}

    if(prob11!=0){
      SMeanX2=Pos[L,(2*mm+i)]*(1-prob11)/prob11*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }

    return(SMeanX2)
  }



  meanXstar4=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMeanXstar4=0}

    if(prob11!=0){
      SMeanXstar4=Pos[L,mmm]^2*(1-prob11)/prob11*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }

    return(SMeanXstar4)
  }





  meanX4=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMeanX4=0}

    if(prob11!=0){
      SMeanX4=Pos[L,(2*mm+i)]^2*(1-prob11)/prob11*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
    }

    return(SMeanX4)
  }

  ####################################################################



  mean3=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMean3=0}

    if(prob11!=0){

      SMean3=(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]++2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
    }

    return(SMean3)
  }



  meanXstar5=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMean5=0}

    if(prob11!=0){

      SMean5=Pos[L,mmm]^2*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
    }

    return(SMean5)
  }





  meanX5=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


    dx1=rep(0,n)
    for(j in 1:n){
      if (Newlabel[j,mmm]!=0){
        if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
          dx1[j]=Di[i,j]^(-delta)
        }
      }
    }
    dx1=replace(dx1,is.infinite(dx1),0)
    dx=sum(dx1)

    prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
    if (is.na(prob11)) {
      prob11 <- 0
    }

    if(prob11==0){SMean5=0}

    if(prob11!=0){

      SMean5=Pos[L,(2*mm+i)]^2*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
    }

    return(SMean5)
  }


  ############################################################################################################################

  estfunNaiv=function(Nlabel,phi,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,D,lambda){



    ##################M-H###############################
    nn=MHiteration
    AreLe0=vv
    IndLe0=ww

    mu=rep(0,mm)
    Sigma11=sigma1^2*solve((lambda1*D+(1-lambda1)*I))
    Sigma1 = make.positive.definite(Sigma11, tol = 1e-3)
    RaE0=mvrnorm(1, mu, Sigma11, tol = 1e-6)


    First=c(AreLe0,RaE0,IndLe0)
    Pos=matrix(0,nrow=nn+1,ncol=((2*mm)+n))
    Pos[1,]=First


    MPHprob3=matrix(0,nn,n)
    Uni3=matrix(0,nn,n)

    MPHprob1=matrix(0,nn,mm)
    Uni1=matrix(0,nn,mm)

    MPHprob2=c()
    Uni2=c()
    for (L in 2:nn){

      AreLe=vv
      IndLe=ww


      RaE=mvrnorm(1, mu, Sigma11, tol = 1e-6)

      Pos[L,(2*mm+1):(2*mm+n)]=IndLe

      Pos[L,1:mm]=AreLe


      Area_sec3=Fy1(Pos[L-1,(2*mm+1):(2*mm+n)],Pos[L-1,1:mm],RaE,alpha1,beta1,beta2,beta3,beta4,delta,time,tau,mm,lambda)
      Area_sec4=Fy1(Pos[L-1,(2*mm+1):(2*mm+n)],Pos[L-1,1:mm],Pos[L-1,(mm+1):(2*mm)],alpha1,beta1,beta2,beta3,beta4,delta,time,tau,mm,lambda)
      Uni2[L]=runif(1,0,1)


      MPHprob2[L]=min(1,prod(Area_sec3/Area_sec4))
      if(Uni2[L]<MPHprob2[L]){
        Pos[L,(mm+1):(2*mm)]=RaE
      }

      if(Uni2[L]>=MPHprob2[L]){
        Pos[L,(mm+1):(2*mm)]=Pos[L-1,(mm+1):(2*mm)]
      }



    }


    ######################################## alpha #######################################
    NewtonAlpha=function(Newalpha1){

      A1=rep(0,time)
      for(t in 1:time){
        A2=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                A2[i]=-as.numeric(SumH(Nlabel,Di,Newalpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        A1[t]=sum(A2)
      }
      SusA1=sum(A1)


      ###################### infectious period ###################

      A3=rep(0,time)
      for(t in 1:time){
        A4=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                SA4=c()
                for(L in 1:nn){
                  SA4[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }
                A4[i]=mean(SA4)*SumH(Nlabel,Di,Newalpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)
              }
            }
          }
        }
        A3[t]=sum(A4)
      }
      InfA3=sum(A3,na.rm=T)
      EndA3=SusA1+InfA3
      ######################### Second Der...########################
      ##################### infectious period ###################

      A5=rep(0,time)
      for(t in 1:time){
        A6=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                SA6=c()
                for(L in 1:nn){
                  SA6[L]=-mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }
                A6[i]=mean(SA6)*SumH(Nlabel,Di,Newalpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2
              }
            }
          }
        }
        A5[t]=sum(A6)
      }
      InfA5=sum(A5,na.rm=T)
      EndA5=EndA3+InfA5

      EstAlpha1=Newalpha1-EndA3/EndA5
      diffAlpha=EstAlpha1-Newalpha1
      result=list(Newalpha1=EstAlpha1,diffAlpha=diffAlpha)
      result
    }

    talpha=1
    Nalpha0=NewtonAlpha(alpha1)
    diffAlpha=Nalpha0$diffAlpha
    Newalpha1=Nalpha0$Newalpha1

    EstAlpha=Newalpha1
    #################################beta1##################################
    NewtonBeat1=function(Newbeta1){
      Beta11=rep(0,time)
      for(t in 1:time){
        Beta12=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta12[i]=-cov1[i]*as.numeric(SumH(Nlabel,Di,EstAlpha,Newbeta1,beta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta11[t]=sum(Beta12)
      }
      SusBeta11=sum(Beta11)


      ###################### infectious period ###################

      Beta13=rep(0,time)
      for(t in 1:time){
        Beta14=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                SB4=c()
                for(L in 1:nn){
                  SB4[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }
                Beta14[i]=cov1[i]*mean(SB4)*SumH(Nlabel,Di,EstAlpha,Newbeta1,beta2,delta,i,mmm,t,tau,lambda)
              }
            }
          }
        }
        Beta13[t]=sum(Beta14)
      }
      InfBeta13=sum(Beta13,na.rm=T)
      EndBeta13=SusBeta11+InfBeta13
      ######################### Second Der...########################
      ##################### infectious period ###################


      Beta112=rep(0,time)
      for(t in 1:time){
        Beta122=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta122[i]=-cov1[i]^2*as.numeric(SumH(Nlabel,Di,EstAlpha,Newbeta1,beta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta112[t]=sum(Beta122)
      }
      SusBeta112=sum(Beta112)




      Beta15=rep(0,time)
      for(t in 1:time){
        Beta16=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                SB6=c()
                for(L in 1:nn){
                  SB6[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }


                SB8=c()
                for(L in 1:nn){
                  SB8[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }


                Beta16[i]=cov1[i]^2*as.numeric(SumH(Nlabel,Di,EstAlpha,Newbeta1,beta2,delta,i,mmm,t,tau,lambda)*mean(SB6))-cov1[i]^2*as.numeric(SumH(Nlabel,Di,EstAlpha,Newbeta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(SB8))
              }
            }
          }
        }
        Beta15[t]=sum(Beta16)
      }
      InfBeta15=sum(Beta15,na.rm=T)
      EndBeta15=SusBeta112+InfBeta15

      EstBeta1New1=Newbeta1-EndBeta13/EndBeta15


      diffBeat1=EstBeta1New1-Newbeta1
      result=list(Newbeta1=EstBeta1New1,diffBeat1=diffBeat1)
      result
    }



    if(is.na(EstAlpha)){
      EstBeta1=NA

    }

    if(!is.na(EstAlpha)){

      tbeta1=1
      Nbeta10=NewtonBeat1(beta1)
      diffBeat1=Nbeta10$diffBeat1
      Newbeta1=Nbeta10$Newbeta1
      EstBeta1=Newbeta1
    }

    #################################beta2##################################
    Newbeta2=beta2
    NewtonBeat2=function(Newbeta2){

      Beta21=rep(0,time)
      for(t in 1:time){
        Beta22=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta22[i]=-cov2[mmm]*as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,Newbeta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta21[t]=sum(Beta22)
      }
      SusBeta21=sum(Beta21)


      ###################### infectious period ###################

      Beta23=rep(0,time)
      for(t in 1:time){
        Beta24=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                SB10=c()
                for(L in 1:nn){
                  SB10[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                Beta24[i]=cov2[mmm]*as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,Newbeta2,delta,i,mmm,t,tau,lambda)*mean(SB10))
              }
            }
          }
        }
        Beta23[t]=sum(Beta24)
      }
      InfBeta23=sum(Beta23,na.rm=T)
      EndBeta23=SusBeta21+InfBeta23

      ######################### Second Der...########################
      ##################### infectious period ###################


      Beta212=rep(0,time)
      for(t in 1:time){
        Beta222=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta222[i]=-cov2[mmm]^2*as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,Newbeta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta212[t]=sum(Beta222)
      }
      SusBeta212=sum(Beta212)




      Beta25=rep(0,time)
      for(t in 1:time){
        Beta26=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){



                SB12=c()
                for(L in 1:nn){
                  SB12[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }


                SB14=c()
                for(L in 1:nn){
                  SB14[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }


                Beta26[i]=cov2[mmm]^2*as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,Newbeta2,delta,i,mmm,t,tau,lambda)*mean(SB12))-cov2[mmm]^2*as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,Newbeta2,delta,i,mmm,t,tau,lambda)^2*mean(SB14))
              }
            }
          }
        }
        Beta25[t]=sum(Beta26)
      }
      InfBeta25=sum(Beta25,na.rm=T)
      EndBeta25=SusBeta212+InfBeta25
      ##################################################################

      EstBeta2New1=Newbeta2-EndBeta23/EndBeta25

      diffBeat2=EstBeta2New1-Newbeta2
      result=list(Newbeta2=EstBeta2New1,diffBeat2=diffBeat2)
      result
    }

    if(is.na(EstBeta1)){
      EstBeta2=NA

    }

    if(!is.na(EstBeta1)){
      tbeta2=1
      Nbeta20=NewtonBeat2(beta2)
      diffBeat2=Nbeta20$diffBeat2
      Newbeta2=Nbeta20$Newbeta2
      EstBeta2=Newbeta2
    }



    if(is.na(EstBeta2)){
      EstBeta3=NA

    }





    if(!is.na(EstBeta2)){


      #################################beta4##################################
      Beta31=rep(0,time)
      for(t in 1:time){
        Beta32=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta32[i]=-as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*meanX1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta31[t]=sum(Beta32)
      }
      SusBeta31=sum(Beta31)


      ###################### infectious period ###################

      Beta33=rep(0,time)
      for(t in 1:time){
        Beta34=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                SXB13=c()
                for(L in 1:nn){
                  SXB13[L]=meanX2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                Beta34[i]=as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*mean(SXB13))
              }
            }
          }
        }
        Beta33[t]=sum(Beta34)
      }
      InfBeta33=sum(Beta33,na.rm=T)
      EndBeta33=SusBeta31+InfBeta33
      ######################### Second Der...########################
      ##################### infectious period ###################


      Beta312=rep(0,time)
      for(t in 1:time){
        Beta322=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){



                Beta322[i]=-as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*meanX3(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta312[t]=sum(Beta322)
      }
      SusBeta312=sum(Beta312)




      Beta35=rep(0,time)
      for(t in 1:time){
        Beta36=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                SXB23=c()
                for(L in 1:nn){
                  SXB23[L]=meanX4(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                SXB53=c()
                for(L in 1:nn){
                  SXB53[L]=meanX5(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                Beta36[i]=as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*mean(SXB23))-as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)^2*mean(SXB53))
              }
            }
          }
        }
        Beta35[t]=sum(Beta36)
      }
      InfBeta35=sum(Beta35,na.rm=T)
      EndBeta35=SusBeta312+InfBeta35

      EstBeta3=beta3-EndBeta33/EndBeta35
      EstBeta3
    }
    ##################################################################










    if(is.na(EstBeta3)){
      EstBeta4=NA

    }





    if(!is.na(EstBeta3)){


      #################################beta4##################################
      Beta41=rep(0,time)
      for(t in 1:time){
        Beta42=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta42[i]=-as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*meanXstar1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta41[t]=sum(Beta42)
      }
      SusBeta41=sum(Beta41)


      ###################### infectious period ###################

      Beta43=rep(0,time)
      for(t in 1:time){
        Beta44=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                SXB2=c()
                for(L in 1:nn){
                  SXB2[L]=meanXstar2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                Beta44[i]=as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*mean(SXB2))
              }
            }
          }
        }
        Beta43[t]=sum(Beta44)
      }
      InfBeta43=sum(Beta43,na.rm=T)
      EndBeta43=SusBeta41+InfBeta43
      ######################### Second Der...########################
      ##################### infectious period ###################


      Beta412=rep(0,time)
      for(t in 1:time){
        Beta422=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){



                Beta422[i]=-as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*meanXstar3(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        Beta412[t]=sum(Beta422)
      }
      SusBeta412=sum(Beta412)




      Beta45=rep(0,time)
      for(t in 1:time){
        Beta46=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                SXB4=c()
                for(L in 1:nn){
                  SXB4[L]=meanXstar4(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                SXB5=c()
                for(L in 1:nn){
                  SXB5[L]=meanXstar5(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }

                Beta46[i]=as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)*mean(SXB4))-as.numeric(SumH(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,delta,i,mmm,t,tau,lambda)^2*mean(SXB5))
              }
            }
          }
        }
        Beta45[t]=sum(Beta46)
      }
      InfBeta45=sum(Beta45,na.rm=T)
      EndBeta45=SusBeta412+InfBeta45

      EstBeta4=beta4-EndBeta43/EndBeta45
      EstBeta4
    }
    ##################################################################



    ################################# Delta #############################

    ################## susceptible period ####################
    NewtonDelta=function(NewDelta){


      S8=rep(0,time)
      for(t in 1:time){
        S7=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                S7[i]=as.numeric(SumHnew1(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,NewDelta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        S8[t]=sum(S7)
      }
      SusNew8=sum(S8)



      ###################### infectious period ###################

      I10=rep(0,time)
      for(t in 1:time){
        I9=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                SXD2=c()
                for(L in 1:nn){
                  SXD2[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }



                I9[i]=-as.numeric(SumHnew1(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,NewDelta,i,mmm,t,tau,lambda)*mean(SXD2))
              }
            }
          }
        }
        I10[t]=sum(I9)
      }
      InfeNew10=sum(I10,na.rm=T)

      EndNew10=SusNew8+InfeNew10
      ######################### Second Der...########################

      S10=rep(0,time)
      for(t in 1:time){
        S9=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                S9[i]=-as.numeric(SumHnew2(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,NewDelta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
              }
            }
          }
        }
        S10[t]=sum(S9)
      }
      SusNew10=sum(S10)

      ###################### infectious period ###################

      I12=rep(0,time)
      for(t in 1:time){
        I11=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                SXD4=c()
                for(L in 1:nn){
                  SXD4[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }


                SXD6=c()
                for(L in 1:nn){
                  SXD6[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
                }


                I11[i]=as.numeric(SumHnew2(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,NewDelta,i,mmm,t,tau,lambda)*mean(SXD4))-as.numeric(SumHnew1(Nlabel,Di,EstAlpha,EstBeta1,EstBeta2,NewDelta,i,mmm,t,tau,lambda)^2*mean(SXD6))
              }
            }
          }
        }
        I12[t]=sum(I11)
      }
      InfeNew12=sum(I12,na.rm=T)

      EndNew12=SusNew10+InfeNew12

      Estdelta1=NewDelta-EndNew10/EndNew12

      diffDelta=Estdelta1-NewDelta
      result=list(NewDelta=Estdelta1,diffDelta=diffDelta)
      result
    }

    if(is.na(EstBeta4)){
      Estdelta=NA

    }

    if(!is.na(EstBeta4)){
      tdelta=1
      Ndelta0=NewtonDelta(delta)
      diffDelta=Ndelta0$diffDelta
      NewDelta=Ndelta0$NewDelta

      Estdelta=NewDelta
    }

    if(is.na(Estdelta)){
      HatSigmmaU=NA
      EstGammau=NA
      MuxStarHat=NA
      SigmaxStarHat=NA
      MuxHat=NA
      SigmaIndHat=NA
    }

    if(!is.na(Estdelta)){
      ##################################### Sigma_U #####################
      ##################################### Sigma_U #####################
      ##################################### Sigma_U #####################
      ##################################### Sigma_U ####################

      Loglik11 <- function(par) {
        lambda1 <- par[1]
        sigma1 <- par[2]

        if (lambda1 <= 0 || lambda1 >= 1 || sigma1 <= 0) {
          return(Inf)
        }

        Log1 <- numeric(nn)
        for (L in 1:nn) {
          mu <- rep(0, mm)
          sigma_matrix <- sigma1^2 * solve(lambda1 * D + (1 - lambda1) * I)

          if (any(is.na(sigma_matrix)) || !all(eigen(sigma_matrix)$values > 0)) {
            return(Inf)
          }

          Log1[L] <- dmvnorm(Pos[L, (mm + 1):(2 * mm)], mean = mu, sigma = sigma_matrix, log = TRUE)
        }
        LLL <- -mean(Log1)
        return(LLL)
      }


      init =c(lambda1,sigma1)
      ml <- optim( init,Loglik11)
      EstU1=ml$par

      HatSigmmaU=EstU1[2]
      EstGammau=EstU1[1]


      #############################################################
      #############################################################

      MuxHat=0
      SigmaIndHat=0

      MuxStarHat=0
      SigmaxStarHat=0

    }
    #############################################################
    result=list(Pos1=Pos,beta1=EstBeta1,beta2=EstBeta2,beta3=EstBeta3,beta4=EstBeta4,alpha1=EstAlpha,delta=Estdelta,sigma1=HatSigmmaU,lambda1=EstGammau)
    result
  }





  LogNaiv=function(Nlabel,Pos1,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,D,lambda){
    nn=MHiteration
    SumH=function(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda){

      SumH1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            SumH1[j]=Di[i,j]^(-delta)
          }
        }
      }
      SumH1=replace(SumH1,is.infinite(SumH1),0)
      SumH2=sum(SumH1)

      SumH3=exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm])*SumH2

      return(SumH3)
    }


    mean1L=function(Pos1,beta3,beta4,mmm,i){
      d1=c()
      for(L in 1:nn){
        d1[L]=exp(beta3*Pos1[L,(2*mm+i)]+beta4*Pos1[L,mmm]+Pos1[L,mmm+mm])
      }
      SMean1=mean(d1)
      return(SMean1)
    }


    mean2L=function(Nlabel,Pos1,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos1[L,(2*mm+i)]+beta4*Pos1[L,mmm]+Pos1[L,(mmm+mm)])*dx)
      if(prob11==0){SMean2=0}

      if(prob11!=0){
        SMean2=log(prob11)
      }

      return(SMean2)
    }



    A1=rep(0,time)
    for(t in 1:time){
      A2=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){
              A2[i]=-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean1L(Pos1,beta3,beta4,mmm,i))
            }
          }
        }
      }
      A1[t]=sum(A2)
    }
    SusA1=sum(A1)


    ###################### infectious period ###################

    A3=rep(0,time)
    for(t in 1:time){
      A4=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if(Nlabel[i]==mmm){
            if(tau[i]==(t+1)){
              SA4=c()
              for(L in 1:nn){
                SA4[L]=mean2L(Nlabel,Pos1,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
              A4[i]=mean(SA4)
            }
          }
        }
      }
      A3[t]=sum(A4)
    }
    InfA3=sum(A3,na.rm=T)
    EndA3L=SusA1+InfA3




    Log1=c()
    for(L in 1:nn){
      Log1[L]=dmvnorm(Pos1[L,(mm+1):(2*mm)], rep(0,mm), sigma1^2*solve((lambda1*D+(1-lambda1)*I)), log = T)
    }
    LLL=mean(Log1)





    loglik1=LLL+EndA3L
    return(loglik1)

  }
  ##################################################################################

  ObsInfNaiv=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,lambda){
    nn=MHiteration

    ##################################################### Inf alpha ###########################################################
    ##################################################### Cov ###########################################################
    alphaS=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        alphaS1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                alphaS1[t]=-SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        alphaS[i,L]=sum(alphaS1)
      }
    }

    alphaI=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        alphaI1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                alphaI1[t]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)
              }
            }
          }
        }
        alphaI[i,L]=sum(alphaI1)
      }
    }

    alphaIS=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        alphaIS[i,L]=alphaS[i,L]+alphaI[i,L]
      }
    }

    SecondAlpha=c()
    for(i in 1:n){
      SecondAlpha[i]=mean(alphaIS[i,]^2)
    }

    FirstAlpha=c()
    for(i in 1:n){
      FirstAlpha[i]=mean(alphaIS[i,])^2
    }

    MAlpha=c()
    for(i in 1:n){
      MAlpha[i]=mean(alphaIS[i,])
    }
    ##################################################################################


    IA5=rep(0,time)
    for(t in 1:time){
      IA6=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){
              ISA6=c()
              for(L in 1:nn){
                ISA6[L]=-mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
              IA6[i]=mean(ISA6)*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2
            }
          }
        }
      }
      IA5[t]=sum(IA6)
    }
    IInfA5=sum(IA5,na.rm=T)
    IEndA5=IInfA5+sum(MAlpha)

    ######################################
    InfAlpha=-IEndA5-sum(SecondAlpha)+sum(FirstAlpha)
    ##################################################################################################################
    ####################################### Info beta1 ##############################
    ####################################### New beta1##############################
    ###############################Cov#######################################
    Beta1S=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta1S1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta1S1[t]=-cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        Beta1S[i,L]=sum(Beta1S1)
      }
    }

    ###################### infectious period ###################
    Beta1I=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta1I1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                Beta1I1[t]=cov1[i]*mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)
              }
            }
          }
        }
        Beta1I[i,L]=sum(Beta1I1)
      }
    }



    Beta1IS=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta1IS[i,L]=Beta1S[i,L]+Beta1I[i,L]
      }
    }

    SecondBeta1=c()
    for(i in 1:n){
      SecondBeta1[i]=mean(Beta1IS[i,]^2)
    }

    FirstBeta1=c()
    for(i in 1:n){
      FirstBeta1[i]=mean(Beta1IS[i,])^2
    }

    MBeta1=c()
    for(i in 1:n){
      MBeta1[i]=mean(Beta1IS[i,])
    }
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    ######################### Second Der...########################
    ##################### infectious period ###################


    IBeta112=rep(0,time)
    for(t in 1:time){
      IBeta122=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){
              IBeta122[i]=-cov1[i]^2*as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
            }
          }
        }
      }
      IBeta112[t]=sum(IBeta122)
    }
    ISusBeta112=sum(IBeta112)




    IBeta15=rep(0,time)
    for(t in 1:time){
      IBeta16=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              ISB6=c()
              for(L in 1:nn){
                ISB6[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              ISB8=c()
              for(L in 1:nn){
                ISB8[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              IBeta16[i]=cov1[i]^2*as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean(ISB6))-cov1[i]^2*as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(ISB8))
            }
          }
        }
      }
      IBeta15[t]=sum(IBeta16)
    }
    IInfBeta15=sum(IBeta15,na.rm=T)
    IEndBeta15=ISusBeta112+IInfBeta15

    ######################################
    InfBeta1=-IEndBeta15-sum(SecondBeta1)+sum(FirstBeta1)
    ###########################################################################################
    ####################################### Info beta2 ##############################
    #######################################Cov##############################
    Beta2S=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta2S1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta2S1[t]=-cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        Beta2S[i,L]=sum(Beta2S1)
      }
    }

    ###################### infectious period ###################
    Beta2I=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta2I1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                Beta2I1[t]=cov2[mmm]*mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)
              }
            }
          }
        }
        Beta2I[i,L]=sum(Beta2I1)
      }
    }



    Beta2IS=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta2IS[i,L]=Beta2S[i,L]+Beta2I[i,L]
      }
    }

    SecondBeta2=c()
    for(i in 1:n){
      SecondBeta2[i]=mean(Beta2IS[i,]^2)
    }

    FirstBeta2=c()
    for(i in 1:n){
      FirstBeta2[i]=mean(Beta2IS[i,])^2
    }

    MBeta2=c()
    for(i in 1:n){
      MBeta2[i]=mean(Beta2IS[i,])
    }
    ########################################################################################
    ########################################################################################
    ######################### Second Der...########################
    ##################### infectious period ###################

    IBeta212=rep(0,time)
    for(t in 1:time){
      IBeta222=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){
              IBeta222[i]=-cov2[mmm]^2*as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
            }
          }
        }
      }
      IBeta212[t]=sum(IBeta222)
    }
    ISusBeta212=sum(IBeta212)




    IBeta25=rep(0,time)
    for(t in 1:time){
      IBeta26=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){


              ISB12=c()
              for(L in 1:nn){
                ISB12[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              ISB14=c()
              for(L in 1:nn){
                ISB14[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              IBeta26[i]=cov2[mmm]^2*as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean(ISB12))-cov2[mmm]^2*as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(ISB14))
            }
          }
        }
      }
      IBeta25[t]=sum(IBeta26)
    }
    IInfBeta25=sum(IBeta25,na.rm=T)
    IEndBeta25=ISusBeta212+IInfBeta25

    ######################################
    InfBeta2=-IEndBeta25-sum(SecondBeta2)+sum(FirstBeta2)
    ###########################################################################################
    ####################################### Info beta3 ##############################
    ####################################### Cov##############################

    Beta3S=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta3S1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta3S1[t]=-SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,(2*mm+i)]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        Beta3S[i,L]=sum(Beta3S1)
      }
    }


    ###################### infectious period ###################

    Beta3I=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta3I1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                Beta3I1[t]=SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanX2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        Beta3I[i,L]=sum(Beta3I1)
      }
    }


    Beta3IS=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta3IS[i,L]=Beta3S[i,L]+Beta3I[i,L]
      }
    }

    SecondBeta3=c()
    for(i in 1:n){
      SecondBeta3[i]=mean(Beta3IS[i,]^2)
    }

    FirstBeta3=c()
    for(i in 1:n){
      FirstBeta3[i]=mean(Beta3IS[i,])^2
    }

    MBeta3=c()
    for(i in 1:n){
      MBeta3[i]=mean(Beta3IS[i,])
    }
    ###############################################################################################
    ###############################################################################################
    #################################beta4##################################
    ##################### infectious period ###################


    IBeta312=rep(0,time)
    for(t in 1:time){
      IBeta322=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){
              IBeta322[i]=-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanX3(Pos,beta3,beta4,mmm,i))
            }
          }
        }
      }
      IBeta312[t]=sum(IBeta322)
    }
    ISusBeta312=sum(IBeta312)




    IBeta35=rep(0,time)
    for(t in 1:time){
      IBeta36=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              ISXB316=c()
              for(L in 1:nn){
                ISXB316[L]=meanX4(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }

              ISXB315=c()
              for(L in 1:nn){
                ISXB315[L]=meanX5(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }

              IBeta36[i]=as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean(ISXB316))-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(ISXB315))
            }
          }
        }
      }
      IBeta35[t]=sum(IBeta36)
    }
    IInfBeta35=sum(IBeta35,na.rm=T)
    IEndBeta35=ISusBeta312+IInfBeta35
    ######################################
    InfBeta3=-IEndBeta35-sum(SecondBeta3)+sum(FirstBeta3)
    ###########################################################################################
    ###########################################################################################
    ####################################### Info beta4 ##############################
    ####################################### Cov##############################

    Beta4S=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta4S1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Beta4S1[t]=-SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,mmm]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        Beta4S[i,L]=sum(Beta4S1)
      }
    }


    ###################### infectious period ###################

    Beta4I=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta4I1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                Beta4I1[t]=SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXstar2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        Beta4I[i,L]=sum(Beta4I1)
      }
    }


    Beta4IS=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        Beta4IS[i,L]=Beta4S[i,L]+Beta4I[i,L]
      }
    }

    SecondBeta4=c()
    for(i in 1:n){
      SecondBeta4[i]=mean(Beta4IS[i,]^2)
    }

    FirstBeta4=c()
    for(i in 1:n){
      FirstBeta4[i]=mean(Beta4IS[i,])^2
    }

    MBeta4=c()
    for(i in 1:n){
      MBeta4[i]=mean(Beta4IS[i,])
    }
    ###############################################################################################
    ###############################################################################################
    #################################beta4##################################
    ##################### infectious period ###################


    IBeta412=rep(0,time)
    for(t in 1:time){
      IBeta422=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){
              IBeta422[i]=-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXstar3(Pos,beta3,beta4,mmm,i))
            }
          }
        }
      }
      IBeta412[t]=sum(IBeta422)
    }
    ISusBeta412=sum(IBeta412)




    IBeta45=rep(0,time)
    for(t in 1:time){
      IBeta46=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              ISXB4=c()
              for(L in 1:nn){
                ISXB4[L]=meanXstar4(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }

              ISXB5=c()
              for(L in 1:nn){
                ISXB5[L]=meanXstar5(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }

              IBeta46[i]=as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean(ISXB4))-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(ISXB5))
            }
          }
        }
      }
      IBeta45[t]=sum(IBeta46)
    }
    IInfBeta45=sum(IBeta45,na.rm=T)
    IEndBeta45=ISusBeta412+IInfBeta45
    ######################################
    InfBeta4=-IEndBeta45-sum(SecondBeta4)+sum(FirstBeta4)
    ###########################################################################################
    ################################################# Info Delta #####################################
    #################################Cov #############################


    DeltaS=matrix(NA,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        DeltaS1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                DeltaS1[t]=SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        DeltaS[i,L]=sum(DeltaS1)
      }
    }



    ###################### infectious period ###################

    DeltaI=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        DeltaI1=rep(0,time)
        for(t in 1:time){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                DeltaI1[t]=-SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        DeltaI[i,L]=sum(DeltaI1)
      }
    }



    DeltaIS=matrix(0,n,nn)
    for(i in 1:n){
      for(L in 1:nn){
        DeltaIS[i,L]=DeltaS[i,L]+DeltaI[i,L]
      }
    }

    SecondDelta=c()
    for(i in 1:n){
      SecondDelta[i]=mean(DeltaIS[i,]^2)
    }

    FirstDelta=c()
    for(i in 1:n){
      FirstDelta[i]=mean(DeltaIS[i,])^2
    }

    MDelta=c()
    for(i in 1:n){
      MDelta[i]=mean(DeltaIS[i,])
    }
    ############################################################################################################
    ############################################################################################################
    ######################### Second Der...########################

    IS10=rep(0,time)
    for(t in 1:time){
      IS9=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if(Nlabel[i]==mmm){
            if (tau[i]>(t+1)|tau[i]== 0){
              IS9[i]=-as.numeric(SumHnew2(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean1(Pos,beta3,beta4,mmm,i))
            }
          }
        }
      }
      IS10[t]=sum(IS9)
    }
    ISusNew10=sum(IS10)

    ###################### infectious period ###################

    II12=rep(0,time)
    for(t in 1:time){
      II11=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if(Nlabel[i]==mmm){
            if(tau[i]==(t+1)){


              ISXD4=c()
              for(L in 1:nn){
                ISXD4[L]=mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              ISXD6=c()
              for(L in 1:nn){
                ISXD6[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              II11[i]=as.numeric(SumHnew2(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean(ISXD4))-as.numeric(SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(ISXD6))
            }
          }
        }
      }
      II12[t]=sum(II11)
    }
    IInfeNew12=sum(II12,na.rm=T)

    IEndNew12=ISusNew10+IInfeNew12
    ######################################
    InfDelta=-IEndNew12-sum(SecondDelta)+sum(FirstDelta)
    ###########################################################################################
    ###########################################################################################
    ###########################################################################################
    ###########################################Inf alpha & beta1################################################

    IAlphaBeta1=rep(0,time)
    for(t in 1:time){
      IAlphaBeta11=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              IAB1=c()
              for(L in 1:nn){
                IAB1[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              IAlphaBeta11[i]=-cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(IAB1)
            }
          }
        }
      }
      IAlphaBeta1[t]=sum(IAlphaBeta11)
    }

    DeriveAB1=sum(IAlphaBeta1,na.rm=T)+sum(MBeta1)

    AlphaBeta1=c()
    for(i in 1:n){
      AlphaBeta1[i]=mean(alphaIS[i,]*Beta1IS[i,])
    }


    InfAlBeta1=-DeriveAB1-sum(AlphaBeta1)+sum(MAlpha*MBeta1)


    ###########################################Inf alpha & beta2################################################

    IAlphaBeta2=rep(0,time)
    for(t in 1:time){
      IAlphaBeta21=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              IAB2=c()
              for(L in 1:nn){
                IAB2[L]=mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }


              IAlphaBeta21[i]=-cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(IAB2)
            }
          }
        }
      }
      IAlphaBeta2[t]=sum(IAlphaBeta21)
    }

    DeriveAB2=sum(IAlphaBeta2,na.rm=T)+sum(MBeta2)

    AlphaBeta2=c()
    for(i in 1:n){
      AlphaBeta2[i]=mean(alphaIS[i,]*Beta2IS[i,])
    }

    InfAlBeta2=-DeriveAB2-sum(AlphaBeta2)+sum(MAlpha*MBeta2)


    ###########################################Inf alpha & beta3################################################


    meanX6=function(Nlabel,Pos,Di,alpha1,beta1,beta3,beta2,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
      if (is.na(prob11)) {
        prob11 <- 0
      }

      if(prob11==0){SMeanX6=0}

      if(prob11!=0){

        SMeanX6=Pos[L,(2*mm+i)]*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
      }

      return(SMeanX6)
    }


    IAlphaBeta3=rep(0,time)
    for(t in 1:time){
      IAlphaBeta31=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              IABx5=c()
              for(L in 1:nn){
                IABx5[L]=meanX6(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }

              IAlphaBeta31[i]=-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(IABx5))
            }
          }
        }
      }
      IAlphaBeta3[t]=sum(IAlphaBeta31)
    }


    DeriveAB3=sum(IAlphaBeta3,na.rm=T)+sum(MBeta3)

    AlphaBeta3=c()
    for(i in 1:n){
      AlphaBeta3[i]=mean(alphaIS[i,]*Beta3IS[i,])
    }

    InfAlBeta3=-DeriveAB3-sum(AlphaBeta3)+sum(MAlpha*MBeta3)

    ###########################################Inf alpha & beta4################################################


    meanXstar6=function(Nlabel,Pos,Di,alpha1,beta1,beta3,beta2,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
      if (is.na(prob11)) {
        prob11 <- 0
      }

      if(prob11==0){SMean6=0}

      if(prob11!=0){

        SMean6=Pos[L,mmm]*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
      }

      return(SMean6)
    }


    IAlphaBeta4=rep(0,time)
    for(t in 1:time){
      IAlphaBeta41=rep(0,n)
      for(i in 1:n){
        for(mmm in 1:mm){
          if (Nlabel[i]==mmm){
            if(tau[i]==(t+1)){

              IAB5=c()
              for(L in 1:nn){
                IAB5[L]=meanXstar6(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }

              IAlphaBeta41[i]=-as.numeric(SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean(IAB5))
            }
          }
        }
      }
      IAlphaBeta4[t]=sum(IAlphaBeta41)
    }


    DeriveAB4=sum(IAlphaBeta4,na.rm=T)+sum(MBeta4)

    AlphaBeta4=c()
    for(i in 1:n){
      AlphaBeta4[i]=mean(alphaIS[i,]*Beta4IS[i,])
    }

    InfAlBeta4=-DeriveAB4-sum(AlphaBeta4)+sum(MAlpha*MBeta4)


    ###########################################Inf alpha & delta ################################################


    IAlphaDelta1=c()
    for(L in 1:nn){
      IAlphaDelta2=rep(0,time)
      for(t in 1:time){
        IAlphaDelta3=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                IAlphaDelta3[i]=SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        IAlphaDelta2[t]=sum(IAlphaDelta3)
      }
      IAlphaDelta1[L]=sum(IAlphaDelta2,na.rm=T)
    }

    DeriveAD=mean(IAlphaDelta1)+sum(MDelta)

    AlphaDelta=c()
    for(i in 1:n){
      AlphaDelta[i]=mean(alphaIS[i,]*DeltaIS[i,])
    }


    InfAlDelta=-DeriveAD-sum(AlphaDelta)+sum(MAlpha*MDelta)


    ###########################################Inf beta1 & delta ################################################


    Delta18=c()
    for(L in 1:nn){
      Delta8=rep(0,time)
      for(t in 1:time){
        Delta7=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                Delta7[i]=cov1[i]*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        Delta8[t]=sum(Delta7)
      }
      Delta18[L]=sum(Delta8)
    }



    ###################### infectious period ###################

    Delta10=c()
    for(L in 1:nn){
      Delta11=rep(0,time)
      for(t in 1:time){
        Delta12=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                Delta12[i]=-cov1[i]*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        Delta11[t]=sum(Delta12)
      }
      Delta10[L]=sum(Delta11,na.rm=T)
    }


    DeltaEnd=c()
    for(L in 1:nn){
      DeltaEnd[L]=Delta18[L]+Delta10[L]
    }






    Delta15=c()
    for(L in 1:nn){
      Delta16=rep(0,time)
      for(t in 1:time){
        Delta17=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                Delta17[i]=cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        Delta16[t]=sum(Delta17,na.rm=T)

      }
      Delta15[L]=sum(Delta16,na.rm=T)
    }

    DeriveBD=mean(Delta15)+mean(DeltaEnd)

    Beta1Delta=c()
    for(i in 1:n){
      Beta1Delta[i]=mean(Beta1IS[i,]*DeltaIS[i,])
    }

    InfBeta1Delta=-DeriveBD-sum(Beta1Delta)+sum(MBeta1*MDelta)


    ###########################################Inf beta2 & delta ################################################


    B2Delta18=c()
    for(L in 1:nn){
      B2Delta8=rep(0,time)
      for(t in 1:time){
        B2Delta7=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B2Delta7[i]=cov2[mmm]*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B2Delta8[t]=sum(B2Delta7)
      }
      B2Delta18[L]=sum(B2Delta8)
    }



    ###################### infectious period ###################

    B2Delta10=c()
    for(L in 1:nn){
      B2Delta11=rep(0,time)
      for(t in 1:time){
        B2Delta12=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                B2Delta12[i]=-cov2[mmm]*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B2Delta11[t]=sum(B2Delta12)
      }
      B2Delta10[L]=sum(B2Delta11,na.rm=T)
    }


    B2DeltaEnd=c()
    for(L in 1:nn){
      B2DeltaEnd[L]=B2Delta18[L]+B2Delta10[L]
    }





    B2Delta15=c()
    for(L in 1:nn){
      B2Delta16=rep(0,time)
      for(t in 1:time){
        B2Delta17=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                B2Delta17[i]=cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B2Delta16[t]=sum(B2Delta17)
      }
      B2Delta15[L]=sum(B2Delta16,na.rm=T)
    }

    DeriveB2D=mean(B2Delta15)+mean(B2DeltaEnd)

    Beta2Delta=c()
    for(i in 1:n){
      Beta2Delta[i]=mean(Beta2IS[i,]*DeltaIS[i,])
    }


    InfBeta2Delta=-DeriveB2D-sum(Beta2Delta)+sum(MBeta2*MDelta)



    ###########################################Inf beta3& delta ################################################

    B3Delta18=c()
    for(L in 1:nn){
      B3Delta8=rep(0,time)
      for(t in 1:time){
        B3Delta7=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B3Delta7[i]=SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,(2*mm+i)]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B3Delta8[t]=sum(B3Delta7)
      }
      B3Delta18[L]=sum(B3Delta8)
    }



    ###################### infectious period ###################

    meanX7=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
      if (is.na(prob11)) {
        prob11 <- 0
      }

      if(prob11==0){SMeanX7=0}

      if(prob11!=0){
        SMeanX7=Pos[L,(2*mm+i)]*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
      }

      return(SMeanX7)
    }



    B3Delta10=c()
    for(L in 1:nn){
      B3Delta11=rep(0,time)
      for(t in 1:time){
        B3Delta12=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                B3Delta12[i]=-SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanX2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B3Delta11[t]=sum(B3Delta12)
      }
      B3Delta10[L]=sum(B3Delta11,na.rm=T)
    }


    B3DeltaEnd=c()
    for(L in 1:nn){
      B3DeltaEnd[L]=B3Delta18[L]+B3Delta10[L]
    }





    B3Delta15=c()
    for(L in 1:nn){
      B3Delta16=rep(0,time)
      for(t in 1:time){
        B3Delta17=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                B3Delta17[i]=SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanX7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B3Delta16[t]=sum(B3Delta17)
      }
      B3Delta15[L]=sum(B3Delta16,na.rm=T)
    }

    DeriveB3D=mean(B3Delta15)+mean(B3DeltaEnd)

    Beta3Delta=c()
    for(i in 1:n){
      Beta3Delta[i]=mean(Beta3IS[i,]*DeltaIS[i,])
    }


    InfBeta3Delta=-DeriveB3D-sum(Beta3Delta)+sum(MBeta3*MDelta)


    ###########################################Inf beta4& delta ################################################

    B4Delta18=c()
    for(L in 1:nn){
      B4Delta8=rep(0,time)
      for(t in 1:time){
        B4Delta7=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B4Delta7[i]=SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,mmm]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B4Delta8[t]=sum(B4Delta7)
      }
      B4Delta18[L]=sum(B4Delta8)
    }



    ###################### infectious period ###################

    meanXstar7=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
      if (is.na(prob11)) {
        prob11 <- 0
      }

      if(prob11==0){SMeanXstar7=0}

      if(prob11!=0){
        SMeanXstar7=Pos[L,mmm]*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
      }

      return(SMeanXstar7)
    }



    B4Delta10=c()
    for(L in 1:nn){
      B4Delta11=rep(0,time)
      for(t in 1:time){
        B4Delta12=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                B4Delta12[i]=-SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXstar2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B4Delta11[t]=sum(B4Delta12)
      }
      B4Delta10[L]=sum(B4Delta11,na.rm=T)
    }


    B4DeltaEnd=c()
    for(L in 1:nn){
      B4DeltaEnd[L]=B4Delta18[L]+B4Delta10[L]
    }





    B4Delta15=c()
    for(L in 1:nn){
      B4Delta16=rep(0,time)
      for(t in 1:time){
        B4Delta17=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){

                B4Delta17[i]=SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*SumHnew1(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXstar7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B4Delta16[t]=sum(B4Delta17)
      }
      B4Delta15[L]=sum(B4Delta16,na.rm=T)
    }

    DeriveB4D=mean(B4Delta15)+mean(B4DeltaEnd)

    Beta4Delta=c()
    for(i in 1:n){
      Beta4Delta[i]=mean(Beta4IS[i,]*DeltaIS[i,])
    }


    InfBeta4Delta=-DeriveB4D-sum(Beta4Delta)+sum(MBeta4*MDelta)



    ########################################### beta1 & beta 2 ################################################################

    B1Beta21=c()
    for(L in 1:nn){
      B1NewBeta21=rep(0,time)
      for(t in 1:time){
        B1NewBeta22=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B1NewBeta22[i]=-cov1[i]*cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B1NewBeta21[t]=sum(B1NewBeta22)
      }
      B1Beta21[L]=sum(B1NewBeta21)
    }

    ###################### infectious period ###################

    B1InfBeta23=c()
    for(L in 1:nn){
      B1NewBeta23=rep(0,time)
      for(t in 1:time){
        B1NewBeta24=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                B1NewBeta24[i]=cov1[i]*cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*mean2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B1NewBeta23[t]=sum(B1NewBeta24)
      }
      B1InfBeta23[L]=sum(B1NewBeta23,na.rm=T)
    }




    New23Beta=c()
    for(L in 1:nn){
      New23Beta[L]=B1Beta21[L]+B1InfBeta23[L]
    }


    ######################### Second Der...########################
    ##################### infectious period ###################

    Beta1Beat2=c()
    for(L in 1:nn){
      IBeta25=rep(0,time)
      for(t in 1:time){
        IBeta26=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                IBeta26[i]=-cov1[i]*cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*mean3(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        IBeta25[t]=sum(IBeta26)
      }
      Beta1Beat2[L]=sum(IBeta25,na.rm=T)
    }


    DeriveBeta1Beat2=mean(Beta1Beat2)+mean(New23Beta)

    Beta1Beat21=c()
    for(i in 1:n){
      Beta1Beat21[i]=mean(Beta1IS[i,]*Beta2IS[i,])
    }

    InfBeta1Beat2=-DeriveBeta1Beat2-sum(Beta1Beat21)+sum(MBeta1*MBeta2)



    ########################################### beta1 & beta 3 ################################################################

    B1Bet1a31=c()
    for(L in 1:nn){
      B1NewBeta31=rep(0,time)
      for(t in 1:time){
        B1NewBeta32=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B1NewBeta32[i]=-cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,(2*mm+i)]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B1NewBeta31[t]=sum(B1NewBeta32)
      }
      B1Bet1a31[L]=sum(B1NewBeta31)
    }

    ###################### infectious period ###################

    B1NewInfBeta33=c()
    for(L in 1:nn){
      B1NewBeta33=rep(0,time)
      for(t in 1:time){
        B1NewBeta34=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                B1NewBeta34[i]=cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanX2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B1NewBeta33[t]=sum(B1NewBeta34)
      }
      B1NewInfBeta33[L]=sum(B1NewBeta33,na.rm=T)
    }




    B1EndBeta33=c()
    for(L in 1:nn){
      B1EndBeta33[L]=B1Bet1a31[L]+B1NewInfBeta33[L]
    }


    ######################### Second Der...########################
    ##################### infectious period ###################

    Beta1Beta3=c()
    for(L in 1:nn){
      IBeta35=rep(0,time)
      for(t in 1:time){
        IBeta36=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                IBeta36[i]=-cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*meanX7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        IBeta35[t]=sum(IBeta36)
      }
      Beta1Beta3[L]=sum(IBeta35,na.rm=T)
    }


    DeriveBeta1Beta3=mean(Beta1Beta3)+mean(B1EndBeta33)

    Beta1Beat31=c()
    for(i in 1:n){
      Beta1Beat31[i]=mean(Beta1IS[i,]*Beta3IS[i,])
    }

    InfBeta3Beta1=-DeriveBeta1Beta3-sum(Beta1Beat31)+sum(MBeta1*MBeta3)



    ########################################### beta1 & beta 4 ################################################################

    B1Bet1a41=c()
    for(L in 1:nn){
      B1NewBeta41=rep(0,time)
      for(t in 1:time){
        B1NewBeta42=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B1NewBeta42[i]=-cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,mmm]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B1NewBeta41[t]=sum(B1NewBeta42)
      }
      B1Bet1a41[L]=sum(B1NewBeta41)
    }

    ###################### infectious period ###################

    B1NewInfBeta43=c()
    for(L in 1:nn){
      B1NewBeta43=rep(0,time)
      for(t in 1:time){
        B1NewBeta44=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                B1NewBeta44[i]=cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXstar2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B1NewBeta43[t]=sum(B1NewBeta44)
      }
      B1NewInfBeta43[L]=sum(B1NewBeta43,na.rm=T)
    }




    B1EndBeta43=c()
    for(L in 1:nn){
      B1EndBeta43[L]=B1Bet1a41[L]+B1NewInfBeta43[L]
    }


    ######################### Second Der...########################
    ##################### infectious period ###################

    Beta1Beta4=c()
    for(L in 1:nn){
      IBeta45=rep(0,time)
      for(t in 1:time){
        IBeta46=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                IBeta46[i]=-cov1[i]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*meanXstar7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        IBeta45[t]=sum(IBeta46)
      }
      Beta1Beta4[L]=sum(IBeta45,na.rm=T)
    }


    DeriveBeta1Beta4=mean(Beta1Beta4)+mean(B1EndBeta43)

    Beta1Beat41=c()
    for(i in 1:n){
      Beta1Beat41[i]=mean(Beta1IS[i,]*Beta4IS[i,])
    }

    InfBeta4Beta1=-DeriveBeta1Beta4-sum(Beta1Beat41)+sum(MBeta1*MBeta4)


    ########################################### beta2 & beta 3 ################################################################

    B2Beta31=c()
    for(L in 1:nn){
      B2NewBeta31=rep(0,time)
      for(t in 1:time){
        B2NewBeta32=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B2NewBeta32[i]=-cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,(2*mm+i)]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B2NewBeta31[t]=sum(B2NewBeta32)
      }
      B2Beta31[L]=sum(B2NewBeta31)
    }

    ###################### infectious period ###################

    B2NewInfBeta33=c()
    for(L in 1:nn){
      B2NewBeta33=rep(0,time)
      for(t in 1:time){
        B2NewBeta34=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                B2NewBeta34[i]=cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanX2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B2NewBeta33[t]=sum(B2NewBeta34)
      }
      B2NewInfBeta33[L]=sum(B2NewBeta33,na.rm=T)
    }




    B2EndBeta33=c()
    for(L in 1:nn){
      B2EndBeta33[L]=B2Beta31[L]+B2NewInfBeta33[L]
    }


    ######################### Second Der...########################
    ##################### infectious period ###################

    Beta2Beta3=c()
    for(L in 1:nn){
      B2IBeta35=rep(0,time)
      for(t in 1:time){
        B2IBeta36=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                B2IBeta36[i]=-cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*meanX7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B2IBeta35[t]=sum(B2IBeta36)
      }
      Beta2Beta3[L]=sum(B2IBeta35,na.rm=T)
    }


    DeriveBeta2Beta3=mean(Beta2Beta3)+mean(B2EndBeta33)

    Beta2Beat31=c()
    for(i in 1:n){
      Beta2Beat31[i]=mean(Beta2IS[i,]*Beta3IS[i,])
    }

    InfBeta3Beta2=-DeriveBeta2Beta3-sum(Beta2Beat31)+sum(MBeta2*MBeta3)




    ########################################### beta2 & beta 4 ################################################################

    B2Beta41=c()
    for(L in 1:nn){
      B2NewBeta41=rep(0,time)
      for(t in 1:time){
        B2NewBeta42=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B2NewBeta42[i]=-cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,mmm]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B2NewBeta41[t]=sum(B2NewBeta42)
      }
      B2Beta41[L]=sum(B2NewBeta41)
    }

    ###################### infectious period ###################

    B2NewInfBeta43=c()
    for(L in 1:nn){
      B2NewBeta43=rep(0,time)
      for(t in 1:time){
        B2NewBeta44=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                B2NewBeta44[i]=cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXstar2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B2NewBeta43[t]=sum(B2NewBeta44)
      }
      B2NewInfBeta43[L]=sum(B2NewBeta43,na.rm=T)
    }




    B2EndBeta43=c()
    for(L in 1:nn){
      B2EndBeta43[L]=B2Beta41[L]+B2NewInfBeta43[L]
    }


    ######################### Second Der...########################
    ##################### infectious period ###################

    Beta2Beta4=c()
    for(L in 1:nn){
      B2IBeta45=rep(0,time)
      for(t in 1:time){
        B2IBeta46=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                B2IBeta46[i]=-cov2[mmm]*SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*meanXstar7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B2IBeta45[t]=sum(B2IBeta46)
      }
      Beta2Beta4[L]=sum(B2IBeta45,na.rm=T)
    }


    DeriveBeta2Beta4=mean(Beta2Beta4)+mean(B2EndBeta43)

    Beta2Beat41=c()
    for(i in 1:n){
      Beta2Beat41[i]=mean(Beta2IS[i,]*Beta4IS[i,])
    }

    InfBeta4Beta2=-DeriveBeta2Beta4-sum(Beta2Beat41)+sum(MBeta2*MBeta4)

    #######################################################################################################


    ########################################### beta3 & beta 4 ################################################################

    B3Beta41=c()
    for(L in 1:nn){
      B3NewBeta41=rep(0,time)
      for(t in 1:time){
        B3NewBeta42=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if (tau[i]>(t+1)|tau[i]== 0){
                B3NewBeta42[i]=-SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*Pos[L,mmm]*Pos[L,(2*mm+i)]*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
              }
            }
          }
        }
        B3NewBeta41[t]=sum(B3NewBeta42)
      }
      B3Beta41[L]=sum(B3NewBeta41)
    }

    ###################### infectious period ###################



    meanXX2=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
      if (is.na(prob11)) {
        prob11 <- 0
      }

      if(prob11==0){SMeanXX2=0}

      if(prob11!=0){
        SMeanXX2=Pos[L,(2*mm+i)]*Pos[L,mmm]*(1-prob11)/prob11*exp(beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])
      }

      return(SMeanXX2)
    }



    meanXX7=function(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda){


      dx1=rep(0,n)
      for(j in 1:n){
        if (Newlabel[j,mmm]!=0){
          if(tau[j]<=t & (tau[j]+lambda[j])>t & tau[j]!=0){
            dx1[j]=Di[i,j]^(-delta)
          }
        }
      }
      dx1=replace(dx1,is.infinite(dx1),0)
      dx=sum(dx1)

      prob11=1-exp(-exp(alpha1+beta1*cov1[i]+beta2*cov2[mmm]+beta3*Pos[L,(2*mm+i)]+beta4*Pos[L,mmm]+Pos[L,mmm+mm])*dx)
      if (is.na(prob11)) {
        prob11 <- 0
      }

      if(prob11==0){SMeanXX7=0}

      if(prob11!=0){
        SMeanXX7=Pos[L,(2*mm+i)]*Pos[L,mmm]*(1-prob11)/prob11^2*exp(2*beta3*Pos[L,(2*mm+i)]+2*beta4*Pos[L,mmm]+2*Pos[L,mmm+mm])
      }

      return(SMeanXX7)
    }






    B3NewInfBeta43=c()
    for(L in 1:nn){
      B3NewBeta43=rep(0,time)
      for(t in 1:time){
        B3NewBeta44=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if(Nlabel[i]==mmm){
              if(tau[i]==(t+1)){
                B3NewBeta44[i]=SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)*meanXX2(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B3NewBeta43[t]=sum(B3NewBeta44)
      }
      B3NewInfBeta43[L]=sum(B3NewBeta43,na.rm=T)
    }




    B3EndBeta43=c()
    for(L in 1:nn){
      B3EndBeta43[L]=B3Beta41[L]+B3NewInfBeta43[L]
    }


    ######################### Second Der...########################
    ##################### infectious period ###################

    Beta3Beta4=c()
    for(L in 1:nn){
      B3IBeta45=rep(0,time)
      for(t in 1:time){
        B3IBeta46=rep(0,n)
        for(i in 1:n){
          for(mmm in 1:mm){
            if (Nlabel[i]==mmm){
              if(tau[i]==(t+1)){


                B3IBeta46[i]=-SumH(Nlabel,Di,alpha1,beta1,beta2,delta,i,mmm,t,tau,lambda)^2*meanXX7(Nlabel,Pos,Di,alpha1,beta1,beta2,beta3,beta4,delta,mmm,i,t,L,tau,lambda)
              }
            }
          }
        }
        B3IBeta45[t]=sum(B3IBeta46)
      }
      Beta3Beta4[L]=sum(B3IBeta45,na.rm=T)
    }


    DeriveBeta3Beta4=mean(Beta3Beta4)+mean(B3EndBeta43)

    Beta3Beat41=c()
    for(i in 1:n){
      Beta3Beat41[i]=mean(Beta3IS[i,]*Beta4IS[i,])
    }

    InfBeta3Beta4=-DeriveBeta3Beta4-sum(Beta3Beat41)+sum(MBeta3*MBeta4)








    A1=c(InfAlpha,InfAlBeta1,InfAlBeta2,InfAlBeta3,InfAlBeta4,InfAlDelta)
    A2=c(InfAlBeta1,InfBeta1,InfBeta1Beat2,InfBeta3Beta1,InfBeta4Beta1,InfBeta1Delta)
    A3=c(InfAlBeta2,InfBeta1Beat2,InfBeta2,InfBeta3Beta2,InfBeta4Beta2,InfBeta2Delta)
    A4=c(InfAlBeta3,InfBeta3Beta1,InfBeta3Beta2,InfBeta3,InfBeta3Beta4,InfBeta3Delta)
    A5=c(InfAlBeta4,InfBeta4Beta1,InfBeta4Beta2,InfBeta3Beta4,InfBeta4,InfBeta4Delta)
    A6=c(InfAlDelta,InfBeta1Delta,InfBeta2Delta,InfBeta3Delta,InfBeta4Delta,InfDelta)
    ObservedInf=rbind(A1,A2,A3,A4,A5,A6)
    Variances=diag(solve(ObservedInf))
    Variances



    CInfAlpha=-IEndA5
    CInfBeta1=-IEndBeta15
    CInfBeta2=-IEndBeta25
    CInfBeta3=-IEndBeta35
    CInfBeta4=-IEndBeta45
    CInfDelta=-IEndNew12
    CInfAlBeta1=-DeriveAB1
    CInfAlBeta2=-DeriveAB2
    CInfAlBeta3=-DeriveAB3
    CInfAlBeta4=-DeriveAB4
    CInfAlDelta=-DeriveAD
    CInfBeta1Delta=-DeriveBD
    CInfBeta2Delta=-DeriveB2D
    CInfBeta3Delta=-DeriveB3D
    CInfBeta4Delta=-DeriveB4D
    CInfBeta1Beat2=-DeriveBeta1Beat2
    CInfBeta3Beta1=-DeriveBeta1Beta3
    CInfBeta4Beta1=-DeriveBeta1Beta4
    CInfBeta3Beta2=-DeriveBeta2Beta3
    CInfBeta4Beta2=-DeriveBeta2Beta4
    CInfBeta3Beta4=-DeriveBeta3Beta4

    CA1=c(CInfAlpha,CInfAlBeta1,CInfAlBeta2,CInfAlBeta3,CInfAlBeta4,CInfAlDelta)
    CA2=c(CInfAlBeta1,CInfBeta1,CInfBeta1Beat2,CInfBeta3Beta1,CInfBeta4Beta1,CInfBeta1Delta)
    CA3=c(CInfAlBeta2,CInfBeta1Beat2,CInfBeta2,CInfBeta3Beta2,CInfBeta4Beta2,CInfBeta2Delta)
    CA4=c(CInfAlBeta3,CInfBeta3Beta1,CInfBeta3Beta2,CInfBeta3,CInfBeta3Beta4,CInfBeta3Delta)
    CA5=c(CInfAlBeta4,CInfBeta4Beta1,CInfBeta4Beta2,CInfBeta3Beta4,CInfBeta4,CInfBeta4Delta)
    CA6=c(CInfAlDelta,CInfBeta1Delta,CInfBeta2Delta,CInfBeta3Delta,CInfBeta4Delta,CInfDelta)
    CObservedInf=rbind(CA1,CA2,CA3,CA4,CA5,CA6)
    CVariances=diag(solve(CObservedInf))
    CVariances
    ######################################Spatial parameters##########################################
    QLoglik11=function(par){

      sigma1=par[[1]]
      lambda1=par[[2]]
      Log1=c()
      for(L in 1:nn){
        Log1[L]=dmvnorm(Pos[L,(mm+1):(2*mm)], rep(0,mm), sigma1^2*solve((lambda1*D+(1-lambda1)*I)), log = T)
      }
      LLL=mean(Log1)
      return(LLL)
    }
    Hes=hessian(QLoglik11,c(sigma1,lambda1))

    VarSpatial=diag(solve(-Hes))

    result=list(CVariances=CVariances,VarSpatial=VarSpatial,ObservedInf=ObservedInf,Variances=Variances)
    result
  }

  ################################################################################

  CovMEA=matrix(0,ITER,mm)
  CovMEInd=matrix(0,ITER,n)
  GePHI=matrix(0,ITER,mm)

  for(iter in 1:ITER){

    I=diag(mm)
    mu=rep(0,mm)
    Sigma0=sigma0^2*solve((lambda0*D+(1-lambda0)*I))
    det(Sigma0)
    GePHI[iter,]=mvrnorm(1, mu, Sigma0, tol = 1e-6)
    phi=GePHI[iter,]

################################################################################################

    LA=c()
    ss=1
    est0=estfunNaiv(Nlabel,phi,Di,alpha0,beta10,beta20,beta30,beta40,delta0,sigma0,lambda0,D,lambda)
    alpha1=est0$alpha1
    delta=est0$delta
    lambda1=est0$lambda1
    sigma1=est0$sigma1
    beta1=est0$beta1
    beta2=est0$beta2
    beta3=est0$beta3
    beta4=est0$beta4
    Pos1=est0$Pos1

    if(is.na(alpha1)|is.na(delta)|is.na(lambda1)|is.na(sigma1)|is.na(beta1)|is.na(beta2)|is.na(beta3)|is.na(beta4)){
      LA[1]=NA
    }
    if(!is.na(alpha1)&&!is.na(delta)&&!is.na(lambda1)&&!is.na(sigma1)&&!is.na(beta1)&&!is.na(beta2)&&!is.na(beta3)&&!is.na(beta4)){

      LA[1]=LogNaiv(Nlabel,Pos1,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,D,lambda)
    }

    repeat{
      ss=ss+1
      if(is.na(LA[1])){
        out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,ss=NA)
        break
      }

      if(ss>25){
        alpha1=NA
        delta=NA
        sigma1=NA
        lambda1=NA
        out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,ss=NA)
        break
      }


      if(!is.na(alpha1)&&!is.na(delta)&&!is.na(lambda1)&&!is.na(sigma1)&&!is.na(beta1)&&!is.na(beta2)&&!is.na(beta3)&&!is.na(beta4)){
        if(sigma1<0 |lambda1>1|lambda1<0|delta<0){
          alpha1=NA
          delta=NA
          sigma1=NA
          lambda1=NA
          out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,ss=ss)
          break
        }
      }
      E0=c(alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1)

      LA


      est=estfunNaiv(Nlabel,phi,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,D,lambda)
      alpha1=est$alpha1
      delta=est$delta
      lambda1=est$lambda1
      sigma1=est$sigma1
      beta1=est$beta1
      beta2=est$beta2
      beta3=est$beta3
      beta4=est$beta4

      Pos1=est$Pos1

      if(!is.na(alpha1)&&!is.na(delta)&&!is.na(lambda1)&&!is.na(sigma1)&&!is.na(beta1)&&!is.na(beta2)&&!is.na(beta3)&&!is.na(beta4)){

        if(sigma1<0 |lambda1>1|lambda1<0|delta<0){
          alpha1=NA
          delta=NA
          sigma1=NA
          lambda1=NA
          out1=list(alpha1=NA,delta=NA,sigma1=NA,lambda1=NA,beta1=NA,ss=ss)


          break
        }

        if(sigma1>=0 & lambda1<=1 & lambda1>=0){

          E1=c(alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1)

          Dist=sqrt(sum((E0-E1)^2))

          LA[ss]=LogNaiv(Nlabel,Pos1,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,D,lambda)

          out1=list(LA=LA,Pos1=Pos1,alpha1=alpha1,beta1=beta1,beta2=beta2,beta3=beta3,beta4=beta4,delta=delta,sigma1=sigma1,lambda1=lambda1,ss=ss,E1=E1,E0=E0)

          if(!abs((LA[ss]-LA[ss-1])/LA[ss])>eps){
            break
          }
        }
      }

      if(is.na(alpha1)|is.na(delta)|is.na(lambda1)|is.na(sigma1)|is.na(beta1)|is.na(beta2)|is.na(beta3)|is.na(beta4)){

        out1=list(alpha1=NA,beta1=NA,beta2=NA,beta3=NA,beta4=NA,delta=NA,sigma1=NA,lambda1=NA)

        break
      }

    }

    NaivPos=Pos1

    NaivEstpar=c(alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1)
    NaivFinalInf=ObsInfNaiv(Nlabel,NaivPos,Di,alpha1,beta1,beta2,beta3,beta4,delta,sigma1,lambda1,lambda)



    if(is.na(sigma1)){
      NaivSimPAR[iter,]=NA
      NaivSimVarPar[iter,]=NA
      NaivCSimVarPar[iter,]=NA
      NaivSimInfPar[iter,]=NA
      NaivSimVarSpatial[iter,]=NA
      NaivSimPos[[iter]]=NA
    }
    if(!is.na(sigma1)){
      NaivSimPAR[iter,]=NaivEstpar
      NaivSimVarPar[iter,]=NaivFinalInf$Variances
      NaivSimInfPar[iter,]=diag(NaivFinalInf$ObservedInf)
      NaivSimVarSpatial[iter,]=(NaivFinalInf$VarSpatial)
      NaivSimPos[[iter]]=NaivPos
      NaivCSimVarPar[iter,]=NaivFinalInf$CVariances

    }

  }
  result=list(NaivSimPAR=NaivSimPAR,NaivSimVarPar=NaivSimVarPar,NaivSimVarSpatial=NaivSimVarSpatial,NaivSimPos=NaivSimPos,NaivCSimVarPar=NaivCSimVarPar)

}
