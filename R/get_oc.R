get_oc=function(shape, m0, mA, hr, frac, ta, tf, c1, c, diff, n, nsim, seed = 2483)
{
  kappa=shape
  m2=mA/hr
  lambda0=log(2)/m0^kappa
  rho0=lambda0^(1/kappa)
  scale0=1/rho0
  H0=function(shape, scale, t){(t/scale)^shape}

  lambda1=log(2)/mA^kappa
  rho1=lambda1^(1/kappa)   # 1/rho is Weibull scale parameter ##
  lambda2=log(2)/m2^kappa
  rho2=lambda2^(1/kappa)
  scale1=1/rho1
  scale2=1/rho2
  scale=c(scale1, scale2)

  s=0
  set.seed(seed)
  n1=ceiling(n*frac)
  n2=ceiling(n*(1-frac))
  tau=ta+tf

  outcome2<-No.success<-n.Subj<-matrix(999, ncol=2, nrow=nsim)

  for (i in 1:nsim)
  {
    for (a in 1:2){
      w=rweibull(n, shape, scale[a])
      u=runif(n, 0, ta)
      x=pmax(0, pmin(w,tau-u))
      delta = as.numeric(w<tau-u)
      O1=sum(delta[1:n1])
      M1=H0(shape,scale0,x[1:n1])
      E1=sum(M1)
      Z1=(E1-O1)/sqrt(E1)

      O=sum(delta)
      M=H0(shape,scale0,x)
      E=sum(M)
      Z=(E-O)/sqrt(E)

      outcome1<-ifelse(Z1>c1, Z, NA)
      No.success[i,a]<-Z
      outcome2[i,a]<-ifelse(outcome1>c, 1, 0)
      n.Subj[i,a] <- ifelse(Z1>c1, n, n1)
    }
  }
  Outcome<-apply(outcome2, 1, sum, na.rm=T)
  # Outcome=0,1,2 both arms nagtive, one arm positive, both arms positive #
  # positve means 2-stage results are positive #
  Outcome[is.na(Outcome)]<-0

  Prob.neg    <-length(Outcome[Outcome==0])/nsim
  Prob.pos    <-length(Outcome[Outcome==2])/nsim
  Prob.negpos <-length(Outcome[Outcome==1])/nsim
  Prob.ArmA<-sum(outcome2[,1], na.rm=TRUE)/nsim
  Prob.ArmB<-sum(outcome2[,2], na.rm=TRUE)/nsim

  mean.Subj<-apply(n.Subj,2,mean, na.rm=T)

  ### After 1st segment
  Prob.select.ArmA<-sum(outcome2[,1][outcome2[,2]==0 | is.na(outcome2[,2])],na.rm=TRUE)/nsim  #only A is selected
  Prob.select.ArmB<-sum(outcome2[,2][outcome2[,1]==0 | is.na(outcome2[,1])],na.rm=TRUE)/nsim  #only B is selected
  Prob.NoArm<-Prob.neg

  ### Selecting an Arm when both arms are positive (2nd Segment)

  No.success.BothArms.select<-No.success[Outcome==2,]

  SSD.SelectArm<-function(x, diff, MOD=FALSE)
  {
    NoArm<-ArmA<-ArmB<-NA
    ArmB<-ifelse(x[2]-x[1] > diff,1,0)
    ArmA<-ifelse(x[2]-x[1] < -diff,1,0)

    ## for equivelence margin
    if(abs(x[2]-x[1]) <= diff) {
      NoArm <- 1
      # ArmA<-ifelse(runif(1,0,1)<0.5,0,1)  # random select a arm for ties A and B
      # ArmB<-1-ArmA
    }
    return(list(ArmA, ArmB, NoArm))
  }

  ### Original SSD with diff=0

  if( length(No.success[Outcome==2,]) > 1 ) {
    SSD.SelectArm.2ndSeg<-matrix(unlist(
      apply(No.success[Outcome==2,],1,SSD.SelectArm, diff)
    ),ncol=3,byrow=T)}   else {SSD.SelectArm.2ndSeg <- matrix(NA,ncol=3)  }

  ProbArmA.2ndSeg<-sum(SSD.SelectArm.2ndSeg[,1],na.rm=TRUE)/nsim
  ProbArmB.2ndSeg<-sum(SSD.SelectArm.2ndSeg[,2],na.rm=TRUE)/nsim
  ProbNoArm.2ndSeg<-sum(SSD.SelectArm.2ndSeg[,3],na.rm=TRUE)/nsim

  Overall.ArmA<-Prob.select.ArmA + ProbArmA.2ndSeg
  Overall.ArmB<-Prob.select.ArmB + ProbArmB.2ndSeg
  Overall.NoArm<-Prob.NoArm + ProbNoArm.2ndSeg

  soln<-data.frame("n"=n,
                   "SSD Arm A"=Overall.ArmA, "SSD Arm B"= Overall.ArmB, "SSD No Arm"=Overall.NoArm,
                   "diff"=diff,"Mean N Arm A"=mean.Subj[1],"Mean N Arm B"=mean.Subj[2])
  return(soln)
}
