optimal.KJ <- function(shape,S0,x0,hr,tf,rate,alpha,beta,dist,data)
{
  calculate_alpha<-function(c2, c1, rho0){
    fun1<-function(z, c1, rho0){
      f<-dnorm(z)*pnorm((rho0*z-c1)/sqrt(1-rho0^2))
      return(f)
    }
    alpha<-integrate(fun1, lower= c2, upper= Inf, c1, rho0)$value
    return(alpha)
  }

  calculate_power<-function(cb, cb1, rho1){
    fun2<-function(z, cb1, rho1){
      f<-dnorm(z)*pnorm((rho1*z-cb1)/sqrt(1-rho1^2))
      return(f)
    }
    pwr<-integrate(fun2,lower=cb,upper=Inf,cb1=cb1,rho1=rho1)$value
    return(pwr)
  }

  fct<-function(zi, ceps=0.0001,alphaeps=0.0001,nbmaxiter=100,dist){
    ta<-as.numeric(zi[1])
    t1<-as.numeric(zi[2])
    c1<-as.numeric(zi[3])

    if (dist=="WB"){
      f0=function(u){(shape/scale0)*(u/scale0)^(shape-1)*exp(-(u/scale0)^shape)}
      s0=function(u){exp(-(u/scale0)^shape)}
      h0=function(u){(shape/scale0)*(u/scale0)^(shape-1)}
      H0=function(u){(u/scale0)^shape}
      s=function(b, u){exp(-b*(u/scale0)^shape)}
      h=function(b,u){b*h0(u)}
      H=function(b,u){b*H0(u)}
      scale0=x0/(-log(S0))^(1/shape)
      scale1=hr
    }


    if (dist=="LN"){
      s0=function(u){1-plnorm(u,scale0,shape)}
      f0=function(u){dlnorm(u,scale0,shape)}
      h0=function(u){f0(u)/s0(u)}
      H0=function(u){-log(s0(u))}
      s=function(b,u) {s0(u)^b}
      h=function(b,u){b*f0(u)/s0(u)}
      H=function(b,u){-b*log(s0(u))}
      scale0=log(x0)-shape*qnorm(1-S0)
      scale1=hr
    }

    if (dist=="LG"){
      s0=function(u){1/(1+(u/scale0)^shape)}
      f0=function(u){(shape/scale0)*(u/scale0)^(shape-1)/(1+(u/scale0)^shape)^2}
      h0=function(u){f0(u)/s0(u)}
      H0=function(u){-log(s0(u))}
      s=function(b,u) {s0(u)^b}
      h=function(b,u){b*f0(u)/s0(u)}
      H=function(b,u){-b*log(s0(u))}
      scale0=x0/(1/S0-1)^(1/shape)
      scale1=hr
    }

    if (dist=="GM"){
      s0=function(u){1-pgamma(u,shape,scale0)}
      f0=function(u){dgamma(u,shape,scale0)}
      h0=function(u){f0(u)/s0(u)}
      H0=function(u){-log(s0(u))}
      s=function(b,u) {s0(u)^b}
      h=function(b,u){b*f0(u)/s0(u)}
      H=function(b,u){-b*log(s0(u))}
      root0=function(t){1-pgamma(x0,shape,t)-S0}
      scale0=uniroot(root0,c(0,10))$root
      scale1=hr
    }

    G=function(t){1-punif(t, tf, ta+tf)}
    g0=function(t){s(scale1,t)*h0(t)*G(t)}
    g1=function(t){s(scale1,t)*h(scale1,t)*G(t)}
    g00=function(t){s(scale1,t)*H0(t)*h0(t)*G(t)}
    g01=function(t){s(scale1,t)*H0(t)*h(scale1,t)*G(t)}
    p0=integrate(g0, 0, ta+tf)$value
    p1=integrate(g1, 0, ta+tf)$value
    p00=integrate(g00, 0, ta+tf)$value
    p01=integrate(g01, 0, ta+tf)$value
    sigma2.1=p1-p1^2+2*p00-p0^2-2*p01+2*p0*p1
    sigma2.0=p0
    om=p0-p1

    G1=function(t){1-punif(t, t1-ta, t1)}
    g0=function(t){s(scale1,t)*h0(t)*G1(t)}
    g1=function(t){s(scale1,t)*h(scale1,t)*G1(t)}
    g00=function(t){s(scale1,t)*H0(t)*h0(t)*G1(t)}
    g01=function(t){s(scale1,t)*H0(t)*h(scale1,t)*G1(t)}
    p0=integrate(g0, 0, ta+tf)$value
    p1=integrate(g1, 0, ta+tf)$value
    p00=integrate(g00, 0, ta+tf)$value
    p01=integrate(g01, 0, ta+tf)$value
    sigma2.11=p1-p1^2+2*p00-p0^2-2*p01+2*p0*p1
    sigma2.01=p0
    om1=p0-p1

    q1=function(t){s0(t)*h0(t)*G1(t)}
    q=function(t){s0(t)*h0(t)*G(t)}
    v1=integrate(q1, 0, ta+tf)$value
    v=integrate(q, 0, ta+tf)$value
    rho0=sqrt(v1/v)
    rho1<-sqrt(sigma2.11/sigma2.1)

    cL<-(-10)
    cU<-(10)
    alphac<-100
    iter<-0
    while ((abs(alphac-alpha)>alphaeps|cU-cL>ceps)&iter<nbmaxiter){
      iter<-iter+1
      c<-(cL+cU)/2
      alphac<-calculate_alpha(c, c1, rho0)
      if (alphac>alpha) {
        cL<-c
      } else {
        cU<-c
      }
    }

    cb1<-sqrt((sigma2.01/sigma2.11))*(c1-(om1*sqrt(rate*t1)/sqrt(sigma2.01)))
    cb<-sqrt((sigma2.0/sigma2.1))*(c-(om*sqrt(rate*ta))/sqrt(sigma2.0))
    pwrc<-calculate_power(cb=cb, cb1=cb1, rho1=rho1)
    res<-c(cL, cU, alphac, 1-pwrc, rho0, rho1, cb1, cb)
    return(res)
  }

  c1<-0 ; rho0<-0; cb1<-0; rho1<-0 ; hz<-c(0,0);
  ceps<-0.001;alphaeps<-0.001;nbmaxiter<-100

  Duration=function(shape,S0,x0,hr,tf,rate,alpha,beta,dist="WB")
  {
    if (dist=="WB"){
      f0=function(u){(shape/scale0)*(u/scale0)^(shape-1)*exp(-(u/scale0)^shape)}
      s0=function(u){exp(-(u/scale0)^shape)}
      h0=function(u){(shape/scale0)*(u/scale0)^(shape-1)}
      H0=function(u){(u/scale0)^shape}
      s=function(b, u){exp(-b*(u/scale0)^shape)}
      h=function(b,u){b*h0(u)}
      H=function(b,u){b*H0(u)}
      scale0=x0/(-log(S0))^(1/shape)
      scale1=hr
    }


    if (dist=="LN"){
      s0=function(u){1-plnorm(u,scale0,shape)}
      f0=function(u){dlnorm(u,scale0,shape)}
      h0=function(u){f0(u)/s0(u)}
      H0=function(u){-log(s0(u))}
      s=function(b,u) {s0(u)^b}
      h=function(b,u){b*f0(u)/s0(u)}
      H=function(b,u){-b*log(s0(u))}
      scale0=log(x0)-shape*qnorm(1-S0)
      scale1=hr
    }

    if (dist=="LG"){
      s0=function(u){1/(1+(u/scale0)^shape)}
      f0=function(u){(shape/scale0)*(u/scale0)^(shape-1)/(1+(u/scale0)^shape)^2}
      h0=function(u){f0(u)/s0(u)}
      H0=function(u){-log(s0(u))}
      s=function(b,u) {s0(u)^b}
      h=function(b,u){b*f0(u)/s0(u)}
      H=function(b,u){-b*log(s0(u))}
      scale0=x0/(1/S0-1)^(1/shape)
      scale1=hr
    }

    if (dist=="GM"){
      s0=function(u){1-pgamma(u,shape,scale0)}
      f0=function(u){dgamma(u,shape,scale0)}
      h0=function(u){f0(u)/s0(u)}
      H0=function(u){-log(s0(u))}
      s=function(b,u) {s0(u)^b}
      h=function(b,u){b*f0(u)/s0(u)}
      H=function(b,u){-b*log(s0(u))}
      root0=function(t){1-pgamma(x0,shape,t)-S0}
      scale0=uniroot(root0,c(0,10))$root
      scale1=hr
    }

    ### need calculate ta for given power for single stage design #####
    root=function(ta){
      tau=ta+tf
      G=function(t){1-punif(t, tf, tau)}
      g0=function(t){s(scale1,t)*h0(t)*G(t)}
      g1=function(t){s(scale1,t)*h(scale1,t)*G(t)}
      g00=function(t){s(scale1,t)*H0(t)*h0(t)*G(t)}
      g01=function(t){s(scale1,t)*H0(t)*h(scale1,t)*G(t)}
      p0=integrate(g0, 0, tau)$value
      p1=integrate(g1, 0, tau)$value
      p00=integrate(g00, 0, tau)$value
      p01=integrate(g01, 0, tau)$value
      s1=sqrt(p1-p1^2+2*p00-p0^2-2*p01+2*p0*p1)
      s0=sqrt(p0)
      om=p0-p1
      rate*ta-(s0*qnorm(1-alpha)+s1*qnorm(1-beta))^2/om^2
    }
    tasingle=uniroot(root, lower=0, upper=50*tf)$root
    nsingle<-ceiling(tasingle*rate)
    tasingle=ceiling(tasingle)
    ans=list(nsingle=nsingle, tasingle=tasingle)
    return(ans)
  }

  ans=Duration(shape,S0,x0,hr,tf,rate,alpha,beta,dist="WB")
  tasingle<-ans$tasingle
  nsingle=ans$nsingle

  Single_stage<-data.frame(nsingle=nsingle,tasingle=tasingle,
                           csingle=qnorm(1-alpha))
  atc0<-data.frame(n=nsingle, t1=tasingle, c1=0.25)
  nbpt<-11
  pascote<-1.26
  cote<-1*pascote
  c1.lim<-nbpt*c(-1,1)
  EnH0<-10000
  iter<-0

  while (iter<nbmaxiter & diff(c1.lim)/nbpt>0.001){
    iter<-iter+1
    cat("iter=",iter,"&EnH0=",round(EnH0, 2),"& Dc1/nbpt=",
        round(diff(c1.lim)/nbpt,5),"\n",sep="")

    if (iter%%2==0) nbpt<-nbpt+1
    cote<-cote/pascote
    n.lim<-atc0$n+c(-1, 1)*nsingle*cote
    t1.lim<-atc0$t1+c(-1, 1)*tasingle*cote
    c1.lim<-atc0$c1+c(-1, 1)*cote
    ta.lim<-n.lim/rate
    t1.lim<-pmax(0, t1.lim)
    n<-seq(n.lim[1], n.lim[2], l=nbpt)
    n<-ceiling(n)
    n<-unique(n)
    ta<-n/rate
    t1<-seq(t1.lim[1], t1.lim[2], l=nbpt)
    c1<-seq(c1.lim[1], c1.lim[2], l=nbpt)
    ta<-ta[ta>0]
    t1<-t1[t1>=0.2*tasingle & t1<=1.2*tasingle]
    z<-expand.grid(list(ta=ta, t1=t1, c1=c1))
    z<-z[z$ta>z$t1,]
    nz<-dim(z)[1]
    z$pap<-pnorm(z$c1)
    z$eta<-z$ta-pmax(0, z$ta-z$t1)*z$pap
    z$enh0<-z$eta*rate
    z<-z[z$enh0<=EnH0,]
    nz1<-dim(z)[1]
    resz<-t(apply(z,1,fct,ceps=ceps,alphaeps=alphaeps,nbmaxiter=nbmaxiter,
                  dist=dist))
    resz<-as.data.frame(resz)
    names(resz)<-c("cL","cU","alphac","betac","rho0","rho1","cb1","cb")
    r<-cbind(z, resz)
    r$pap<-pnorm(r$c1)
    r$eta<-r$ta-pmax(0, r$ta-r$t1)*r$pap
    r$etar<-r$eta*rate
    r$tar<-r$ta*rate
    r$c<-r$cL+(r$cU-resz$cL)/2
    r$diffc<-r$cU-r$cL
    r$diffc<-ifelse(r$diffc<=ceps, 1, 0)
    r<-r[1-r$betac>=1-beta,]
    r<-r[order(r$enh0),]
    r$n<-r$ta*rate

    if (dim(r)[1]>0) {
      atc<-r[, c("ta", "t1", "c1", "n")][1,]

      r1<-r[1,]
      if (r1$enh0<EnH0) {
        EnH0<-r1$enh0
        atc0<-atc
      }
    } else {
      atc<-data.frame(ta=NA, t1=NA, c1=NA, n=NA)
    }
    atc$iter<-iter
    atc$enh0<-r$enh0[1]
    atc$EnH0<-EnH0
    atc$tai<-ta.lim[1]
    atc$tas<-ta.lim[2]
    atc$ti<-t1.lim[1]
    atc$ts<-t1.lim[2]
    atc$ci<-c1.lim[1]
    atc$cs<-c1.lim[2]
    atc$cote<-cote
    if (iter==1) {
      atcs<-atc
    } else {
      atcs<-rbind(atcs, atc)
    }
  }
  atcs$i<-1:dim(atcs)[1]
  a<-atcs[!is.na(atcs$n),]

  p<-t(apply(a,1,fct,ceps=ceps,alphaeps=alphaeps,nbmaxiter=nbmaxiter,
             dist=dist))
  p<-as.data.frame(p)
  names(p)<-c("cL","cU","alphac","betac","rho0","rho1","cb1","cb")
  p$i<-a$i

  res<-merge(atcs, p, by="i", all=T)
  res$pap<-round(pnorm(res$c1),4)  ## stopping prob of stage 1 ##
  res$ta<-res$n/rate
  res$Enh0<-(res$ta-pmax(0, res$ta-res$t1)*res$pap)*rate
  res$diffc<-res$cU-res$cL
  res$c<-round(res$cL+(res$cU-res$cL)/2,4)
  res$diffc<-ifelse(res$diffc<=ceps, 1, 0)
  res<-res[order(res$enh0),]
  des<-round(res[1,],4)

  param<-data.frame(shape=shape,S0=S0,hr=hr,alpha=alpha,beta=beta,
                    rate=rate, x0=x0, tf=tf)
  Two_stage<- data.frame(n1= ceiling(des$t1*param$rate), c1=des$c1,
                         n=ceiling(des$ta*param$rate), c=des$c, t1=des$t1,
                         MTSL=des$ta+param$tf, ESS=des$Enh0, PS=des$pap)
  param<-data.frame(S0=S0,hr=hr,rate=rate)
  DESIGN<-list(param=param,Two_stage=Two_stage)
  return(DESIGN)
}
