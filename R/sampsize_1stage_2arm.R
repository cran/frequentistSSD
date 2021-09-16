sampsize_1stage_2arm=function(kappa,alpha,beta,m0,m1,ta,tf,delta){
  f=function(s){pnorm(s)*dnorm(s)}
  root1=function(c){alpha-integrate(f,c,Inf)$value}
  c=uniroot(root1, lower=0, upper=999)$root

  lambda0=log(2)/m0^kappa
  lambda1=log(2)/m1^kappa
  lambda2=lambda1*delta
  delta1=lambda1/lambda0
  delta2=lambda2/lambda0
  S0=function(t){exp(-lambda0*t^kappa)}
  h0=function(t){kappa*lambda0*t^(kappa-1)}
  S1=function(t){exp(-lambda1*t^kappa)}
  h1=function(t){kappa*lambda1*t^(kappa-1)}
  S2=function(t){exp(-lambda2*t^kappa)}
  h2=function(t){kappa*lambda2*t^(kappa-1)}
  tau=ta+tf
  G=function(t){1-punif(t, tf, tau)}
  f0=function(t){S0(t)*h0(t)*G(t)}
  f1=function(t){S1(t)*h1(t)*G(t)}
  f2=function(t){S2(t)*h2(t)*G(t)}
  p0=integrate(f0, 0,tau)$value
  p1=integrate(f1, 0,tau)$value
  p2=integrate(f2, 0,tau)$value
  root2=function(n){
    mu1=-sqrt(n*(p1+p0)/2)*log(delta1)
    mu2=-sqrt(n*(p2+p0)/2)*log(delta2)
    g=function(s){pnorm(s-mu1,mean=0,sd=1)*dnorm(s-mu2, mean=0,sd=1)}
    1-beta-integrate(g,c,Inf)$value
  }
  n=uniroot(root2, lower=10, upper=500)$root
  message("For each arm, required sample size:")
  return(ceiling(n))
}
