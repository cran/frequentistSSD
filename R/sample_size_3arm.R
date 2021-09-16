sample_size_3arm=function(kappa, m0, mA, mB, delta, ta, tf, P, diff)
{
  lambda0=log(2)/m0^kappa
  lambdaA=log(2)/mA^kappa
  lambdaB=log(2)/mB^kappa
  mT=mA/delta^(1/kappa)
  lambdaT=log(2)/mT^kappa
  deltaA=lambdaA/lambda0
  deltaB=lambdaB/lambda0
  deltaT=lambdaT/lambda0

  tau=ta+tf
  S1=function(t){exp(-lambdaA*t^kappa)}
  h1=function(t){kappa*lambdaA*t^(kappa-1)}
  S2=function(t){exp(-lambdaB*t^kappa)}
  h2=function(t){kappa*lambdaB*t^(kappa-1)}
  S3=function(t){exp(-lambdaT*t^kappa)}
  h3=function(t){kappa*lambdaT*t^(kappa-1)}
  G=function(t){1-punif(t, tf, tau)}
  f1=function(t){S1(t)*h1(t)*G(t)}
  f2=function(t){S2(t)*h2(t)*G(t)}
  f3=function(t){S3(t)*h3(t)*G(t)}

  p=(integrate(f1, 0,tau)$value+integrate(f2, 0,tau)$value+integrate(f3, 0,tau)$value)/3

  root=function(n){
    muA=-sqrt(n*p)*log(deltaA)
    muB=-sqrt(n*p)*log(deltaB)
    muT=-sqrt(n*p)*log(deltaT)
    a=muT-muA
    b=muT-muB

    integrand = function(x, a, b){pnorm(a+x-diff)*pnorm(b+x-diff)*dnorm(x)}
    ans=P-integrate(integrand, lower=-Inf, upper=Inf, a, b)$value
  }
  n=uniroot(root, lower=1, upper=9999)$root
  message("For each arm, required sample size:")
  return(ceiling(n))
}
