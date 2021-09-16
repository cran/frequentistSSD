sample_size=function(kappa,mA,hr,ta,tf,diff,P)
{
  allocation <- .5
  lambda1=log(2)^(1/kappa)/mA
  lambda2=lambda1*hr^(1/kappa)
  tau=ta+tf
  S1=function(t){exp(-(lambda1*t)^kappa)}
  h1=function(t){kappa*lambda1^kappa*t^(kappa-1)}
  S2=function(t){exp(-(lambda2*t)^kappa)}
  h2=function(t){kappa*lambda2^kappa*t^(kappa-1)}
  G=function(t){1-punif(t, tf, tau)}
  f1=function(t){S1(t)*h1(t)*G(t)}
  f2=function(t){S2(t)*h2(t)*G(t)}
  theta0=-log(hr)
  p0=(integrate(f1, 0,tau)$value+integrate(f2, 0,tau)$value)/2
  zP=qnorm(P)
  n=ceiling((zP+diff/sqrt(2))^2/(allocation*(1-allocation)*theta0^2*p0)/2)
  message("For each arm, required sample size:")
  return(n)
}
