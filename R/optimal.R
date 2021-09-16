optimal <- function(shape, m0, ms, tf, ta, tot_size, dist){
  x0<-m0
  lambda0 <- log(2)^(1/shape)/m0
  S0 <- exp(-(lambda0*m0)^shape)
  lambda1 <- log(2)^(1/shape)/ms
  hr <- (lambda1/lambda0)^shape
  rate <- tot_size/ta
  alpha_seq <- seq(0.05, 0.2, 0.01)
  beta_seq <- seq(0.05, 0.2, 0.01)
  for(i in 1:length(alpha_seq)){
    for(j in 1:length(beta_seq)){
      design <- optimal.KJ(shape,S0,x0,hr,tf,rate,alpha_seq[i],beta_seq[j],dist="WB")
      design$Two_stage$alpha <- alpha_seq[i]
      design$Two_stage$beta <- beta_seq[j]
      design$Two_stage <- design$Two_stage[, c(9, 10, 1:8)]
      if(design$Two_stage$n == tot_size){
        return(design)
      }
    }
  }
  beta_seq_two <- seq(0.06, 0.99, 0.01)
  for(k in 1:length(beta_seq_two)){
    design <- optimal.KJ(shape,S0,x0,hr,tf,rate,0.2,beta_seq_two[k],dist="WB")
    design$Two_stage$alpha <- 0.2
    design$Two_stage$beta <- beta_seq_two[k]
    design$Two_stage <- design$Two_stage[, c(9, 10, 1:8)]
    if(design$Two_stage$n == tot_size){
      return(design)
    }
  }
}
