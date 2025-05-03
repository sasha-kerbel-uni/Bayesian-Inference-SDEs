library("MASS")

GBM_exact <- function(end_time=1, delta_t, x0=1, mu, sigma_sq){
  X <- vector(length = end_time/delta_t)
  X[1] <- x0
  for (i in 2:length(X)){
    X[i] <- rlnorm(1,
                   log(X[i-1])+(mu-sigma_sq/2)*delta_t,
                   sqrt(sigma_sq*delta_t))
  }
  return(X)
}

lprior <- function(psi){
  #psi= (log(mu), log(sigma))
  return(sum(dnorm(psi,0,1,log=TRUE)))
}

llike_exact <- function(x,psi,delta_t){
  #x = (x0,x_deltat,x_2deltat,...,x_ndeltat)
  return(sum(dlnorm(x[2:length(x)],
                    log(x[1:(length(x) - 1)])+(exp(psi[1])-0.5*exp(2*psi[2]))*delta_t,
                    exp(psi[2])*sqrt(delta_t),
                    log = TRUE)))
}

llike_em <- function(x,psi,delta_t){
  return(sum(dnorm(x[2:length(x)],
                   x[1:(length(x) - 1)]+exp(psi[1])*x[1:(length(x) - 1)]*delta_t,
                   exp(psi[2])*x[1:(length(x) - 1)]*sqrt(delta_t),
                   log = TRUE)))
}

mh <- function(exact_likelihood,synthetic_data,delta_t,N,V){
  mat <- matrix(0,ncol=2,nrow=N)
  psi <- c(0,0) #initalise psi
  mat[1,] <- psi
  nacceptances <- 0
  llikepsi<- -1e6
  for (i in 2:N){
    can <- mvrnorm(n=1,psi,V)     #propose psi* (i.e can)
    
    if (exact_likelihood){
      llikecan <- llike_exact(synthetic_data,can,delta_t)
    }
    else{
      llikecan <- llike_em(synthetic_data,can,delta_t)
    }
    
    lpriorcan <- lprior(can)
    lpriorpsi <- lprior(psi)
    
    laprob <- llikecan+lpriorcan-llikepsi-lpriorpsi
    
    u <- runif(1)
    if(log(u) < laprob){
      psi <- can #accept 
      llikepsi<-llikecan
      nacceptances <- nacceptances + 1
    }
    mat[i,] <- psi
  }
  acceptance_rate <- nacceptances/(N-1)
  print(acceptance_rate)
  return(mat)
  #return(acceptance_rate)
}

#simulate data
set.seed(1)
data<-GBM_exact(end_time=10, 0.1, 1, 0.2, 0.1) #sigma=0.316
plot(ts(data,start=0,deltat=0.1), xlab = "t", ylab = expression(X[t]))
axis(2, at = 1, labels = 1, las = 1)

c_dg_bridge <- function(delta_t=1,delta_tau=0.1,xstart=1,xend,sigma_sq,u){
  rho = 0.95
  m <- round(delta_t/delta_tau)
  x <- numeric(length = m+1)
  x[1] <- xstart
  x[m+1] <- xend
  
  for (i in 1:(m-1)){ 
    mu <- x[i]+(xend-x[i])/(m+1-i)
    var <- (m-i)*sigma_sq*x[i]^2*delta_tau/(m+1-i)
    x[i+1] <- abs(mu+sqrt(var)*u[i])
    #reflecting barrier at x=0
  }
  
  return(x)
}

log_em_prod <- function(x,mu,sigma_sq,delta_tau){
  m <- length(x)-1
  return (sum(dnorm(x[2:(m+1)],
                    x[1:m]+mu*x[1:m]*delta_tau,
                    x[1:m]*sqrt(sigma_sq*delta_tau), 
                    log = TRUE)))
}

log_dg_prop <- function(x,sigma_sq,delta_tau){
  m <- length(x)-1
  total = 0
  for (i in 1:(m-1)){
    total <- total + dnorm(x[i+1],
                           x[i]+(x[m+1]-x[i])/(m+1-i),
                           x[i]*sqrt(sigma_sq*delta_tau*(m-i)/(m+1-i)),
                           log = TRUE)
    
  }
  return(total)
}

w_log_scale <- function(bridge,mu,sigma_sq,delta_tau){
  
  em <- log_em_prod(bridge,mu,sigma_sq,delta_tau)
  
  dg_prop <- log_dg_prop(bridge,sigma_sq,delta_tau)
  
  return(em - dg_prop)
}

c_weighted_bridge_resampling <- function(delta_t,K,delta_tau,xstart,xend,mu,sigma_sq,u){
  #K is number of bridges simulated
  m <- round(delta_t/delta_tau)
  proposals <- matrix(NA,nrow = m+1,ncol = K)
  weights <- numeric(length = K)
  for (i in 1:K){
    proposals[,i] <- c_dg_bridge(delta_t,delta_tau,xstart,xend,sigma_sq,u[i,])
    weights[i] <- exp(w_log_scale(proposals[,i],mu,sigma_sq,delta_tau))
  }
  resampled_indices <- sample(1:K,size = K, replace = TRUE, prob = weights)
  resampled <- proposals[,resampled_indices]
  avg_unnormalised_weight <- mean(weights)
  return(list(proposals,resampled,avg_unnormalised_weight))
}

c_estimated_log_likelihood <- function(x,mu,sigma_sq,delta_t,K,delta_tau,u)
{
  n <- length(x) #n is the number of observations 
  density_estimate <- 0
  for (i in 1:(n-1)){
    out<-c_weighted_bridge_resampling(delta_t,K,delta_tau,x[i],x[i+1],mu,sigma_sq,u[i,,])
    density_estimate <- density_estimate +  log(out[[3]])
  }
  return(density_estimate)
}

cpmmh <- function(synthetic_data,N,V,delta_t,K,delta_tau,rho){
  mat <- matrix(0,ncol=3,nrow=N)
  psi <- c(log(0.2),log(sqrt(0.1))) #initalise psi
  mat[1,] <- c(psi,0)
  m <- delta_t/delta_tau
  n <- length(synthetic_data)
  u <- array(rnorm(n*K*(m-1),0,1), dim = c(n,K,(m-1)))
  nacceptances <- 0
  llikepsi <- -1e12 #initialised so first proposal is always accepted
  for (i in 2:N){
    can <- mvrnorm(n=1,psi,V)
    ucan <- array(rnorm(n*K*(m-1),rho*u,sqrt(1-rho^2)),dim = c(n,K,(m-1))) #CHECK THIS IS UPDATING UCAN ARRAY CORRECTLY
    lpriorcan <- lprior(can)
    lpriorpsi <- lprior(psi)
    llikecan <-  
      c_estimated_log_likelihood(synthetic_data,exp(can[1]),exp(2*can[2]),delta_t,K,delta_tau,ucan) ###Interpreting psi = c(log(mu),log(sigma)) 
    laprob <- llikecan+lpriorcan-llikepsi-lpriorpsi
    w <- runif(1)
    if(log(w) < laprob){
      psi <- can #accept 
      u <- ucan
      llikepsi <- llikecan
      nacceptances <- nacceptances + 1
    }
    mat[i,] <- c(psi,llikecan)
  }
  acceptance_rate <- nacceptances/(N-1)
  print(acceptance_rate)
  return(mat)
}

c_var_log_est <- function(K,trials = 100){
  results <- replicate(trials, c_estimated_log_likelihood(data2,0.2,0.1,1,K,0.1,array(rnorm(length(data2)*K*(10-1),0,1), dim = c(length(data2),K,(10-1)))))
  return(var(results))
}

c_var_log_est(5,1000)

set.seed(1)
outcpmmh<-cpmmh(data2,1000,diag(rep(0.05,2)),1,5,0.1,0.95)
plot(ts(exp(outcpmmh[,1:2])))
var(outcpmmh[,3])

#tune
Vtune<-var(outcpmmh[,1:2])
#Re-run
set.seed(1)
outcpmmh2<-cpmmh(data2,5000,Vtune,1,50,0.1,0.95)

time_taken_cpmmh <- system.time({
  outcpmmh2<-cpmmh(data2,5000,Vtune,1,50,0.1,0.95)
})

effectiveSize(outcpmmh2)/time_taken_cpmmh[3]

output_dir <- "C:\\Users//HP\\OneDrive - Durham University\\Course\\Project IV\\Plots"
pdf(file=file.path(output_dir, "GBM PMMH posteriors.pdf"), onefile=F, width = 5, height = 5) 

plot(ts(exp(outcpmmh2)))
#Compare output across schemes
par(mfrow=c(2,1),mar = c(3, 3, 2, 1), oma = c(0,0,0,0), mgp = c(1.5, 0.5, 0))
plot(density(exp(out2[50:10000,1])), xlab = expression(mu),main ="", ylim =c(0,10), lwd =2)
lines(density(exp(outEM2[50:10000,1])),col=2)
lines(density(exp(outpmmh2[,1])),col=3)
lines(density(exp(outcpmmh2[,1])),col=4)
abline(v=0.2)
plot(density(exp(out2[50:10000,2])), xlab = expression(sigma),main = "", ylim = c(0,15), lwd =2)
lines(density(exp(outEM2[50:10000,2])),col=2)
lines(density(exp(outpmmh2[,2])),col=3)
lines(density(exp(outcpmmh2[,2])),col=4)
abline(v=sqrt(0.1))

dev.off()

