eu_payoff <- function(ST,K){
  return (max(ST-K,0))
}

STsim <- function(end,deltat,init,theta)
{
  n<-end/deltat+1   
  mat<-matrix(0,ncol=2,nrow=n)  
  mat[1,]<-init
  for (i in 2:n)
  {
    #v_process
    mat[i,2] <- rnorm(1,
                      mat[i-1,2]+theta[3]*(theta[4]-mat[i-1,2])*deltat,
                      theta[5]*sqrt(deltat))
    #x_process
    mat[i,1] <- rnorm(1,
                      mat[i-1,1]+theta[1]*mat[i-1,1]*deltat,
                      theta[2]*mat[i-1,1]*exp(0.5*mat[i-1,2])*sqrt(deltat))
  }
  return(mat[n,1])
}

#t is expiry time
MC_price <- function(M,r,t,theta,S0,K){
  price <- 0
  for (i in 1:M){
    #using theta1=mu-sigma^2
    ST <- STsim(t,0.01,c(S0,theta[4]),c(theta[1]-0.5*exp(theta[4]),theta[2:5])) #think about impact of discretisation can I get away with more than 0.01
    price <- price + eu_payoff(ST,K)
  }
  price <- (1/M)*exp(-r*t)*price
  return(price)
}


bsm_formula_eu <- function(r,t,sigma,S0,K){
  d1 <- (log(S0/K)+(r+0.5*sigma^2)*t)/(sigma*sqrt(t))
  d2 <- d1-sigma*sqrt(t)
  price <- S0*pnorm(d1)-exp(-r*t)*K*pnorm(d2)
  return(price)
}


set.seed(1)
round(bsm_formula_eu(0.01-0.5*exp(-3.5),100,exp(0.5*-3.5),1,1.05),4)#analytic formula
MC_price(1000,0.01-0.5*exp(-3.5),100,c(0.01,1,0,-3.5,0),1,1.05) #MC with BSM
MC_price(1000,0.01-0.5*exp(-3.5),100,c(0.01,1,1.0,-3.5,0.5),1,1.05) #MC with SV


multi_theta_MC_price <- function(M,r,t,theta_samples,S0,K){
  prices <- numeric(length = nrow(theta_samples))
  for (i in 1:nrow(theta_samples)){
    prices[i] <- MC_price(M,r,t,theta_samples[i,],S0,K)
  }
  return(prices)
}

thinned_theta_samples <- unscaled_pf[90*1:100,]

set.seed(1) #didn't seed on previous run
prices <- multi_theta_MC_price(1000,0.01-0.5*exp(-3.5),100,thinned_theta_samples,1,1.05)

hist(prices, breaks = nbins(prices, 0.08), main = "", labels = FALSE)
axis(1, at = seq(0, 1.2, by = 0.1))
summary(prices)
round(quantile(prices, prob = c(0.025,0.975)),4)

nbins <- function(x, bin_width){
  return(seq(min(x) - bin_width,
             max(x) + bin_width,
             by = bin_width))
}

SVsim <- function(end,deltat,init,theta)
{
  n<-end/deltat+1   
  mat<-matrix(0,ncol=2,nrow=n)  
  mat[1,]<-init
  for (i in 2:n)
  {
    #v_process
    mat[i,2] <- rnorm(1,
                      mat[i-1,2]+theta[3]*(theta[4]-mat[i-1,2])*deltat,
                      theta[5]*sqrt(deltat))
    #x_process
    mat[i,1] <- rnorm(1,
                      mat[i-1,1]+theta[1]*mat[i-1,1]*deltat,
                      theta[2]*mat[i-1,1]*exp(0.5*mat[i-1,2])*sqrt(deltat))
  }
  return(mat)
}

MC_price_asian <- function(M,r,t,theta,S0,K){
  price <- 0
  for (i in 1:M){
    #using theta1=mu-sigma^2
    X <- SVsim(t,0.01,c(S0,theta[4]),c(theta[1]-0.5*exp(theta[4]),theta[2:5]))[,1] #think about impact of discretisation can I get away with more than 0.01
    price <- price + eu_payoff(mean(X),K) #eu_payoff but with mean(X) rather than XT
  }
  price <- (1/M)*exp(-r*t)*price
  return(price)
}

MC_price_asian(1000,0.01-0.5*exp(-3.5),100,c(0.01,1,1.0,-3.5,0.5),1,1.05) 


multi_theta_MC_price_asian <- function(M,r,t,theta_samples,S0,K){
  prices <- numeric(length = nrow(theta_samples))
  for (i in 1:nrow(theta_samples)){
    prices[i] <- MC_price_asian(M,r,t,theta_samples[i,],S0,K)
  }
  return(prices)
}

set.seed(1)
asian_prices <- multi_theta_MC_price_asian(1000,0.01-0.5*exp(-3.5),100,thinned_theta_samples,1,1.05)

hist(asian_prices, breaks = nbins(asian_prices, 0.03), main = "", labels = FALSE, xlab = "prices")
axis(1, at = seq(0, 0.6, by = 0.05))
summary(asian_prices)
quantile(asian_prices, prob = c(0.025,0.975))
