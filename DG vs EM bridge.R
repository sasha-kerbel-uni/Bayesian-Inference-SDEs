dg_bridge <- function(delta_t=1,delta_tau=0.1,xstart=1,xend,sigma_sq){
  m <- round(delta_t/delta_tau)
  x <- numeric(length = m+1)
  x[1] <- xstart
  x[m+1] <- xend
  
  for (i in 1:(m-1)){ 
    mu <- x[i]+(xend-x[i])/(m+1-i)
    var <- (m-i)*sigma_sq*x[i]^2*delta_tau/(m+1-i)
    x[i+1] <- abs(rnorm(1,mu,sqrt(var)))
    #reflecting barrier at x=0
  }
  
  return(x)
}

vec=dg_bridge(1,0.01,1,5,0.1)
plot(ts(vec,start=0,deltat=0.01))

em_bridge <- function(delta_t=1, delta_tau=0.1, xstart=1,xend, mu, sigma_sq){
  X <- numeric(length = (delta_t/delta_tau)+1)
  X[1] <- xstart
  for (i in 2:(length(X)-1)){
    X[i] <- abs(rnorm(1,
                      X[i-1]+mu*X[i-1]*delta_tau,
                      X[i-1]*sqrt(sigma_sq*delta_tau)))
    #Adding absolute value is effectively a reflecting barrier at 0
  }
  X[delta_t/delta_tau+1] <- xend
  return(X)
}


#Plots
plot_bridges <- function(delta_t,delta_tau,xstart,xend,mu,sigma_sq){
  par(mfrow = c(1,2))
  set.seed(1)
  plot(ts(dg_bridge(delta_t,delta_tau ,xstart, xend, sigma_sq)),ylim = c(0,6), main = "MDB", ylab = expression(X[t]), xaxt = "n")
  axis(1, at = 1:6, labels = seq(0,delta_t, by = delta_tau))
  for (i in 2:50){
    set.seed(i)
    lines(ts(dg_bridge(delta_t,delta_tau,xstart, xend, sigma_sq)))
  }
  
  set.seed(1)
  
  plot(ts(em_bridge(delta_t,delta_tau,xstart, xend, mu, sigma_sq)),ylim = c(0,6), main = "Myopic EM", ylab = expression(X[t]), xaxt = "n")
  axis(1, at = 1:6, labels = seq(0,delta_t, by = delta_tau))
  for (i in 2:50){
    set.seed(i)
    lines(ts(em_bridge(delta_t,delta_tau,xstart, xend, mu, sigma_sq)))
  }
}

plot_bridges(delta_t=1,delta_tau=0.2,xstart=1,xend=5,mu=0.1,sigma_sq=0.2)

