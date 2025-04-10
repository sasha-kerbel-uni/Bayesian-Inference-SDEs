library(MASS)
#create synthetic data
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

set.seed(1)
out<-SVsim(252,0.01,c(5,-3.5),c(0.01,1,1.0,-3.5,0.5))

par(mfrow = c(2,1),mar = c(0, 4, 0, 2), oma = c(3, 1, 2, 0))
plot(ts(out[,1],start=0,deltat=0.01), ylab = expression(X[t]), xaxt = "n", yaxt = "n", ylim = c(0,25))
axis(2, at=seq(0, 20, by=5))
plot(ts(out[,2],start=0,deltat=0.01), ylab = expression(V[t]), xlab = "t")
#thin to get data at integer times
data<-out[1+(0:252)*100,] #length 253

Bridgesim <- function(end,deltat,init,Xend,theta)
{
  m<-end/deltat  
  mat<-matrix(0,ncol=2,nrow=m+1)  
  mat[1,]<-init
  mat[m+1,1]<-Xend
  llikeQ<-0 #loglike under D&G bridge
  llikeP<-0 #loglike under target
  for (i in 1:m)
  {
    #forward simulate v_process
    mat[i+1,2] <- rnorm(1,
                        mat[i,2]+theta[3]*(theta[4]-mat[i,2])*deltat,
                        theta[5]*sqrt(deltat))
  }
  for (i in 1:(m-1))
  { 
    #simulating x_process with dg_bridge
    mu <- mat[i,1]+(Xend-mat[i,1])/(m+1-i) #not multiplying by delta_tau
    var <- (m-i)*theta[2]*mat[i,1]*exp(0.5*mat[i,2])*deltat/(m+1-i) 
    mat[i+1,1] <- abs(rnorm(1,mu,sqrt(var)))
    #reflecting barrier at x=0
    llikeQ<-llikeQ+dnorm(mat[i+1,1],mu,sqrt(var),log=TRUE)
    llikeP<-llikeP+dnorm(mat[i+1,1],mat[i,1]+theta[1]*mat[i,1]*deltat,  
                         theta[2]*mat[i,1]*exp(0.5*mat[i,2])*sqrt(deltat),log=TRUE)
  }
  llikeP<-llikeP+dnorm(mat[m+1,1],mat[m,1]+theta[1]*mat[m,1]*deltat,  
                       theta[2]*mat[m,1]*exp(0.5*mat[m,2])*sqrt(deltat),log=TRUE)
  return(list(mat,llikeP,llikeQ))
}

particle_filter <- function(xobs,delta_t,delta_tau,theta,K,v0){
  R <- length(xobs)-1
  weights <- numeric(length = R) #estimates for transition densities
  sample_weights <- numeric(length = K) #used for resampling v_samples at each step
  v_samples <- rep(v0,K)
  for (t in 1:R){
    for (j in 1:K){
      bridgesim <- Bridgesim(delta_t,delta_tau,c(xobs[t],v_samples[j]),xobs[t+1],theta)
      v_samples[j] <- bridgesim[[1]][nrow(bridgesim[[1]]),2] #Want the final row
      sample_weights[j] <-  exp(bridgesim[[2]]-bridgesim[[3]])
    }
    v_samples <- sample(v_samples,K,replace = TRUE, prob = sample_weights)
    weights[t] <- mean(sample_weights)
  }
  #returns log likelihood estimate (which is log of product of weights)
  return(sum(log(weights)))
}

particle_filter(data[,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),20,-3.5)

lprior <- function(psi){
  return(dunif(psi[1],min=-4.7,max=-4.5,log=TRUE)+ #log(0.01)
           dunif(psi[2],min=-0.1,max=0.1,log=TRUE)+  #log(1)
           dunif(psi[3],min=-0.1,max=0.1,log=TRUE)+ #log(1)
           dunif(psi[4],min=-3.6,max=-3.4,log=TRUE)+ 
           dunif(psi[5],min=-0.8,max=-0.6,log=TRUE)) #log(0.5)
}

pmmh <- function(xobs,N,V,delta_t,K,delta_tau,v0){
  mat <- matrix(0,ncol=6,nrow=N)
  initial_theta <- c(0.01,1,1.0,-3.5,0.5)
  psi <- c(log(initial_theta[1]),log(initial_theta[2]),log(initial_theta[3]),initial_theta[4],log(initial_theta[5])) #initalise psi
  mat[1,] <- c(psi,0)
  nacceptances <- 0
  llikepsi <- -1e12 #initialised so first proposal is always accepted
  for (i in 2:N){
    can <- mvrnorm(n=1,psi,V)
    lpriorcan <- lprior(can)
    lpriorpsi <- lprior(psi)
    llikecan <- particle_filter(xobs,delta_t,delta_tau,
                                c(exp(can[1]),exp(can[2]),exp(can[3]),can[4],exp(can[5])),
                                K,-3.5)
    laprob <- llikecan+lpriorcan-llikepsi-lpriorpsi
    #print(llikecan)
    u <- runif(1)
    if(log(u) < laprob){
      psi <- can #accept 
      llikepsi <- llikecan
      nacceptances <- nacceptances + 1
    }
    mat[i,] <- c(psi,llikecan)
  }
  acceptance_rate <- nacceptances/(N-1)
  print(acceptance_rate)
  return(mat)
}

pf_var_log_estimator <- function(K, trials = 100){
    results <- replicate(trials, particle_filter(data[1:100,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),K,-3.5))
    return(var(results))
}

pf_var_log_estimator(40,100)
#For n=253, K=20 -> 18, K=50 -> 9.8, K =100 -> 5.9, K = 150 -> 4.3
#For n=100, K=100 -> 0.75, K=50 -> 1.1, K=40 -> 1.5


set.seed(1)
pmmhout_pf <- pmmh(data[1:100,1],N=100,V=diag(rep(0.0001,5)),1,K=40,0.1,-3.5)
var(pmmhout_pf[,6]) 
plot(ts(exp(pmmhout_pf)))

set.seed(2)
pmmhout2_pf <- pmmh(data[1:100,1],N=1000,V=var(pmmhout_pf[,1:5]),1,K=40,0.1,-3.5)
var(pmmhout2_pf[,6])
plot(ts(exp(pmmhout2_pf)))

set.seed(3)
pmmhout3_pf <- pmmh(data[1:100,1],N=10000,V=var(pmmhout2_pf[,1:5]),1,K=40,0.1,-3.5)
#Took 2 hours to run
var(pmmhout3_pf[,6])
pmmhout3_burnin_pf <- pmmhout3_pf[1001:10000,]
unscaled_pf <- array(c(exp(pmmhout3_burnin_pf[,c(1:3)]),pmmhout3_burnin_pf[,4],exp(pmmhout3_burnin_pf[,5])),dim = c(9000,5))
colnames(unscaled_pf) <- c("theta1","theta2","theta3","theta4","theta5")

time_taken_pf <- system.time({
  pmmhout3_pf <- pmmh(data[1:100,1],N=10000,V=var(pmmhout2_pf[,1:5]),1,K=40,0.1,-3.5)
})
time_taken_pf
effectiveSize(unscaled_pf)/time_taken_pf[3]

round(effectiveSize(unscaled_pf)/time_taken_25[3],4)


round(apply(unscaled_pf, MARGIN = 2,mean),4)
round(apply(unscaled_pf, MARGIN = 2,sd),4)

round(apply(unscaled_pf,MARGIN =2, quantile, probs = c(0.025,0.975)),4)

output_dir <- "C:\\Users\\kerbe\\OneDrive - Durham University\\Course\\Project IV\\Plots"
pdf(file=file.path(output_dir, "SV_PF_PMMH_%d.pdf"), onefile=F, width = 3, height = 2) 


par(mfrow = c(1,1),mar = c(3, 4, 1, 1), mgp = c(2,1,0))
plot(ts(unscaled_pf[,1]),xlab = "t",ylab = expression(theta[1]),main = "")
plot(density(unscaled_pf[,1], bw = 0.0004), xlab = expression(theta[1]), main = "")
abline(v=0.01, col = "red")
plot(ts(unscaled_pf[,2]),xlab = "t",ylab = expression(theta[2]),main = "")
plot(density(unscaled_pf[,2], bw = 0.03),xlab = expression(theta[2]), main = "")
abline(v=1,col="red")
plot(ts(unscaled_pf[,3]),xlab = "t",ylab = expression(theta[3]),main = "")
plot(density(unscaled_pf[,3], bw = 0.06),xlab = expression(theta[3]), main = "")
abline(v=1,col="red")
plot(ts(unscaled_pf[,4]),xlab = "t",ylab = expression(theta[4]),main = "")
plot(density(unscaled_pf[,4], bw = 0.03),xlab = expression(theta[4]), main = "")
abline(v=-3.5,col="red")
plot(ts(unscaled_pf[,5]),xlab = "t",ylab = expression(theta[5]),main = "")
plot(density(unscaled_pf[,5], bw = 0.02),xlab = expression(theta[5]), main = "")
abline(v=0.5,col="red")

dev.off()
