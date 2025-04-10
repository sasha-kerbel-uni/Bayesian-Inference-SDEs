library(MASS)
library(coda)

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

par(mfrow = c(2,1),mar = c(0, 3, 0, 2), oma = c(3, 1, 2, 0))
plot(ts(out[,1],start=0,deltat=0.01), ylab = expression(X[t]), xaxt = "n", yaxt = "n", ylim = c(0,25))
axis(2, at=seq(0, 20, by=5))
plot(ts(out[,2],start=0,deltat=0.01), ylab = expression(V[t]), xlab = "t")
mtext("t", side = 1, line = 2, outer = TRUE)
#thin to get data at integer times
data<-out[1+(0:252)*100,] #length 253



par(mfrow=c(3,1))
plot(ts(data[,1],start=0,deltat=1), col = "red", type = "p", pch =20, xlab = "t", ylab = expression(X[t]))

#par(mfrow=c(2,1))
plot(ts(out[,1],start=0,deltat=0.01),xlab = "t", ylab = expression(X[t]), type = "p")
points(ts(data[,1],start=0,deltat=1), col = "red", pch = 20)
plot(ts(exp(out[,2]),start=0,deltat=0.01), xlab = "t", ylab = expression(exp(V[t])), type = "p")

#forward simulate V process then bridge X process
#Xend is obs at end point (could be inter-obs time)
#calculates likelihood under bridge and target (Euler-Maruyama)
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

par(mfrow = c(1,1))

test<-Bridgesim(1,0.2,c(5,-3.5),4.989,c(0.01,1,1.0,-3.5,0.5))
plot(ts(test[[1]][,1],start=0,deltat=0.2),ylim=c(4,7))
for(i in 1:100){
test<-Bridgesim(1,0.2,c(5,-3.5),4.989,c(0.01,1,1.0,-3.5,0.5))
lines(ts(test[[1]][,1],start=0,deltat=0.2),col=i)  
}

plot(ts(test[[1]][,2],start=0,deltat=0.2),ylim=c(-4.5,-2.5))
for(i in 1:100){
  test<-Bridgesim(1,0.2,c(5,-3.5),4.989,c(0.01,1,1.0,-3.5,0.5))
  lines(ts(test[[1]][,2],start=0,deltat=0.2),col=i)  
}

#p^_e(x^0|theta)
#end is the end time that the bridge should span
log_est_likelihood <- function(xobs,end,deltat,theta,K,v0)
{
  n=length(xobs)-1
  m=end/deltat
  wts <- rep(0,K)
  for (k in 1:K) #loop over importance samples
  {
   y0<-c(xobs[1],v0)  
   for(i in 1:n) #loop over obs intervals
   {
    out<-Bridgesim(end,deltat,y0,xobs[i+1],theta)
    y0<-c(xobs[i+1],out[[1]][m+1,2]) #updating start for next bridge
    wts[k]<-wts[k]+out[[2]]-out[[3]]
   }
   wts[k] <- exp(wts[k]) #changing weights back to original scale (after being log-scaled)
  }
  return(log(mean(wts)))
}

log_est_likelihood(data[,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),20,-3.5)

simulate_trajectory <- function(xobs,end,deltat,theta,v0){
  n=length(xobs)-1
  m=1/deltat
  mat <- matrix(NA,nrow = n*m+1,ncol = 2)
  y0<-c(xobs[1],v0)  
  mat[1,] <- y0
  for(i in 1:n) #loop over obs intervals
  {
    out<-Bridgesim(end,deltat,y0,xobs[i+1],theta)
    mat[((i-1)*m+2):((i*m+1)),] <- out[[1]][2:(m+1),]
    y0<-c(xobs[i+1],out[[1]][m+1,2]) #updating start for next bridge
  }
  return(mat)
}

#Plotting "true" trajectories against simulations
set.seed(1)
test_traj <- simulate_trajectory(data[,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),-3.5)
par(mfrow = c(2,1),mar = c(0, 3, 0, 2), oma = c(3, 1, 2, 0), mgp = c(2, 0.5, 0))
plot(ts(data[,2], deltat=1), lwd=2, ylab = expression(V[t]), xaxt = "n", xlab = "", las =1)
lines(ts(test_traj[,2], deltat=0.2), col = "red")

plot(ts(data[,1], deltat=1), lwd =2, xlab = "t", ylab = expression(X[t]), las =1)
#plot(ts(out,deltat=0.01))
lines(ts(test_traj[,1], deltat=0.2), col = "red")
mtext("t", side = 1, line = 1)  # Adjust 'line' value to control distance





par(mfrow = c(2,1),mar = c(0, 4, 0, 2), oma = c(3, 1, 2, 0), mgp = c(2, 0.5, 0))
plot(ts(test_traj[,2], deltat=0.2, start =0), col = rgb(91, 45, 135,maxColorValue = 255),xlab = "", ylab = expression(V[t]),xaxt= "n", las =1)
plot(ts(test_traj[,1], deltat=0.2, start = 0), col = "red", xlab = "t", ylab = expression(X[t]), las =1)
points(ts(data[,1],start=0,deltat=1), col = "black", pch =20)
mtext("t", side = 1, line = 1)



#How many bridges?
vec<-rep(0,200)
for(i in 1:200){
  vec[i]<-log_est_likelihood(data[1:25,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),10,-3.5) 
}
var(vec)

required_bridges_sequential <- function(n, Kstart=0){
  variance <- 1e9 #big number (just needs to be bigger than while condition)
  K <- Kstart
  while(variance>1.5){
    K <- K+1
    vec<-rep(0,100)
    for(i in 1:100){
      vec[i]<-log_est_likelihood(data[1:n,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),K,-3.5) 
    }
    variance <- var(vec)
  }
  return(K)
}

required_bridges_sequential(10)

compute_variance <- function(K, n, trials=100) {
  results <- replicate(trials, log_est_likelihood(data[1:n,1],1,0.2,c(0.01,1,1.0,-3.5,0.5),K,-3.5))
  return(var(results))
}

compute_variance(50,100)
compute_variance(500,100)

required_bridges_binary_search <- function(n, K_min=1, K_max=500, tolerance = 0.5){
  while (K_min < K_max){
    K_mid <- floor((K_min + K_max) / 2)  # Midpoint for binary search
    var_mid <- compute_variance(K_mid, n)  # Compute variance at K_mid
    
    if(var_mid<1.5){
      K_max <- K_mid
    }
    else{
      K_min <- K_mid+1
    }
    
    if (abs(K_max-K_min)<tolerance){ #choose tolerance <1 to find the "exact" K - although will vary as it is stochastic
      break
    }
  }
  return(K_min)
}

required_bridges_binary_search(50,1,60)
compute_variance(53,50)
compute_variance(500,60)

set.seed(1)
x <- numeric(6)
x[1] <- required_bridges_sequential(10,1)
x[2] <- required_bridges_sequential(20,1)
x[3] <- required_bridges_sequential(30,1)
x[4] <- required_bridges_sequential(40,9) 
x[5] <- required_bridges_sequential(50,38)
x[6] <- required_bridges_sequential(60,200)
print(compute_variance(200,n=60,100))
#not technically increasing as there is stochasicity
x=c(2,5,9,24,45,201)
grid <- c(10,20,30,40,50,60)

par(mfrow = c(1,1),mar = c(4, 5, 2, 2), mgp = c(3,1,0))
plot(grid,x,type = "l",xlab = "n", ylab = "K")

set.seed(1)
K <- seq(2,40,by=4)
Kvariance <- numeric(length = length(K))
for (i in 1:length(K)){
  k <- K[i]
  Kvariance[i] <- compute_variance(k,n=50,trials = 200)
}


par(mfrow = c(1,1),mar = c(4, 5, 2, 2), mgp = c(3,1,0))
plot(K,Kvariance, xlab = "K", type = "l",ylab = expression(Var(log(hat(f) * (x ~ "|" ~ theta)))))
Kvariance

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
    llikecan <- log_est_likelihood(xobs,delta_t,delta_tau,
                                   c(exp(can[1]),exp(can[2]),exp(can[3]),can[4],exp(can[5])),
                                   K,v0)
    laprob <- llikecan+lpriorcan-llikepsi-lpriorpsi
    print(llikecan)
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

set.seed(1)
pmmhout_25 <- pmmh(data[1:25,1],N=100,V=diag(rep(0.0001,5)),1,K=10,0.1,-3.5)
var(pmmhout_25[,6]) #17
plot(ts(exp(pmmhout_25)))

set.seed(2)
pmmhout2_25 <- pmmh(data[1:25,1],N=1000,V=var(pmmhout_25[,1:5]),1,K=10,0.1,-3.5)
var(pmmhout2_25[,6])
plot(ts(exp(pmmhout2_25)))

set.seed(3)
pmmhout3_25 <- pmmh(data[1:25,1],N=10000,V=var(pmmhout2_25[,1:5]),1,K=10,0.1,-3.5)
var(pmmhout3_25[,6])
pmmhout3_burnin_25 <- pmmhout3_25[1001:10000,]
unscaled_25 <- array(c(exp(pmmhout3_burnin_25[,c(1:3)]),pmmhout3_burnin_25[,4],exp(pmmhout3_burnin_25[,5])),dim = c(9000,5))
colnames(unscaled_25) <- c("theta1","theta2","theta3","theta4","theta5")

time_taken_25 <- system.time({
  pmmhout3_25 <- pmmh(data[1:25,1],N=10000,V=var(pmmhout2_25[,1:5]),1,K=10,0.1,-3.5)
})
time_taken_25
round(effectiveSize(unscaled_25)/time_taken_25[3],4)

str(unscaled_25)
round(apply(unscaled_25, MARGIN = 2,sd),4)

round(apply(unscaled_25,MARGIN =2, quantile, probs = c(0.025,0.975)),4)
#final plots
output_dir <- "C:\\Users\\kerbe\\OneDrive - Durham University\\Course\\Project IV\\Plots"
pdf(file=file.path(output_dir, "SV_PMMH_25_%d.pdf"), onefile=F, width = 3, height = 2) 


par(mfrow = c(1,1),mar = c(3, 4, 1, 1), mgp = c(2,1,0))
plot(ts(unscaled_25[,1]),xlab = "t",ylab = expression(theta[1]),main = "")
plot(density(unscaled_25[,1], bw = 0.0004), xlab = expression(theta[1]), main = "")
abline(v=0.01, col = "red")
plot(ts(unscaled_25[,2]),xlab = "t",ylab = expression(theta[2]),main = "")
plot(density(unscaled_25[,2], bw = 0.03),xlab = expression(theta[2]), main = "")
abline(v=1,col="red")
plot(ts(unscaled_25[,3]),xlab = "t",ylab = expression(theta[3]),main = "")
plot(density(unscaled_25[,3], bw = 0.06),xlab = expression(theta[3]), main = "")
abline(v=1,col="red")
plot(ts(unscaled_25[,4]),xlab = "t",ylab = expression(theta[4]),main = "")
plot(density(unscaled_25[,4], bw = 0.03),xlab = expression(theta[4]), main = "")
abline(v=-3.5,col="red")
plot(ts(unscaled_25[,5]),xlab = "t",ylab = expression(theta[5]),main = "")
plot(density(unscaled_25[,5], bw = 0.02),xlab = expression(theta[5]), main = "")
abline(v=0.5,col="red")

dev.off()



#FOR N=100 BAD EXAMPLE
set.seed(1)
pmmhout_100 <- pmmh(data[1:100,1],N=100,V=diag(rep(0.0001,5)),1,K=10,0.1,-3.5)
var(pmmhout_100[,6]) #17
plot(ts(exp(pmmhout_100)))

set.seed(2)
pmmhout2_100 <- pmmh(data[1:100,1],N=1000,V=var(pmmhout_100[,1:5]),1,K=10,0.1,-3.5)
var(pmmhout2_100[,6])
plot(ts(exp(pmmhout2_100)))

set.seed(3)
pmmhout3_100 <- pmmh(data[1:100,1],N=10000,V=var(pmmhout2_100[,1:5]),1,K=10,0.1,-3.5)
var(pmmhout3_100[,6])
pmmhout3_burnin_100 <- pmmhout3_100[1001:10000,]
unscaled_100 <- array(c(exp(pmmhout3_burnin_100[,c(1:3)]),pmmhout3_burnin_100[,4],exp(pmmhout3_burnin_100[,5])),dim = c(9000,5))
colnames(unscaled_100) <- c("theta1","theta2","theta3","theta4","theta5")

time_taken_100 <- system.time({
  pmmhout3_100 <- pmmh(data[1:100,1],N=10000,V=var(pmmhout2_100[,1:5]),1,K=10,0.1,-3.5)
})

acc_prob_100_manual <- 0.03040304

#final plots
output_dir <- "C:\\Users\\kerbe\\OneDrive - Durham University\\Course\\Project IV\\Plots"
pdf(file=file.path(output_dir, "SV_PMMH_100_BAD_%d.pdf"), onefile=F, width = 3, height = 2) 

par(mfrow = c(1,1),mar = c(3, 4, 1, 1), mgp = c(2,1,0))
plot(ts(unscaled_100[,1]),xlab = "t",ylab = expression(theta[1]),main = "")
plot(ts(unscaled_100[,2]),xlab = "t",ylab = expression(theta[2]),main = "")
plot(ts(unscaled_100[,3]),xlab = "t",ylab = expression(theta[3]),main = "")
plot(ts(unscaled_100[,4]),xlab = "t",ylab = expression(theta[4]),main = "")
plot(ts(unscaled_100[,5]),xlab = "t",ylab = expression(theta[5]),main = "")

dev.off()



