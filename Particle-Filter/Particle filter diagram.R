library(MASS)

delta_tau <- 0.2

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

first_bridge_sim <- function(K){
  trajectories <- array(NA, dim = c(1/delta_tau+1,2,K))
  weights <- numeric(length = K)
  for (k in 1:K){
    bridgesim <- Bridgesim(1, #delta_t
                           delta_tau, #delta_tau
                           c(data[1,1],data[1,2]), #c(xobs[t],v_samples[j])
                           data[2,1], #xobs[t+1]
                           c(0.01,1,1.0,-3.5,0.5)) #theta
    trajectories[,,k] <- bridgesim[[1]]
    weights[k] <- exp(bridgesim[[2]]-bridgesim[[3]])
    }
  return(list(trajectories,weights))
}

set.seed(3)
K=10
first_bridge <- first_bridge_sim(K)
#data generated according to SVsim

#plotting x process
par(mfrow = c(1,2),mar = c(4, 5, 2, 2), mgp = c(3,1,0))
plot(ts(first_bridge[[1]][,1,1],start=0,deltat=delta_tau),ylim=c(4,7), xlim = c(0,2), xlab ="t", ylab = expression(X[t]))
points(0,data[1,1], pch = 4)
points(1,data[2,1], pch = 4)
points(2,data[3,1], pch = 4)
for (i in 1:9){
  lines(ts(first_bridge[[1]][,1,i],start=0,deltat=delta_tau), col = i+1)
}
#mtext("t", side = 1, line = 2, outer = FALSE)

#plotting v process
plot(ts(first_bridge[[1]][,2,1],start=0,deltat=delta_tau), ylim = c(-4.5,-2), xlim = c(0,2), xlab = "t", ylab = expression(V[t]))
points(0,data[1,2])
for (i in 1:9){
  lines(ts(first_bridge[[1]][,2,i],start=0,deltat=delta_tau), col = i+1)
  points(1,first_bridge[[1]][1/delta_tau+1,2,i], cex = 3*first_bridge[[2]][i]/sqrt(sum(first_bridge[[2]]^2)))
}

second_bridge_sim <- function(first_bridge,K){
  v_samples <- sample(first_bridge[[1]][nrow(first_bridge[[1]]),2,],K,replace = TRUE, prob = first_bridge[[2]])
  v_samples_indices <- match(v_samples,first_bridge[[1]][nrow(first_bridge[[1]]),2,])
  trajectories <- array(NA, dim = c(1/delta_tau+1,2,K))
  weights <- numeric(length = K)
  for (k in 1:K){
    bridgesim <- Bridgesim(1, #delta_t
                           delta_tau, #delta_tau
                           c(data[2,1],v_samples[k]), #c(xobs[t],v_samples[j])
                           data[3,1], #xobs[t+1]
                           c(0.01,1,1.0,-3.5,0.5)) #theta
    trajectories[,,k] <- bridgesim[[1]]
    weights[k] <- exp(bridgesim[[2]]-bridgesim[[3]])
  }
  return(list(trajectories,weights, v_samples_indices))
}

set.seed(3)
second_bridge <- second_bridge_sim(first_bridge,100)

#plotting x process
plot(ts(first_bridge[[1]][,1,1],start=0,deltat=delta_tau),ylim=c(4,7), xlim = c(0,2),xlab ="t", ylab = expression(X[t]))
lines(ts(second_bridge[[1]][,1,1],start=1,deltat=delta_tau))
points(0,data[1,1], pch = 4)
points(1,data[2,1], pch = 4)
points(2,data[3,1], pch = 4)
for (i in 1:(K-1)){
  lines(ts(first_bridge[[1]][,1,i],start=0,deltat=delta_tau), col = i+1)
  lines(ts(second_bridge[[1]][,1,i],start=1,deltat=delta_tau), col = second_bridge[[3]][i]+1)
  #need to make color of 2nd bridge match the one which is is continuing
}

#plotting v process
plot(ts(first_bridge[[1]][,2,1],start=0,deltat=delta_tau), ylim = c(-4.5,-2), xlim = c(0,2), xlab ="t", ylab = expression(V[t]))
lines(ts(second_bridge[[1]][,2,1],start=1,deltat=delta_tau))
points(0,data[1,2])
for (i in 1:(K-1)){
  lines(ts(first_bridge[[1]][,2,i],start=0,deltat=delta_tau), col = i+1)
  points(1,first_bridge[[1]][1/delta_tau+1,2,i], cex = 3*first_bridge[[2]][i]/sqrt(sum(first_bridge[[2]]^2)))
  lines(ts(second_bridge[[1]][,2,i],start=1,deltat=delta_tau), col = second_bridge[[3]][i]+1)
  points(2,second_bridge[[1]][1/delta_tau+1,2,i], cex = 10*second_bridge[[2]][i]/sqrt(sum(second_bridge[[2]]^2)))
}

#final plots
output_dir <- "C:\\Users\\qnkp75\\OneDrive - Durham University\\Course\\Project IV\\Plots"

#png(file=file.path(output_dir, "particle filter diagram_%d.png"), width = 2000, height = 2000) 
pdf(file=file.path(output_dir, "particle filter diagram_%d.pdf"), onefile=F, width = 5, height = 5) 

par(mfrow = c(1,1),mar = c(3, 4, 1, 1), mgp = c(2,1,0))

plot(ts(first_bridge[[1]][,1,1],start=0,deltat=delta_tau),ylim=c(4,6.5), xlim = c(0,2), xlab ="t", ylab = expression(X[t]))
points(0,data[1,1], pch = 4)
points(1,data[2,1], pch = 4)
points(2,data[3,1], pch = 4)
for (i in 1:9){
  lines(ts(first_bridge[[1]][,1,i],start=0,deltat=delta_tau), col = i+1)
}

plot(ts(first_bridge[[1]][,2,1],start=0,deltat=delta_tau), ylim = c(-4.2,-2.9), xlim = c(0,2), xlab = "t", ylab = expression(V[t]))
points(0,data[1,2])
for (i in 1:9){
  lines(ts(first_bridge[[1]][,2,i],start=0,deltat=delta_tau), col = i+1)
  points(1,first_bridge[[1]][1/delta_tau+1,2,i], pch = 1,cex = 5*first_bridge[[2]][i]/sqrt(sum(first_bridge[[2]]^2)))
}

plot(ts(first_bridge[[1]][,1,1],start=0,deltat=delta_tau),ylim=c(4,6.5), xlim = c(0,2),xlab ="t", ylab = expression(X[t]))
lines(ts(second_bridge[[1]][,1,1],start=1,deltat=delta_tau))
points(0,data[1,1], pch = 4)
points(1,data[2,1], pch = 4)
points(2,data[3,1], pch = 4)
for (i in 1:(K-1)){
  lines(ts(first_bridge[[1]][,1,i],start=0,deltat=delta_tau), col = i+1)
  lines(ts(second_bridge[[1]][,1,i],start=1,deltat=delta_tau), col = second_bridge[[3]][i]+1)
  #need to make color of 2nd bridge match the one which is is continuing
}

plot(ts(first_bridge[[1]][,2,1],start=0,deltat=delta_tau), ylim = c(-4.2,-2.9), xlim = c(0,2), xlab ="t", ylab = expression(V[t]))
lines(ts(second_bridge[[1]][,2,1],start=1,deltat=delta_tau))
points(0,data[1,2])
for (i in 1:(K-1)){
  lines(ts(first_bridge[[1]][,2,i],start=0,deltat=delta_tau), col = i+1)
  points(1,first_bridge[[1]][1/delta_tau+1,2,i], pch = 1, cex = 5*first_bridge[[2]][i]/sqrt(sum(first_bridge[[2]]^2)))
  lines(ts(second_bridge[[1]][,2,i],start=1,deltat=delta_tau), col = second_bridge[[3]][i]+1)
  points(2,second_bridge[[1]][1/delta_tau+1,2,i], pch = 1, cex = 15*second_bridge[[2]][i]/sqrt(sum(second_bridge[[2]]^2)))
}

dev.off()
