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
#seed 5 looks better
set.seed(5)
data<-GBM_exact(end_time=50, 0.1, 1, 0.2, 0.1) #sigma=0.316
plot(ts(data,start=0,deltat=0.1), xlab = "t", ylab = expression(X[t]))
axis(2, at = 1, labels = 1, las = 1)

set.seed(1)
out_pilot<-mh(TRUE,data,0.1,1000,diag(rep(0.04,2)))
#out <- mh(TRUE,data,0.1,10000,2.38^2/2*var(out_pilot))
out <- mh(TRUE,data,0.1,10000,var(out_pilot))
plot(`colnames<-`(ts(exp(out[50:10000,])), c("\u03BC","\u03C3")), main ="")
outEM_pilot<-mh(FALSE,data,0.1,1000,diag(rep(0.04,2)))
#outEM<-mh(FALSE,data,0.1,10000,2.38^2/2*var(outEM_pilot))
outEM<-mh(FALSE,data,0.1,10000,var(outEM_pilot))

plot(`colnames<-`(ts(exp(outEM[50:10000,])), c("\u03BC","\u03C3")), main ="")
par(mfrow=c(2,1))
plot(density(exp(out[50:10000,1])), ylab = expression(mu), main = "",ylim = c(0,5))
lines(density(exp(outEM[50:10000,1])),col=2)
abline(v=0.2)
plot(density(exp(out[50:10000,2])), ylab = expression(sigma), main = "")
lines(density(exp(outEM[50:10000,2])),col=2)
abline(v=sqrt(0.1))


#Try deltat=1
# set.seed(1)
# data2<-GBM_exact(end_time=10, 1, 1, 0.2, 0.1) #sigma=0.316
# plot(ts(data2,start=0,deltat=1))

data2 <- data[10*(0:50)] #thinning the original data rather than generating a new dataset
plot(ts(data2,start=0,deltat=1))

set.seed(1)
out2_pilot<-mh(TRUE,data2,1,1000,diag(rep(0.04,2)))
#out2<-mh(TRUE,data2,1,10000,2.38^2/2*var(out2_pilot))
out2<-mh(TRUE,data2,1,10000,var(out2_pilot))
plot(`colnames<-`(ts(exp(out2[50:10000,])), c("\u03BC","\u03C3")), main ="")
outEM2_pilot<-mh(FALSE,data2,1,1000,diag(rep(0.04,2)))
#outEM2<-mh(FALSE,data2,1,10000,2.38^2/2*var(outEM2_pilot))
outEM2<-mh(FALSE,data2,1,10000,var(outEM2_pilot))
plot(`colnames<-`(ts(exp(outEM2[50:10000,])), c("\u03BC","\u03C3")), main ="")
par(mfrow=c(2,1))
plot(density(exp(out2[50:10000,1])),ylab = expression(mu), main = "")
lines(density(exp(outEM2[50:10000,1])),col=2)
plot(density(exp(out2[50:10000,2])), ylab = expression(sigma), main = "")
lines(density(exp(outEM2[50:10000,2])),col=2)

output_dir <- "C:\\Users\\HP\\OneDrive - Durham University\\Course\\Project IV\\Plots"
pdf(file=file.path(output_dir, "RWM posteriors for varying discretisations.pdf"), onefile=F, width = 5, height = 5) 

par(mfrow=c(2,2),mar = c(3, 3, 2, 1), oma = c(0,0,0,0), mgp = c(1.5, 0.5, 0))
plot(density(exp(out[1000:10000,1])), xlab = expression(mu), main = expression(Delta * t == 0.1), ylab = "")
lines(density(exp(outEM[1000:10000,1])),col=2)
abline(v=0.2)
plot(density(exp(out[1000:10000,2])), xlab = expression(sigma), main = expression(Delta * t == 0.1), ylab ="")
lines(density(exp(outEM[1000:10000,2])),col=2)
abline(v=sqrt(0.1))
plot(density(exp(out2[1000:10000,1])),xlab = expression(mu), main = expression(Delta * t == 1),ylab ="")
lines(density(exp(outEM2[1000:10000,1])),col=2)
abline(v=0.2)
plot(density(exp(out2[1000:10000,2])), xlab = expression(sigma), main = expression(Delta * t == 1),ylab ="")
lines(density(exp(outEM2[1000:10000,2])),col=2)
abline(v=sqrt(0.1))

dev.off()
