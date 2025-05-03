#N(0,1) target
#Random walk Metroplis

mh<-function(iters=11000,v=1,sigma_sq=1) 
{
  vec=vector("numeric", iters)
  x=0
  targ_old=target(x,sigma_sq)
  vec[1]=x
  for (i in 2:iters) {
    innov=rnorm(1,0,v)
    can=x+innov
    targ_new=target(can,sigma_sq)
    aprob=targ_new/targ_old
    u=runif(1)
    if (u < aprob) { 
      x=can
      targ_old=targ_new
    }
    vec[i]=x
  }
  return(vec)
}

mh_correlated<-function(iters=11000,v=1,sigma_sq=1,rho=0.95,u0=0) 
{
  vec=vector("numeric", iters)
  x=0
  u=u0
  targ_old=correlated_target(x,sigma_sq,u)
  vec[1]=x
  for (i in 2:iters) {
    innov=rnorm(1,0,v)
    can=x+innov
    u_can = rnorm(1,rho*u,sqrt(1-rho^2))
    targ_new=correlated_target(can,sigma_sq,u_can)
    aprob=targ_new/targ_old
    w=runif(1)
    if (w < aprob) { 
      x=can
      u=u_can
      targ_old=targ_new
    }
    vec[i]=x
  }
  return(vec)
}

#Exact target density evaluation
# target<-function(x)
# {
#   return(dnorm(x))
# }

#Target density evaluation with lognormal noise
target <- function(x,sigma_sq)
{
  u = rlnorm(1,-1/(2*sigma_sq),sqrt(sigma_sq))
  #Chosen logmean = -1/(2*sigma_sq) to ensure constant expectation but keep variance parameterised
  return(dcauchy(x, location = 0, scale = 1)*u)
}

#Adjusted for correlated pmmh
correlated_target <- function(x,sigma_sq,u)
{
  out <- dcauchy(x, location = 0, scale =1)*exp(-sigma_sq/2+sqrt(sigma_sq)*u)
  return(out)
  #leaving out g(u) which is assumed to be generated from N(0,1) so cancels in ratio?
}


#Run sampler - argument is a list of the desired variances of the noise
run_sampler <- function(variances){
  #par(mfrow=c(length(variances),4))
  xlim <- c(-10,10)
  breaks <- breaks <- c(-50, -10, seq(-10, 10, by = 0.25), 10, 50)
  for (i in 1:length(variances)){
    out = mh(sigma_sq=variances[i])[1001:11000]
    plot(ts(out),main = bquote(sigma^2 == .(variances[i])), cex.lab=5,cex.axis = 5, cex.main = 5,ylim = c(-20,25), xlab = "Iteration", ylab = expression(theta))
    hist(out,breaks = breaks, xlim = xlim,freq=FALSE,main = "uncorrelated",cex.main = 4,cex.lab=2,cex.axis = 2, ylim = c(0,0.9), xlab = expression(theta))
    lines(seq(-4,4,0.001),dnorm(seq(-4,4,0.001)),type="l")
    
    out = mh_correlated(sigma_sq=variances[i])[1001:11000]
    plot(ts(out),main = bquote(sigma^2 == .(variances[i])),cex.lab=5,cex.axis = 5, cex.main = 5, ylim = c(-20,25), xlab = "Iteration", ylab = expression(theta))
    hist(out,breaks = breaks, xlim = xlim,freq=FALSE,main = "correlated",cex.main = 4,cex.lab=2,cex.axis = 2, ylim = c(0,0.9), xlab = expression(theta))
    lines(seq(-4,4,0.001),dnorm(seq(-4,4,0.001)),type="l")
  }
}
#previously cex.main = 5 and cex.lab = 3

output_dir <- "C:\\Users//HP\\OneDrive - Durham University\\Course\\Project IV\\Plots\\Toy example new"
#output_dir <- "C:\\Users\\qnkp75\\OneDrive - Durham University\\Course\\Project IV\\Plots"

#png(file=file.path(output_dir, "pmmh toy example_%d.png"), width = 1000, height = 1000) 
pdf(file=file.path(output_dir, "pmmh toy example_%d.pdf"), onefile = F,width = 6, height = 6) 
set.seed(5)
#seed 2 okay
par(mfrow = c(1,1),mar=c(5, 5, 6, 2),cex.axis = 2, cex.lab = 2, cex.main = 4 )
run_sampler(c(0.01,5,20))
dev.off()

par(mfrow=c(3,1),cex.axis=1.0,cex=1.4,mgp=c(2.1,1,0))
plot(density(mh(sigma_sq=0.01)[1001:11000],adjust=2.0),main="",xlab=expression(log~theta[1]),ylab="")
plot(density(mh(sigma_sq=0.01)[1001:11000],adjust=2.0),main="",xlab=expression(log~theta[1]),ylab="")
plot(density(mh(sigma_sq=0.01)[1001:11000],adjust=2.0),main="",xlab=expression(log~theta[1]),ylab="")
