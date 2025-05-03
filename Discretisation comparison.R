set.seed(123)

discretise_time <- function(end_time,delta_t){
  return(seq(from = 0, to = end_time, by = delta_t))
}

GBM_exact <- function(end_time=1, delta_t, x0=1, mu, sigma_sq){
  X <- vector(length = end_time/delta_t+1)
  if (x0 > 0 ){
    X[1] <- x0
  }
  else{
    warning("can't intialise with x0 <= 0")
  }
  for (i in 2:length(X)){
    X[i] <- rlnorm(1,
                   log(X[i-1])+(mu-sigma_sq/2)*delta_t,
                   sqrt(sigma_sq*delta_t))
                   #use sd not variance
  }
  return(X)
}

GBM_EM <- function(end_time=1, delta_t, x0=1, mu, sigma_sq){
  X <- vector(length = end_time/delta_t+1)
  if (x0 != 0 ){
    X[1] <- x0
  }
  else{
    warning("can't intialise with x0 = 0 as can't evaluate log of 0")
  }
  for (i in 2:length(X)){
    X[i] <- abs(rnorm(1,
                  X[i-1]+mu*X[i-1]*delta_t,
                  sqrt(sigma_sq*X[i-1]*delta_t)))
    #Adding absolute value is effectively a reflecting barrier at 0
  }
  return(X)
}

matrix_of_runs <- function(is_exact,iters,end_time, delta_t,x0,mu,sigma_sq){
  mat <- matrix(nrow=iters,ncol = end_time/delta_t+1)
  for (i in 1:iters){
    if (is_exact){
      mat[i,] <- GBM_exact(end_time,delta_t,x0,mu,sigma_sq)
    }
    else{
      mat[i,] <- GBM_EM(end_time,delta_t,x0,mu,sigma_sq)
    }
  }
  return(mat)
}


exact_mat20 <- matrix_of_runs(TRUE,1000,20,0.001,1,0.1,0.2)
em_0.1 <- matrix_of_runs(FALSE,1000,20,0.1,1,0.1,0.2)
em_1 <- matrix_of_runs(FALSE,1000,20,1,1,0.1,0.2)
em_2 <- matrix_of_runs(FALSE,1000,20,2,1,0.1,0.2)
em_4 <- matrix_of_runs(FALSE,1000,20,4,1,0.1,0.2)

#output_dir <- "C:\\Users\\qnkp75\\OneDrive - Durham University\\Course\\Project IV\\Plots"
output_dir <- "C:\\Users\\HP\\OneDrive - Durham University\\Course\\Project IV\\Plots"
pdf(file=file.path(output_dir, "Discretisation comparison.pdf"), onefile=F, width = 5, height = 5) 


par(mfrow = c(1,1),mar = c(3.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0), mgp = c(2, 0.5, 0))
mean_traj_exact<- exp(0.1*0.1*0:200)
lower_quantile_exact <- apply(exact_mat20[,1+0:200*100],2,quantile, probs = 0.05)
upper_quantile_exact <- apply(exact_mat20[,1+0:200*100],2,quantile, probs = 0.95)
plot(ts(mean_traj_exact, start = 0,deltat = 0.1),ylim = c(0,35),col = "black", xlab = "t", ylab = expression(X[t]), main = "", xlim = c(0,20), lwd = 2)
lines(ts(lower_quantile_exact,start = 0,deltat = 0.1),col = "black", lty =2, lwd =2)
lines(ts(upper_quantile_exact, start = 0, deltat = 0.1),col = "black",lty =2, lwd =2)
mat_list <- list(
  em_0.1 = em_0.1, 
  em_1 = em_1, 
  em_2 = em_2, 
  em_4 = em_4
)

for (name in names(mat_list)) {
  col <- switch(name,
                "em_0.1"= "blue",
                "em_1"="cyan",
                "em_2" ="orange",
                "em_4" = "red")
  delta_t <- switch(name,
                    "em_0.1"= 0.1,
                    "em_1"=1,
                    "em_2" =2,
                    "em_4" = 4)
  lines(ts(apply(mat_list[[name]], 2, mean), start = 0, deltat = delta_t),col = col)
  lines(ts(apply(mat_list[[name]], 2, quantile, probs = 0.05),start = 0,deltat = delta_t),col = col,lty=2)
  lines(ts(apply(mat_list[[name]], 2, quantile, probs = 0.95),start = 0,deltat = delta_t),col = col,lty=2)
}
legend("topleft", 
       legend = c("Exact", 
                  expression(EM ~ Delta * t == 0.1), 
                  expression(EM ~ Delta * t == 1), 
                  expression(EM ~ Delta * t == 2), 
                  expression(EM ~ Delta * t == 4)), 
       col = c("black", "blue", "cyan", "orange", "red"), 
       lty = 1,
       text.width = max(strwidth(c("Exact", "EM Δt = 0.1", "EM Δt = 1", "EM Δt = 2", "EM Δt = 4"))))

dev.off()
