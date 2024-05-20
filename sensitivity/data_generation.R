require(here)
require(tidyverse)
require(mvtnorm)
require(nimble)
require(parallel)

set.seed(1)

## Functions
gc_gc_cov <- nimbleFunction(     
  run = function(dists = double(2),time = double(2), pars = double(1),
                 sigma2 = double(0),tau2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        temp1 <- (1 + (time[i,j]/pars[4])^(2*pars[1]))
        temp <- sigma2/temp1^(pars[5] + 1) * (1 + (dists[i,j]/pars[3])/temp1)^(-pars[2])
        
        result[i, j] <- temp
        result[j, i] <- temp
      }
    }
    for(i in 1:(n)){
      
      temp1 <- (1 + (time[i,i]/pars[4])^(2*pars[1]))
      temp <- sigma2/temp1^(pars[5] + 1) * (1 + (dists[i,i]/pars[3])/temp1)^(-pars[2])
      
      result[i, i] <- temp + tau2
    }
    return(result)
  })

fn_optim_gc <- function(par, y, dist_mat, time_mat, tau2){
  cov_mat = gc_gc_cov(dists = dist_mat,
                      time = time_mat, 
                      pars = c(1,2,exp(par[2]),exp(par[3]),1),
                      sigma2 = exp(par[1]),
                      tau2 = tau2_true)
  dmvnorm(y, rep(0,n), cov_mat, log = TRUE)
}

## Space(network)-time domain data
nt <- 3 # number of time measurements at each location

# Matrix of space distances: geodesic and Euclidean
dist_mat_s <- as.matrix(read.csv(here("publication_code",
                                      "simulation",
                                      "dist_mat_sim.csv"), header = FALSE)) 
dist_mat_euc_s <- as.matrix(read.csv(here("publication_code",
                                          "simulation",
                                          "dist_mat_euc_sim.csv"), header = FALSE)) 

dist_mat_s <- dist_mat_s[1:5,1:5] # submatrices to ease computations
dist_mat_euc_s <- dist_mat_euc_s[1:5,1:5]

ns <- nrow(dist_mat_s) # number of space points

# Expanded matrix of space distances, 
# with replicates corresponding to the different time shots
dist_mat <- dist_mat_s %x% matrix(1,nt,nt) 
dist_mat_euc <- dist_mat_euc_s %x% matrix(1,nt,nt) 


n <- ns * nt # total number of data points

# Matrix of time distances
time_mat <- as.matrix(rdist(runif(n, 0,1)))

## Parameters
tau2_true <- 0.1
sig2_true <- 0.9
cov_pars_true <- c(1, 2, 20, 0.2, 1)

cov_mat = gc_gc_cov(dists = dist_mat,
                    time = time_mat, 
                    pars = cov_pars_true,
                    sigma2 = sig2_true,
                    tau2 = tau2_true)

## Data generation
y <- c(scale(mvrnorm(1,rep(0,n),cov_mat), scale = FALSE))


## Model fit (MLE)
gc_network <- optim(par = c(log(0.9), log(20), log(.2)),
                    fn = fn_optim_gc,
                    control = list(fnscale = -1),
                    y = y,
                    dist_mat = dist_mat,
                    time_mat = time_mat,
                    tau2 = tau2_true)






