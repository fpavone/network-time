require(here)
require(tidyverse)
require(mvtnorm)
require(MetricGraph)
require(sp)
require(nimble)
require(parallel)

rm(list=ls())

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

fn_optim_gc <- function(par,y,dist_mat,time_mat,tau2){
  cov_mat = gc_gc_cov(dists = dist_mat,time = time_mat, 
                      pars = c(1,2,exp(par[2]),exp(par[3]),1),
                      sigma2 = exp(par[1]),
                      tau2 = tau2)
  dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
}
## Network generation
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 50)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines <- sp::SpatialLines(list(Lines(list(line1),ID="1"),
                               Lines(list(line2),ID="2"),
                               Lines(list(line3),ID="3"),
                               Lines(list(line4),ID="4")))


graph <- metric_graph$new(lines = Lines)
graph$plot()

## Plotting covariance function of Porcu et al.
# - fix locations/mesh
# - fix time shots
# - compute covariance
# - plot covariance


# Transoform mesh to data
graph$build_mesh(h = 0.1)
graph$plot(mesh=TRUE)

graph$mesh
lon <- graph$mesh$V[, 1]
lat <- graph$mesh$V[, 2]
mesh_in_space <- cbind(lon, lat)
points <- SpatialPointsDataFrame(coords = mesh_in_space,
                                 data = data.frame(y = rep(0,nrow(mesh_in_space))))

graph$add_observations(points)
graph$plot(data = "y", vertex_size = 0)


ns <- nrow(dist_mat_s) # number of space points

n.obs <- 5
obs.loc <- cbind(sample(1:graph$nE, n.obs, replace = TRUE), runif(n.obs))
obs.lonlat <- graph$coordinates(PtE = obs.loc, normalized = TRUE)
obs <- rnorm(n.obs)
points <- SpatialPointsDataFrame(coords = obs.lonlat,
                                 data = data.frame(y = obs))
graph$add_observations(points)
graph$plot(data = "y", vertex_size = 0)

graph$coordinates(XY = matrix(c(0, 0.5), 1,2)) 
#SOMETHING DOES NOT WORK PROPERLY



########################################
## Generate data and fit Porcu et al. ##
########################################

graph$plot()
graph$clear_observations()

n.obs <- 10
data <- data.frame(edge_number = sample(1:graph$nE, n.obs, replace = TRUE),
                   distance_on_edge = runif(n.obs),
                   y = rep(0, n.obs)) # fake observations, only for locations
graph$add_observations(data = data, normalized = TRUE)
graph$plot(data = "y", vertex_size = 0)

graph$compute_geodist()
graph$geo_dist
data_gdist <- graph$geo_dist[[1]][-c(1:4),-c(1:4)]

# hyperparameters
tau2_true <- 0.1
sig2_true <- 0.9
cov_pars_true <- c(1, 2, 20, 0.2, 1)

nt <- 3 # number of time instances
ns <- nrow(data)

dist_mat <- data_gdist %x% matrix(1,nt,nt) # space distances
n <- ns * nt # total number of data points

time_mat <- as.matrix(rdist(runif(n, 0,1))) # time distances

cov_mat <- gc_gc_cov(dists = dist_mat,
                     time = time_mat, 
                     pars = cov_pars_true,
                     sigma2 = sig2_true,
                     tau2 = tau2_true)

dmvnorm(rnorm(n),rep(0,n),cov_mat,log = TRUE)
y <- rnorm(n)
gc_network <- optim(par = c(log(0.9), log(20), log(.2)),
                    fn = fn_optim_gc,
                    control = list(fnscale = -1),
                    y = y,
                    dist_mat = dist_mat,
                    time_mat = time_mat,
                    tau2 = tau2_true)

##################################################
## Sample data from Whittle-Matern SPDE and fit ##
##################################################
graph$clear_observations()

n.obs <- 10
nt <- 3
PtE <- cbind(sample(1:graph$nE, size = n.obs, replace = TRUE), # edge number
             runif(n.obs)) # distance on edge

sigma <- 1.3
alpha <- 1
range <- 0.2
u_t <- lapply(1:nt, \(x)(sample_spde(sigma = sigma, 
                   range = range,
                   graph = graph, PtE = PtE))) # sample different time shots
u_t_mat <- matrix(unlist(u_t),3,10,byrow=TRUE)
u_t_vec <- as.numeric(u_t_mat) # it should be in the correct order

data_1 <- data.frame(edge_number = PtE[,1],
                     distance_on_edge =  PtE[,2],
                     u_1 = u_t[[1]])
# I just add data for t=1. I need location on graph 
graph$add_observations(data = data_1, normalized = TRUE)
graph$plot(data = "u_1", vertex_size = 0)

# Computing geodesic distance (it's possible to compute also the resistance distance)
graph$compute_geodist()
data_gdist <- graph$geo_dist[[1]][-c(1:4),-c(1:4)]

# hyperparameters
tau2_true <- 0.1
sig2_true <- 0.9
cov_pars_true <- c(1, 2, 20, 0.2, 1)

ns <- n.obs

dist_mat <- data_gdist %x% matrix(1,nt,nt) # space distances
n <- ns * nt # total number of data points

time_mat <- as.matrix(rdist(runif(n, 0,1))) # time distances

cov_mat <- gc_gc_cov(dists = dist_mat,
                     time = time_mat, 
                     pars = cov_pars_true,
                     sigma2 = sig2_true,
                     tau2 = tau2_true)

dmvnorm(u_t_vec,rep(0,n),cov_mat,log = TRUE)

# need to flat u_t in order to create y vector with the correct ordering


gc_network <- optim(par = c(log(0.9), log(20), log(.2)),
                    fn = fn_optim_gc,
                    control = list(fnscale = -1),
                    y = u_t_vec,
                    dist_mat = dist_mat,
                    time_mat = time_mat,
                    tau2 = tau2_true)

gc_network$par



