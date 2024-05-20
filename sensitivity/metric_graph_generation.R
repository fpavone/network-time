require(here)
require(tidyverse)
require(mvtnorm)
require(MetricGraph)
require(sp)

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


n_obs <- 100

data <- tibble(edge_number = sample(1:graph$nE, 
                                    size = n_obs, 
                                    replace = TRUE),
               distance_on_edge = runif(n_obs),
               y = rnorm(n_obs))

graph$add_observations(data = data, normalized = TRUE)

graph$plot(data = "y")

data2 <- tibble(edge_number = c(1,2,3,4),
               distance_on_edge = rep(0.1,4),
               z = rnorm(4))
graph$add_observations(data = data2, normalized = TRUE)
graph$plot(data = "z")


## Sampling the process
sigma <- 1.3
range <- 0.15 # range parameter
sigma_e <- 0.1

n.obs <- 200
obs.loc <- cbind(sample(1:graph$nE, n.obs, replace=TRUE), 
                 runif(n.obs))
u <- sample_spde(range = range, 
                 sigma = sigma, 
                 alpha = 1,
                 graph = graph, 
                 PtE = obs.loc, 
                 method = "Q")
graph$add_observations(data = data.frame(edge_number = obs.loc[, 1],
                                         distance_on_edge = obs.loc[, 2],
                                         u = u),
                       normalized = TRUE)

graph$plot(data = "u")

graph$build_mesh(h = 0.01)
graph$plot(mesh=TRUE)


C <- spde_covariance(c(1, 0.01), range = range, sigma = sigma, alpha = 1,
                     graph = graph)
graph$plot_function(C, vertex_size = 0, edge_width = 0.5)

spde_covariance(P, kappa, tau, range, sigma, alpha, graph)


precisions <- spde_precision(kappa = range, tau = 1/sigma^2, alpha = 1, graph, BC = 1, build = TRUE)

1/diag(precisions)

# look at variance at node? how does it change with network topology changes?
# this can be done also with the prior

# try to simulate data and fit Porcu et al. model
# try to use Porcu et al. covariance to simulate data on the graph
# the link is in computing the distances on the graph, in the same way they 
# are expected by Porcu et al. Is it possible?



graph$compute_geodist()
graph$geo_dist

graph$clear_observations()
graph$compute_geodist()
graph$geo_dist

graph$add_observations(data = data.frame(edge_number = c(1,1,3,2,4),
                                         distance_on_edge = runif(5),
                                         u = rnorm(5)),
                       normalized = TRUE)
graph$compute_resdist()
graph$res_dist





