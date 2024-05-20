require(here)
require(tidyverse)
require(patchwork)
require(mvtnorm)
require(MetricGraph)
require(sp)
require(nimble)
require(parallel)
require(igraph)
theme_set(
  theme_light() +
    theme(
      strip.background = element_rect(color = 'gray', fill = 'white'),
      strip.text.x = element_text(color = 'black'),
      strip.text.y = element_text(color = 'black')
    )
)

rm(list=ls())

set.seed(342345)

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


fn_optim_gc <- function(par, y, dist_mat, time_mat, n){
  cov_mat <- gc_gc_cov(dists = dist_mat,time = time_mat, 
                      pars = c(1,2,exp(par[2]),exp(par[3]),1),
                      sigma2 = exp(par[1]),
                      tau2 = exp(par[4]))
  dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
}


######################################################### 
## We investigate here the sensitivity of the process  ##
## to variation of its values                          ##
#########################################################

## Network topology
line1 <- Line(rbind(c(0,0),c(1,0)))
line2 <- Line(rbind(c(0,0),c(0,1)))
line3 <- Line(rbind(c(0,1),c(-1,1)))
theta <- seq(from=pi,to=3*pi/2,length.out = 50)
line4 <- Line(cbind(sin(theta),1+ cos(theta)))
Lines <- sp::SpatialLines(list(Lines(list(line1),ID="1"),
                               Lines(list(line2),ID="2"),
                               Lines(list(line3),ID="3"),
                               Lines(list(line4),ID="4")))

nodes_coord <- tibble(x = c(0,1,0,-1),
                      y = c(0,0,1,1),
                      node = c(1,2,3,4))
graph <- metric_graph$new(lines = Lines)
graph$plot()

graph$clear_observations()

## Generate data
# locations: network nodes + n_s random locations
n_obs_on_edges <- 10
nv <- 4
ns <- n_obs_on_edges + nv
nt <- 6
PtE <- cbind(c(1,1,3,4),c(0,1,0,1))
PtE <- rbind(PtE,
             cbind(sample(1:graph$nE, size = ns - 4, replace = TRUE), # edge number
                   runif(min = 0.05,
                         max = 0.95,
                         n = ns - 4))) # distance on edge

# Simulate data with Whittle-Matern SPDE
sigma <- 1.3
alpha <- 1
range <- 0.2
u_t <- lapply(1:nt, \(x)(sample_spde(sigma = sigma, 
                                     range = range,
                                     graph = graph, PtE = PtE))) # sample different time shots
u_t_mat <- matrix(unlist(u_t), nt, ns, byrow=TRUE)
u_t_vec <- as.numeric(u_t_mat) # it should be in the correct order

data_1 <- data.frame(edge_number = PtE[,1],
                     distance_on_edge =  PtE[,2],
                     u_1 = u_t[[1]])
graph$add_observations(data = data_1, normalized = TRUE)
graph$plot(data = "u_1", vertex_size = 0)


data_gdist <- graph$compute_geodist_PtE(PtE = PtE, 
                                        normalized = TRUE, 
                                        include_vertices = FALSE) # vertices already included bc in the observations
dist_mat <- data_gdist %x% matrix(1,nt,nt) # space distances
n <- ns * nt # total number of data points

time_mat <- as.matrix(dist(rep(1:nt,ns))) #as.matrix(dist(rep(1,n))) # as.matrix(dist(runif(n, 0,1))) # time distances

## Model fit (space-time)
tau2_true <- 0.1
sig2_true <- 0.9
cov_pars_true <- c(1, 2, 20, 0.2, 1)

gc_network <- optim(par = rnorm(4), #c(log(0.9), log(20), log(.2), log(1)),
                    fn = fn_optim_gc,
                    control = list(fnscale = -1),
                    y = u_t_vec,
                    dist_mat = dist_mat,
                    time_mat = time_mat,
                    n = n)

gc_network$par

## Isotropic GP fit
res_exp <- graph_lme(u_1 ~ -1, graph = graph, model = list(type = "isoCov"))
res_exp$coeff$random_effects
res_exp$loglik


## Witthle-Matern fit
res_WM <- graph_lme(u_1 ~ -1, graph = graph, model = 'WM1')
res_WM$coeff$random_effects
res_WM$loglik

## Laplacian fit
# res_GL <- graph_lme(u_1 ~ -1, graph = graph, model = list(type = "graphLaplacian"))
# res_GL$coeff

# write code: change value in vertex (magnitude), and save parameters
# should data be only on vertexes? no, probably it's fine if it is not
# compare res_exp and res_WM
# can we make a network with all edges of length 1? it seems yes, by modifying the filed graph$edge_length


# network nodes are the first nv entries
y <- u_t[[1]] # original data
res_sens <- tibble(node = numeric(),
                   h = numeric(),
                   loglik = numeric(),
                   tau = numeric(),
                   kappa = numeric(),
                   model = character())
for(v in 1:nv){
  for(h in seq(-2,2,by=0.1)){
    y_h <- y
    y_h[v] <- y[v] + h
    data_h <- data.frame(edge_number = PtE[,1],
                         distance_on_edge =  PtE[,2],
                         y_h = y_h)
    graph$clear_observations()
    graph$add_observations(data = data_h, normalized = TRUE)
    # Isotropic GP fit
    res_exp <- graph_lme(y_h ~ -1, graph = graph, model = list(type = "isoCov"))
    res_exp$coeff$random_effects
    # Witthle-Matern fit
    res_WM <- graph_lme(y_h ~ -1, graph = graph, model = 'WM1')
    res_WM$coeff$random_effects
    # Saving
    res_sens <- add_row(res_sens,
                        node = rep(v,2),
                        h = rep(h,2),
                        loglik = c(res_exp$loglik,
                                   res_WM$loglik),
                        tau = c(res_exp$coeff$random_effects[1],
                                res_WM$coeff$random_effects[1]),
                        kappa = c(res_exp$coeff$random_effects[2],
                                  res_WM$coeff$random_effects[2]),
                        model = c('isoCov', 
                                  'WM1'))
  }
}


res_sens %>%
  filter(model == "isoCov") %>%
  pivot_longer(cols = c('loglik','tau','kappa'), 
               names_to = 'parameter', 
               values_to = 'value') %>%
  ggplot(aes(x = h, y = value, group = node)) + 
  geom_line(aes(color = factor(node)), alpha = 0.8) +
  facet_grid(rows = vars(parameter),
             cols = vars(model),
             scale = "free") +
  guides(color = "none") +
  labs(color = "") +
  theme(axis.title.y = element_blank()) -> plot_sens_isocov

res_sens %>%
  filter(model == "WM1") %>%
  pivot_longer(cols = c('loglik','tau','kappa'), 
               names_to = 'parameter', 
               values_to = 'value') %>%
  ggplot(aes(x = h, y = value, group = node)) + 
  geom_line(aes(color = factor(node)), alpha = 0.8) +
  facet_grid(rows = vars(parameter),
             cols = vars(model),
             scale = "free") +
  labs(color = "") + 
  theme(axis.title.y = element_blank()) -> plot_sens_wm1

graph$plot(vertex_size = 2) + 
  geom_label(data = nodes_coord,
             aes(x = x, y = y, label = node), size = 3) +
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank()) -> plot_net


(plot_sens_isocov + plot_sens_wm1) / plot_net +
  plot_layout(heights = c(1,0.3))

#############################
#### A different network ####
# sample random point
# sample random connections
## Network topology
V <- rbind(c(0, 0),
           c(1, 0),
           c(1, 1),
           c(0, 1),
           c(-1, 1),
           c(-1, 0),
           c(0, -2))
E <- rbind(c(1, 2),
           c(2, 3),
           c(3, 4),
           c(4, 5),
           c(5, 6),
           c(6, 1),
           c(4, 1),
           c(1, 7))
graph2 <- metric_graph$new(V = V, E = E)
graph2$plot()
graph2$compute_geodist()
graph2$geo_dist
graph2$compute_resdist()
graph2$res_dist


graph3 <- graph2$clone(deep = TRUE)
graph3$edge_lengths <- rep(1,8)
graph3$plot()
graph3$compute_geodist()
graph3$geo_dist
graph3$compute_resdist()
graph3$res_dist



graph2$A(obs_to_vert = TRUE)


net_2 <- igraph::graph_from_data_frame(d = tibble(source = graph2$E[,1],
                                                  target = graph2$E[,2],
                                                  weight = graph2$edge_lengths),
                                       vertices = 1:nrow(graph2$V),
                                       directed = FALSE)


net_3 <- igraph::graph_from_data_frame(d = tibble(source = graph3$E[,1],
                                                  target = graph3$E[,2],
                                                  weight = graph3$edge_lengths),
                                       vertices = 1:nrow(graph3$V),
                                       directed = FALSE)

igraph::closeness(net_2)
igraph::closeness(net_3)

igraph::betweenness(net_2)
igraph::betweenness(net_3)

graph_pems <- metric_graph$new(lines = pems$lines)
graph_pems$plot()


###########################################
#### Graph examples (from Borovitskiy) ####
###########################################

PtE_node <- function(node, E) # function to create PtE for a given node
{
  n <- nrow(E)
  from_to <- 0
  for(i in 1:n){
    if(sum(E[i,] == node) > 0){
      if(E[i,1] == node) from_to <- 1 else from_to <- 2
      break
    }
  }
  return(c(i, from_to - 1))
}

PtE_allnodes <- function(tot_nodes, E)
{
  PtE <- matrix(NA, tot_nodes, 2)
  for(i in 1:tot_nodes){
    PtE[i,] <- PtE_node(i,E)
  }
  return(PtE)
}

## Complete Graph (10 total nodes)
tot_nodes <- 10
V <- matrix(0, tot_nodes, 2)
for(i in 1:tot_nodes){
  V[i,] <- c(cos(2*pi/tot_nodes*(i-1)), sin(2*pi/tot_nodes*(i-1)))
}
E <- t(combn(tot_nodes,2))
comp_graph <- metric_graph$new(V = V, E = E)
comp_graph$plot()

comp_graph_all1 <- comp_graph$clone(deep = TRUE)
comp_graph_all1$edge_lengths <- rep(1,length(comp_graph_all1$edge_lengths))

net_comp <- graph_from_data_frame(d = tibble(source = comp_graph$E[,1],
                                             target = comp_graph$E[,2],
                                             weight = comp_graph$edge_lengths),
                                  vertices = 1:nrow(comp_graph$V),
                                  directed = FALSE)
comp_meta <- tibble(node = 1:tot_nodes,
                    closeness = closeness(net_comp),
                    betweenness = betweenness(net_comp),
                    degree = degree(net_comp))

net_comp_all1 <- graph_from_data_frame(d = tibble(source = comp_graph_all1$E[,1],
                                                  target = comp_graph_all1$E[,2],
                                                  weight = comp_graph_all1$edge_lengths),
                                       vertices = 1:nrow(comp_graph_all1$V),
                                       directed = FALSE)
comp_all1_meta <- tibble(node = 1:tot_nodes,
                         closeness = closeness(net_comp_all1),
                         betweenness = betweenness(net_comp_all1),
                         degree = degree(net_comp_all1))

# Data generation with Whittle-Matern SPDE
# Observations only on the edges
sigma <- 1.3
alpha <- 1
range <- 2
PtE_comp <- PtE_allnodes(tot_nodes,E)

# comp_graph$build_mesh(h = 0.01)
# u_comp_mesh <- sample_spde(range = range, 
#                            sigma = sigma, 
#                            alpha = alpha,
#                            graph = comp_graph, 
#                            type = "mesh", 
#                            method = "Q")
# comp_graph$plot_function(u_comp_mesh, 
#                          vertex_size = 0,
#                          edge_width = 0.5)
# C_comp <- spde_covariance(P = c(5, 0.99), 
#                           range = range, 
#                           sigma = sigma, 
#                           alpha = alpha,
#                           graph = comp_graph)
# comp_graph$plot_function(C_comp, vertex_size = 0, edge_width = 0.5)

u_comp <- sample_spde(sigma = sigma, 
                      range = range,
                      graph = comp_graph, 
                      PtE = PtE_comp)
data_comp <- data.frame(edge_number = PtE_comp[,1],
                        distance_on_edge =  PtE_comp[,2],
                        u = u_comp)
comp_graph$add_observations(data = data_comp, normalized = TRUE)
comp_graph$plot(data = "u", vertex_size = 0)

# 
# comp_graph_all1$build_mesh(h = 0.01)
# C_comp_all1 <- spde_covariance(P = c(5, 0.99), 
#                           range = range, 
#                           sigma = sigma, 
#                           alpha = alpha,
#                           graph = comp_graph_all1)
# comp_graph_all1$plot_function(C_comp_all1/max(C_comp_all1), vertex_size = 0, edge_width = 0.5)

u_comp_all1 <- sample_spde(sigma = sigma, 
                           range = range,
                           graph = comp_graph_all1, 
                           PtE = PtE_comp)
data_comp_all1 <- data.frame(edge_number = PtE_comp[,1],
                             distance_on_edge =  PtE_comp[,2],
                             u = u_comp_all1)
comp_graph_all1$add_observations(data = data_comp_all1, normalized = TRUE)
comp_graph_all1$plot(data = "u", vertex_size = 0)

## Star Graph (1 internal + 9 external nodes)
tot_nodes <- 10
V <- matrix(0, tot_nodes, 2)
for(i in 1:(tot_nodes-1)){
  V[i+1,] <- c(cos(2*pi/(tot_nodes-1)*(i-1)), sin(2*pi/(tot_nodes-1)*(i-1)))
}
E <- matrix(1, tot_nodes-1, 2)
E[,2] <- 1 + 1:(tot_nodes-1)
star_graph <- metric_graph$new(V = V, E = E)
star_graph$plot()

net_star <- graph_from_data_frame(d = tibble(source = star_graph$E[,1],
                                             target = star_graph$E[,2],
                                             weight = star_graph$edge_lengths),
                                  vertices = 1:nrow(star_graph$V),
                                  directed = FALSE)
star_meta <- tibble(node = 1:tot_nodes,
                    closeness = closeness(net_star),
                    betweenness = betweenness(net_star),
                    degree = degree(net_star))

# Data generation with Whittle-Matern SPDE
# Observations only on the edges
sigma <- 1.3
alpha <- 1
range <- 2
PtE_star <- PtE_allnodes(tot_nodes,E)

# star_graph$build_mesh(h = 0.01)
# u_star_mesh <- sample_spde(range = range, 
#                            sigma = sigma, 
#                            alpha = alpha,
#                            graph = star_graph, 
#                            type = "mesh", 
#                            method = "Q")
# star_graph$plot_function(u_star_mesh, 
#                          vertex_size = 0,
#                          edge_width = 0.5)
# C_star <- spde_covariance(P = c(5, 0.99), 
#                           range = range, 
#                           sigma = sigma, 
#                           alpha = alpha,
#                           graph = star_graph)
# star_graph$plot_function(C_star, vertex_size = 0, edge_width = 0.5)

u_star <- sample_spde(sigma = sigma, 
                      range = range,
                      graph = star_graph, 
                      PtE = PtE_star)
data_star <- data.frame(edge_number = PtE_star[,1],
                        distance_on_edge =  PtE_star[,2],
                        u = u_star)
star_graph$add_observations(data = data_star, normalized = TRUE)
star_graph$plot(data = "u", vertex_size = 0)


# Random Graph
# Maybe to use only _all1 to guarantee that it is an Eucledean graph.
tot_nodes <- 10
tot_edges <- 25
V <- matrix(0, tot_nodes, 2)
for(i in 1:tot_nodes){
  V[i,] <- runif(2, min = -1, max = 1)
}
E_full <- t(combn(tot_nodes,2))
E <- E_full[sample(1:nrow(E_full), tot_edges, replace = FALSE),]
rand_graph <- metric_graph$new(V = V, E = E)
rand_graph$plot()

net_rand <- graph_from_data_frame(d = tibble(source = rand_graph$E[,1],
                                             target = rand_graph$E[,2],
                                             weight = rand_graph$edge_lengths),
                                  vertices = 1:nrow(rand_graph$V),
                                  directed = FALSE)
rand_meta <- tibble(node = 1:tot_nodes,
                    closeness = closeness(net_rand),
                    betweenness = betweenness(net_rand),
                    degree = degree(net_rand))

rand_graph_all1 <- rand_graph$clone(deep = TRUE)
rand_graph_all1$edge_lengths <- rep(1,length(rand_graph_all1$edge_lengths))

net_rand_all1 <- graph_from_data_frame(d = tibble(source = rand_graph_all1$E[,1],
                                                  target = rand_graph_all1$E[,2],
                                                  weight = rand_graph_all1$edge_lengths),
                                       vertices = 1:nrow(rand_graph_all1$V),
                                       directed = FALSE)
rand_all1_meta <- tibble(node = 1:tot_nodes,
                         closeness = closeness(net_rand_all1),
                         betweenness = betweenness(net_rand_all1),
                         degree = degree(net_rand_all1))

# Data generation with Whittle-Matern SPDE
# Observations only on the edges
sigma <- 1.3
alpha <- 1
range <- 2
PtE_rand <- PtE_allnodes(tot_nodes,E)

# rand_graph$build_mesh(h = 0.01)
# u_rand_mesh <- sample_spde(range = range, 
#                            sigma = sigma, 
#                            alpha = alpha,
#                            graph = rand_graph, 
#                            type = "mesh", 
#                            method = "Q")
# rand_graph$plot_function(u_rand_mesh, 
#                          vertex_size = 0,
#                          edge_width = 0.5)
# C_rand <- spde_covariance(P = c(5, 0.001), 
#                           range = range, 
#                           sigma = sigma, 
#                           alpha = alpha,
#                           graph = rand_graph)
# rand_graph$plot_function(C_rand, vertex_size = 0, edge_width = 0.5)


u_rand <- sample_spde(sigma = sigma, 
                      range = range,
                      graph = rand_graph, 
                      PtE = PtE_rand)
data_rand <- data.frame(edge_number = PtE_rand[,1],
                        distance_on_edge =  PtE_rand[,2],
                        u = u_rand)
rand_graph$add_observations(data = data_rand, normalized = TRUE)
rand_graph$plot(data = "u", vertex_size = 0)

u_rand_all1 <- sample_spde(sigma = sigma, 
                           range = range,
                           graph = rand_graph_all1, 
                           PtE = PtE_rand)
data_rand_all1 <- data.frame(edge_number = PtE_rand[,1],
                             distance_on_edge =  PtE_rand[,2],
                             u = u_rand_all1)
rand_graph_all1$add_observations(data = data_rand_all1, normalized = TRUE)
rand_graph_all1$plot(data = "u", vertex_size = 0)



## Running experiment
run_experiment <- function(graph, u, PtE, graph_name)
{
  y <- u # original data
  nv <- length(u) # NOTE: ONLY IF OBSERVATIONS ARE ONLY ON VERTEX
  res_sens <- tibble(node = numeric(),
                     h = numeric(),
                     loglik = numeric(),
                     tau = numeric(),
                     kappa = numeric(),
                     model = character(),
                     graph_name = character())
  for(v in 1:nv){
    for(h in seq(-2,2,by=0.1)){
      y_h <- y
      y_h[v] <- y[v] + h
      data_h <- data.frame(edge_number = PtE[,1],
                           distance_on_edge =  PtE[,2],
                           y_h = y_h)
      graph$clear_observations()
      graph$add_observations(data = data_h, normalized = TRUE)
      # Isotropic GP fit
      res_exp <- graph_lme(y_h ~ -1, graph = graph, model = list(type = "isoCov"))
      # res_exp$coeff$random_effects
      # Witthle-Matern fit
      res_WM <- graph_lme(y_h ~ -1, graph = graph, model = 'WM1')
      # res_WM$coeff$random_effects
      # Saving
      res_sens <- add_row(res_sens,
                          node = rep(v,2),
                          h = rep(h,2),
                          loglik = c(res_exp$loglik,
                                     res_WM$loglik),
                          tau = c(res_exp$coeff$random_effects[1],
                                  res_WM$coeff$random_effects[1]),
                          kappa = c(res_exp$coeff$random_effects[2],
                                    res_WM$coeff$random_effects[2]),
                          model = c('isoCov', 
                                    'WM1'),
                          graph_name = rep(graph_name, 2))
    }
  }
  return(res_sens)
}

graph_list <- list(comp_graph, comp_graph_all1, star_graph, rand_graph_all1)
u_list <- list(u_comp, u_comp_all1, u_star, u_rand_all1)
PtE_list <- list(PtE_comp, PtE_comp, PtE_star, PtE_rand)
graph_name_list <- list('Complete', 'Complete_all1', 'Star', 'Random_all1')

results <- mcmapply(graph = graph_list,
                    u = u_list,
                    PtE = PtE_list,
                    graph_name = graph_name_list,
                    FUN = run_experiment,
                    mc.cores = 2,
                    SIMPLIFY = FALSE)

results %>%
  map_dfr(.f = ~ .) -> data_results

# save.image(file = here("sensitivity", "sensitivity.Rdata"))
# load(here("sensitivity", "sensitivity.Rdata"))

rand_meta %>%
  mutate(closeness_cum = cume_dist(closeness),
         betweenness_med = betweenness/median(betweenness),
         bet2 = (betweenness - median(betweenness))/diff(range(betweenness)))
  
  
  
comp_meta %>%  
  mutate(graph_name = "Complete") %>%
  bind_rows(comp_all1_meta %>%
              mutate(graph_name = "Complete_all1")) %>%
  bind_rows(star_meta %>%
              mutate(graph_name = "Star")) %>%
  bind_rows(rand_all1_meta %>%
              mutate(graph_name = "Random_all1")) %>%
  mutate(closeness_trasf =  (closeness - median(closeness))/diff(range(closeness)),
         betweenness_trasf =  (betweenness - median(betweenness))/diff(range(betweenness)),
         degree_trasf = (degree - median(degree))/diff(range(degree)),
         .by = graph_name) %>%
  replace_na(replace = list(degree_trasf = 0.5,
                            betweenness_trasf = 0.5,
                            closeness_trasf = 0.5))-> net_meta



data_results %>%
  full_join(net_meta, by = c('node', 'graph_name')) %>%
  pivot_longer(cols = c('loglik','tau','kappa'),
               names_to = 'parameter') %>%
  filter(model == 'WM1') %>%
  mutate_at(vars(value),  \(x){ifelse(x<20,x,20)}) %>%
  ggplot(aes(x = h, y = value, group = node)) +
  geom_line(aes(color = betweenness_trasf), alpha = 0.7) +
  scale_color_gradient2(low = '#4575b4', high = '#b2182b', mid = 'gray') +
  facet_grid(cols = vars(graph_name), 
             rows = vars(parameter),
             scale = 'free') +
  theme(axis.title = element_blank()) -> plot1
(comp_graph$plot(data = "u", vertex_size = 0, data_size = 4) |
comp_graph_all1$plot(data = "u", vertex_size = 0, data_size = 4) |
rand_graph_all1$plot(data = "u", vertex_size = 0, data_size = 4) |
star_graph$plot(data = "u", vertex_size = 0, data_size = 4) +
  plot_layout(guides = "collect") ) -> plot2

plot1 / plot2
# CHECK IF ORDERING ARE ALL CONSISTENT
# data only on nodes?
# should we average over different simulations?

