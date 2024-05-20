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

source(here('sensitivity','helper_fun.R'))

## Complete Graph (10 total nodes)
tot_nodes <- 10
V_comp <- matrix(0, tot_nodes, 2)
for(i in 1:tot_nodes){
  V_comp[i,] <- c(cos(2*pi/tot_nodes*(i-1)), sin(2*pi/tot_nodes*(i-1)))
}
E_comp <- t(combn(tot_nodes,2))
comp_graph <- metric_graph$new(V = V_comp, E = E_comp)
# comp_graph$plot()

comp_graph_all1 <- comp_graph$clone(deep = TRUE)
comp_graph_all1$edge_lengths <- rep(1,length(comp_graph_all1$edge_lengths))

## Star Graph (1 internal + 9 external nodes)
tot_nodes <- 10
V_star <- matrix(0, tot_nodes, 2)
for(i in 1:(tot_nodes-1)){
  V_star[i+1,] <- c(cos(2*pi/(tot_nodes-1)*(i-1)), sin(2*pi/(tot_nodes-1)*(i-1)))
}
E_star <- matrix(1, tot_nodes-1, 2)
E_star[,2] <- 1 + 1:(tot_nodes-1)
star_graph <- metric_graph$new(V = V_star, E = E_star)
# star_graph$plot()

# Random Graph
# Maybe to use only _all1 to guarantee that it is an Eucledean graph.
tot_nodes <- 10
tot_edges <- 25
V_rand <- matrix(0, tot_nodes, 2)
for(i in 1:tot_nodes){
  V_rand[i,] <- runif(2, min = -1, max = 1)
}
E_full <- t(combn(tot_nodes,2))
E_rand <- E_full[sample(1:nrow(E_full), tot_edges, replace = FALSE),]
rand_graph <- metric_graph$new(V = V_rand, E = E_rand)
# rand_graph$plot()

rand_graph_all1 <- rand_graph$clone(deep = TRUE)
rand_graph_all1$edge_lengths <- rep(1, length(rand_graph_all1$edge_lengths))

## Experiment
run_experiment <- function(graph, PtE, tau, kappa, alpha, nsim)
{
  res <- matrix(NA, nrow = nsim, ncol = 4)
  colnames(res) <- c('isoCov_tau',
                     'isoCov_kappa',
                     'WM1_tau',
                     'WM1_kappa')
  for(i in 1:nsim){
    u <- sample_spde(tau = tau,
                     kappa = kappa,
                     graph = graph,
                     alpha = alpha,
                     sigma_e = 0,
                     PtE = PtE)
    data <- data.frame(edge_number = PtE[,1],
                       distance_on_edge =  PtE[,2],
                       u = u)
    graph$clear_observations()
    graph$add_observations(data = data, 
                           normalized = TRUE)
    res_exp <- graph_lme(u ~ -1, graph = graph, model = list(type = "isoCov"))
    res_WM <- graph_lme(u ~ -1, graph = graph, model = 'WM1')
    res[i,] <- c(res_exp$coeff$random_effects, res_WM$coeff$random_effects)
  }
  return(res)
}

graph_list <- list(comp_graph, comp_graph_all1, star_graph, rand_graph_all1)
names(graph_list) <- c('Complete', 'Complete_all1', 'Star', 'Random_all1')
# n_obs_edge <- 10
# n_obs <- 100
# PtE_list <- list(PtE_obs_edge(n_obs_edge, E_comp),
#                  PtE_obs_edge(n_obs_edge, E_comp),
#                  PtE_obs_edge(n_obs_edge, E_star),
#                  PtE_obs_edge(n_obs_edge, E_rand))
# PtE_list <- list(PtE_obs_tot(n_obs, comp_graph$E, comp_graph$edge_lengths),
#                  PtE_obs_tot(n_obs, comp_graph_all1$E, comp_graph_all1$edge_lengths),
#                  PtE_obs_tot(n_obs, star_graph$E, star_graph$edge_lengths),
#                  PtE_obs_tot(n_obs, rand_graph_all1$E, rand_graph_all1$edge_lengths))


# mapply(graph = graph_list,
#        PtE = PtE_list,
#        FUN = function(graph,PtE){graph$compute_geodist_PtE(PtE = PtE, include_vertices = FALSE)})

n_obs <- tot_nodes + 10*sapply(graph_list,
                               FUN = \(x){ceiling(sum(x$edge_lengths))})
names(n_obs) <- c('Complete','Complete_all1','Star','Random_all1')
PtE_list <- list(PtE_obs_nodes_edge(n_obs['Complete'], comp_graph$E, comp_graph$edge_lengths),
                 PtE_obs_nodes_edge(n_obs['Complete_all1'], comp_graph_all1$E, comp_graph_all1$edge_lengths),
                 PtE_obs_nodes_edge(n_obs['Star'], star_graph$E, star_graph$edge_lengths),
                 PtE_obs_nodes_edge(n_obs['Random_all1'], rand_graph_all1$E, rand_graph_all1$edge_lengths))
# PtE_list2 <- lapply(PtE_list,
#                     function(x){
#                       index <- sort.int(x[,1], index.return = TRUE)$ix
#                       return(x[index,])
#                     })
# - tau: precision. Smaller -> higher variance
# - kappa: length scale. Smaller -> larger length scale
# - alpha: fixed to 1 (nu = 1/2)

tau_true <- 1
kappa_true <- 1
alpha <- 1
nsim <- 20

results <- mcmapply(graph = graph_list,
                    PtE = PtE_list,
                    FUN = run_experiment,
                    MoreArgs = list(tau = tau_true,
                                    kappa = kappa_true,
                                    alpha = alpha,
                                    nsim = nsim),
                    mc.cores = 3,
                    SIMPLIFY = FALSE,
                    USE.NAMES = TRUE)

# comment <- '
# Observations on edges, on each edge n_obs_edge observations.
# '
# comment <- '
# Observations on edges,tot obs.
# '
comment <- 'Testing PtE_obs_nodes_edge'

save.image(here('sensitivity',
                paste('stability',
                      # paste('_nobs', n_obs, sep = ''),
                      '_test_PtE_obs_nodes_edge',
                      # '_NOPAR_',
                      format(Sys.Date(),'%d%m%y'),
                      '.Rdata', sep = '')))



