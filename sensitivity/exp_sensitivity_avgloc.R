require(here)
require(tidyverse)
require(patchwork)
require(mvtnorm)
require(MetricGraph)
require(sp)
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

# args <- commandArgs(trailingOnly = TRUE)
# n_obs_notnodes <- as.numeric(args[1])

options(mc.cores = 10)

tot_nodes <- 10

## Complete Graph (10 total nodes)
V_comp <- matrix(0, tot_nodes, 2)
for(i in 1:tot_nodes){
  V_comp[i,] <- c(cos(2*pi/tot_nodes*(i-1)), sin(2*pi/tot_nodes*(i-1)))
}
E_comp <- t(combn(tot_nodes,2))
comp_graph <- metric_graph$new(V = V_comp, E = E_comp)
# comp_graph$plot()

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

## Star Graph (1 internal + 9 external nodes)
# tot_nodes <- 10
V_star <- matrix(0, tot_nodes, 2)
for(i in 1:(tot_nodes-1)){
  V_star[i+1,] <- c(cos(2*pi/(tot_nodes-1)*(i-1)), sin(2*pi/(tot_nodes-1)*(i-1)))
}
E_star <- matrix(1, tot_nodes-1, 2)
E_star[,2] <- 1 + 1:(tot_nodes-1)
star_graph <- metric_graph$new(V = V_star, E = E_star)
# star_graph$plot()

net_star <- graph_from_data_frame(d = tibble(source = star_graph$E[,1],
                                             target = star_graph$E[,2],
                                             weight = star_graph$edge_lengths),
                                  vertices = 1:nrow(star_graph$V),
                                  directed = FALSE)
star_meta <- tibble(node = 1:tot_nodes,
                    closeness = closeness(net_star),
                    betweenness = betweenness(net_star),
                    degree = degree(net_star))

# Random Graph
# Maybe to use only _all1 to guarantee that it is an Eucledean graph.
# tot_nodes <- 10
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

net_rand_all1 <- graph_from_data_frame(d = tibble(source = rand_graph_all1$E[,1],
                                                  target = rand_graph_all1$E[,2],
                                                  weight = rand_graph_all1$edge_lengths),
                                       vertices = 1:nrow(rand_graph_all1$V),
                                       directed = FALSE)
rand_all1_meta <- tibble(node = 1:tot_nodes,
                         closeness = closeness(net_rand_all1),
                         betweenness = betweenness(net_rand_all1),
                         degree = degree(net_rand_all1))


# Power-law graph
g <- sample_graph_igraph(n_nodes = tot_nodes,
                         degree_dist = \(x){x^{-1}})

powerlaw_graph <- igraph2metric(g)
powerlaw_graph_all1 <- powerlaw_graph$clone(deep = TRUE)
powerlaw_graph_all1$edge_lengths <- rep(1, length(powerlaw_graph$edge_lengths))

net_powerlaw_all1 <- graph_from_data_frame(d = tibble(source = powerlaw_graph_all1$E[,1],
                                                      target = powerlaw_graph_all1$E[,2],
                                                      weight = powerlaw_graph_all1$edge_lengths),
                                           vertices = 1:nrow(powerlaw_graph_all1$V),
                                           directed = FALSE)
powerlaw_all1_meta <- tibble(node = 1:tot_nodes,
                             closeness = closeness(net_powerlaw_all1),
                             betweenness = betweenness(net_powerlaw_all1),
                             degree = degree(net_powerlaw_all1))

## Experiment

perturbation_node <- function(v, h, u, graph, PtE)
{
  u_h <- u
  u_h[v] <- u[v] + h # vertexes are the first nv observations
  data <- data.frame(edge_number = PtE[,1],
                     distance_on_edge =  PtE[,2],
                     u_h = u_h)
  graph_tmp <- graph$clone(deep = TRUE)
  graph_tmp$add_observations(data = data, 
                             normalized = TRUE)

  err_exp <- try(res_exp <- graph_lme(u_h ~ -1, 
                                      graph = graph_tmp, 
                                      model = list(type = "isoCov")))
  err_WM <- try(res_WM <- graph_lme(u_h ~ -1, 
                                    graph = graph_tmp, 
                                    model = 'WM'))
  
  if(inherits(err_exp, "try-error")){
    coeff_exp <- rep(NA,2)
  } else coeff_exp <- res_exp$coeff$random_effects
  if(inherits(err_WM, "try-error")){
    coeff_WM <- rep(NA,3)
  } else coeff_WM <- res_WM$coeff$random_effects

  res <- c(v,
           h,
           coeff_exp, 
           coeff_WM)
  names(res) <- c('node',
                  'h',
                  'isoCov_tau',
                  'isoCov_kappa',
                  'WM_alpha',
                  'WM_tau',
                  'WM_kappa')
  return(res)
}

single_simulation <- function(iter, graph, n_obs, tau, kappa, alpha)
{
  PtE <- PtE_obs_nodes_edge(n_obs, graph$E, graph$edge_lengths)
  h_seq <- seq(-1.6, 1.6, by = 0.2) #seq(-2, 2, by = 0.2)
  nv <- nrow(graph$V)
  v_h_grid <- expand.grid(v = 1:nv, h = h_seq)
  u <- sample_spde(tau = tau,
                   kappa = kappa,
                   graph = graph,
                   alpha = alpha,
                   sigma_e = 0,
                   PtE = PtE)
  res <- mapply(v = v_h_grid[,1],
                h = v_h_grid[,2],
                FUN = perturbation_node,
                MoreArgs = list(u = u,
                                graph = graph,
                                PtE = PtE),
                SIMPLIFY = FALSE,
                USE.NAMES = TRUE)
  res2 <- map_dfr(res, .f = ~ .) %>%
    mutate(iter = iter)
  return(res2)
}



run_experiment <- function(graph, n_obs, tau, kappa, alpha, nsim)
{
  
  res <- mcmapply(iter = 1:nsim,
                  FUN = single_simulation,
                  MoreArgs = list(graph = graph,
                                  n_obs = n_obs,
                                  tau = tau,
                                  kappa = kappa,
                                  alpha = alpha),
                  SIMPLIFY = FALSE,
                  USE.NAMES = TRUE)
  print(res)
  res2 <- map_dfr(res, .f = ~ .)
  print('Graph done.')
  return(res2)
}



graph_list <- list(#comp_graph, 
                   comp_graph_all1, 
                   star_graph, 
                   rand_graph_all1,
                   powerlaw_graph_all1)
names(graph_list) <- c(#'Complete', 
                       'Complete_all1', 
                       'Star', 
                       'Random_all1',
                       'Powerlaw_all1')

# One observation per network node, plus more observations on the edges
# Detail: on the edges, roughly one observation each unit length
n_obs_list <- tot_nodes + 10*sapply(graph_list,
                                    FUN = \(x){ceiling(sum(x$edge_lengths))})
names(n_obs_list) <- names(graph_list)

tau_true <- 1
kappa_true <- 1
alpha <- 1
nsim <- 40

results <- mapply(graph = graph_list,
                  n_obs = n_obs_list,
                  FUN = run_experiment,
                  MoreArgs = list(tau = tau_true,
                                  kappa = kappa_true,
                                  alpha = alpha,
                                  nsim = nsim),
                  SIMPLIFY = FALSE,
                  USE.NAMES = TRUE)

save.image(here('sensitivity',
                paste('sensitivity_',
                      'obsperunit_',
                      'avgloc_',
                      # 'LOCAL_',
                      format(Sys.Date(),'%d%m%y'),
                      '.Rdata', sep = '')))


