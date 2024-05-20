require(here)
require(tidyverse)
require(ggh4x)
require(patchwork)
require(mvtnorm)
require(MetricGraph)
require(sp)
require(parallel)
require(igraph)
require(reshape2)
theme_set(
  theme_light() +
    theme(
      strip.background = element_rect(color = 'gray', fill = 'white'),
      strip.text.x = element_text(color = 'black'),
      strip.text.y = element_text(color = 'black')
    )
)

source(here('sensitivity','helper_fun.R'))
source(here('sensitivity','basic_network.R'))

args <- commandArgs(trailingOnly = TRUE)
dist <- args[1]
if(dist == 'power'){
  beta <- as.numeric(args[2])
}

set.seed(43523)
n_nodes <- 50

if(dist == 'constant'){
  degree_dist <- \(x){(x + 1 - x)}
} else if(dist == 'exp'){
  degree_dist <- \(x){exp(-x)}
} else if(dist == 'power'){
  degree_dist <- \(x){(x^{-beta})}
}

g <- sample_graph_igraph(n_nodes = n_nodes,
                         degree_dist = degree_dist)

graph <- igraph2metric(g)

graph_all1 <- graph$clone(deep = TRUE)
graph_all1$edge_lengths <- rep(1, length(graph$edge_lengths))
g_all1 <- graph_from_data_frame(d = tibble(source = graph_all1$E[,1],
                                           target = graph_all1$E[,2],
                                           weight = graph_all1$edge_lengths),
                                vertices = 1:nrow(graph_all1$V),
                                directed = FALSE)

graph_exp <- graph_all1
g_exp <- g_all1

graph_exp$compute_resdist(obs = FALSE)
graph_exp$compute_geodist(obs = FALSE)
dist_res <- as.matrix(graph_exp$res_dist$`__vertices`)
dist_geo <- as.matrix(graph_exp$geo_dist$`__vertices`)


H <- diag(n_nodes) - (1/n_nodes)*rep(1,n_nodes)%*%t(rep(1,n_nodes))
#only for unitary edges nets (therefore we don't do inverse on distance)
adj <- as.matrix(as_adj(g_exp))
# diag(adj) <- colSums(adj)
A_moran <- H %*% adj %*% H
B_moran <- H

# Setup for Moran's I with ball of res-dist = 1
adj_res <- 1/dist_res 
diag(adj_res) <- 0
for(i in 1:(n_nodes-1)){
  for(j in (i+1):n_nodes){
    if(adj_res[i,j] < 1){
      adj_res[i,j] <- adj_res[j,i] <- 0
    }
  }
}
A_moran_res <- H %*% adj_res %*% H
B_moran_res <- H



true_edgelengths <- graph_exp$edge_lengths
n_nodes <- nrow(graph_exp$V)
PtE_nodes <- PtE_allnodes(n_nodes, graph_exp$E)
graph_exp$edge_lengths <- 1/true_edgelengths
graph_exp$compute_laplacian(obs = FALSE)
L <- as.matrix(graph_exp$Laplacian$`__vertices`)

loc_quadratic <- function(z, Arow, Brow, Q, B)
{
  res <- as.numeric(2/(t(z)%*%B%*%z)*t(Arow - Q*Brow)%*%z)
  return(res)
}
asym_quadratic <- function(A, B)
{
  res <- diag(A)/diag(B)
  return(res)
}

sigma_seq <- c(0.5,1,2)
kappa_seq <- c(0.1,1,10)
nsim <- 100

loc_inf_S <- array(NA, dim = c(nsim, 
                               n_nodes, 
                               length(sigma_seq),
                               length(kappa_seq)),
                   dimnames = list(simulation = 1:nsim,
                                   node = 1:n_nodes,
                                   sigma = sigma_seq,
                                   kappa = kappa_seq))
loc_inf_S_teo <- array(NA, dim = c(n_nodes, 
                                   length(sigma_seq),
                                   length(kappa_seq)),
                       dimnames = list(node = 1:n_nodes,
                                       sigma = sigma_seq,
                                       kappa = kappa_seq))
S <- array(NA, dim = c(nsim,
                       length(sigma_seq),
                       length(kappa_seq)),
           dimnames = list(simulation = 1:nsim,
                           sigma = sigma_seq,
                           kappa = kappa_seq))



loc_inf_moran <- array(NA, dim = c(nsim, 
                                   n_nodes, 
                                   length(sigma_seq),
                                   length(kappa_seq)),
                       dimnames = list(simulation = 1:nsim,
                                       node = 1:n_nodes,
                                       sigma = sigma_seq,
                                       kappa = kappa_seq))
loc_inf_moran_teo <- array(NA, dim = c(n_nodes, 
                                       length(sigma_seq),
                                       length(kappa_seq)),
                           dimnames = list(node = 1:n_nodes,
                                           sigma = sigma_seq,
                                           kappa = kappa_seq))
I_moran <- array(NA, dim = c(nsim,
                             length(sigma_seq),
                             length(kappa_seq)),
                 dimnames = list(simulation = 1:nsim,
                                 sigma = sigma_seq,
                                 kappa = kappa_seq))


loc_inf_moran_res <- array(NA, dim = c(nsim, 
                                       n_nodes, 
                                       length(sigma_seq),
                                       length(kappa_seq)),
                           dimnames = list(simulation = 1:nsim,
                                           node = 1:n_nodes,
                                           sigma = sigma_seq,
                                           kappa = kappa_seq))
loc_inf_moran_res_teo <- array(NA, dim = c(n_nodes, 
                                           length(sigma_seq),
                                           length(kappa_seq)),
                               dimnames = list(node = 1:n_nodes,
                                               sigma = sigma_seq,
                                               kappa = kappa_seq))
I_moran_res <- array(NA, dim = c(nsim,
                                 length(sigma_seq),
                                 length(kappa_seq)),
                     dimnames = list(simulation = 1:nsim,
                                     sigma = sigma_seq,
                                     kappa = kappa_seq))


loc_inf_raylegh <- array(NA, dim = c(nsim, 
                                     n_nodes, 
                                     length(sigma_seq),
                                     length(kappa_seq)),
                         dimnames = list(simulation = 1:nsim,
                                         node = 1:n_nodes,
                                         sigma = sigma_seq,
                                         kappa = kappa_seq))
loc_inf_raylegh_teo <- array(NA, dim = c(n_nodes, 
                                         length(sigma_seq),
                                         length(kappa_seq)),
                             dimnames = list(node = 1:n_nodes,
                                             sigma = sigma_seq,
                                             kappa = kappa_seq))
L_raylegh <- array(NA, dim = c(nsim,
                               length(sigma_seq),
                               length(kappa_seq)),
                   dimnames = list(simulation = 1:nsim,
                                   sigma = sigma_seq,
                                   kappa = kappa_seq))



dist_mat <- graph_exp$compute_resdist_PtE(PtE = PtE_nodes,
                                          normalized = TRUE,
                                          include_vertices = FALSE)

Sigma_array <- array(NA, dim = c(length(sigma_seq),
                                 length(kappa_seq),
                                 n_nodes,
                                 n_nodes),
                     dimnames = list(sigma = sigma_seq,
                                     kappa = kappa_seq,
                                     node = 1:n_nodes,
                                     node = 1:n_nodes))

for(sigma in sigma_seq){
  for(kappa in kappa_seq){
    # Sigma <- sigma^2*exp(-kappa*D)
    Sigma <- as.matrix(exp_covariance(dist_mat, c(sigma, kappa)))
    Sigma_array[as.character(sigma),
                as.character(kappa),
                ,] <- Sigma
    
    loc_inf_S_teo[,
                  as.character(sigma),
                  as.character(kappa)] <- 4*diag(L%*%Sigma%*%t(L))
    # for ratio of quadratic forms is not correct formular
    loc_inf_moran_teo[,
                      as.character(sigma),
                      as.character(kappa)] <- 4*diag(A_moran%*%Sigma%*%t(A_moran))
    loc_inf_moran_res_teo[,
                          as.character(sigma),
                          as.character(kappa)] <- 4*diag(A_moran_res%*%Sigma%*%t(A_moran_res))
    loc_inf_raylegh_teo[,
                        as.character(sigma),
                        as.character(kappa)] <- 4*diag(L%*%Sigma%*%t(L))
    
    for(i in 1:nsim){
      u_nodes <- t(chol(Matrix::forceSymmetric(Sigma)))%*%rnorm(n_nodes)
      
      S[as.character(i),
        as.character(sigma),
        as.character(kappa)] <- as.numeric(t(u_nodes) %*% L %*% u_nodes)
      loc_inf_S[as.character(i),
                ,
                as.character(sigma),
                as.character(kappa)] <- apply(L,
                                              MARGIN = 1,
                                              FUN = function(x, u_nodes){as.numeric(2*t(x) %*% u_nodes)},
                                              u_nodes = u_nodes,
                                              simplify = TRUE)
      
      I_moran[as.character(i),
              as.character(sigma),
              as.character(kappa)] <- as.numeric(t(u_nodes) %*% A_moran %*% u_nodes)/as.numeric(t(u_nodes) %*% B_moran %*% u_nodes)
      loc_inf_moran[as.character(i),
                    ,
                    as.character(sigma),
                    as.character(kappa)] <- mapply(Arow = as.list(data.frame(t(A_moran))),
                                                   Brow = as.list(data.frame(t(B_moran))),
                                                   FUN = loc_quadratic,
                                                   MoreArgs = list(z = u_nodes,
                                                                   Q = I_moran[i],
                                                                   B = B_moran),
                                                   USE.NAMES = TRUE)
      
      I_moran_res[as.character(i),
                  as.character(sigma),
                  as.character(kappa)] <- as.numeric(t(u_nodes) %*% A_moran_res %*% u_nodes)/as.numeric(t(u_nodes) %*% B_moran_res %*% u_nodes)
      loc_inf_moran_res[as.character(i),
                        ,
                        as.character(sigma),
                        as.character(kappa)] <- mapply(Arow = as.list(data.frame(t(A_moran_res))),
                                                       Brow = as.list(data.frame(t(B_moran_res))),
                                                       FUN = loc_quadratic,
                                                       MoreArgs = list(z = u_nodes,
                                                                       Q = I_moran_res[i],
                                                                       B = B_moran_res),
                                                       USE.NAMES = TRUE)
      
      
      L_raylegh[as.character(i),
                as.character(sigma),
                as.character(kappa)] <- as.numeric(t(u_nodes) %*% L %*% u_nodes)/as.numeric(t(u_nodes) %*% u_nodes)
      loc_inf_raylegh[as.character(i),
                      ,
                      as.character(sigma),
                      as.character(kappa)] <- mapply(Arow = as.list(data.frame(t(L))),
                                                     Brow = as.list(data.frame(t(diag(n_nodes)))),
                                                     FUN = loc_quadratic,
                                                     MoreArgs = list(z = u_nodes,
                                                                     Q = L_raylegh[i],
                                                                     B = diag(n_nodes)),
                                                     USE.NAMES = TRUE)
      
    }
  }
}

asym_moran <- asym_quadratic(A_moran,B_moran)
names(asym_moran) <- 1:n_nodes
asym_moran_res <- asym_quadratic(A_moran_res,B_moran_res)
names(asym_moran_res) <- 1:n_nodes
asym_raylegh <- asym_quadratic(L,diag(n_nodes))
names(asym_raylegh) <- 1:n_nodes

score_sens <- diag(as.matrix(adj %*% exp(-dist_mat) %*% t(adj)))

meta_data <- tibble(node = names(degree(g_exp)),
                    Degree_res = apply(dist_res,
                                       MARGIN = 1,
                                       FUN = \(x)(sum(x<=1)-1)),
                    Degree = degree(g_exp),
                    Betweenness = betweenness(g_exp, directed = FALSE),
                    Closeness = closeness(g_exp),
                    Harmonic = harmonic_centrality(g_exp),
                    Authority = authority_score(g_exp)$vector,
                    PageRank = page_rank(g_exp, directed = FALSE)$vector,
                    Triangles = sapply(1:n_nodes, FUN = count_triangles, graph = g_exp),
                    Eigen_centrality = eigen_centrality(g_exp, directed = FALSE)$vector,
                    Subgraph_centrality = subgraph_centrality(g_exp, diag = FALSE),
                    Sensitivity_score_adj = score_sens)


as_tibble(melt(loc_inf_moran)) %>%
  mutate(what = "local Moran's I") %>%
  bind_rows(melt(loc_inf_moran_res) %>%
              mutate(what = "local Moran's I (resistance)")) %>%
  bind_rows(melt(loc_inf_S) %>%
              mutate(what = "local S (Laplacian)")) %>%
  bind_rows(melt(loc_inf_raylegh) %>%
              mutate(what = "local Raylegh Laplacian")) %>%
  mutate_at(vars(node), as.character) -> data_plot

as_tibble(melt(loc_inf_moran_teo)) %>%
  mutate(what = "local Moran's I") %>%
  bind_rows(melt(loc_inf_moran_res_teo) %>%
              mutate(what = "local Moran's I (resistance)")) %>%
  bind_rows(melt(loc_inf_S_teo) %>%
              mutate(what = "local S (Laplacian)")) %>%
  bind_rows(melt(loc_inf_raylegh_teo) %>%
              mutate(what = "local Raylegh Laplacian")) %>%
  mutate_at(vars(node), as.character) -> data_plot_teo


# data_plot %>%
#   group_by(node, sigma, kappa, what) %>%
#   summarize(est_s2 = var(value)) %>%
#   full_join(data_plot_teo, by = c("node", "sigma", "kappa", "what")) %>%
#   filter(sigma == 1, kappa == 1,
#          value < 20) %>%
#   ggplot(aes(x = value, y = est_s2)) +
#   geom_point(aes(color = )) +
#   facet_grid(rows = vars(what), scale = 'free')

as_tibble(melt(L_raylegh)) %>%
  mutate(what = 'local Raylegh Laplacian') %>%
  bind_rows(as_tibble(melt(I_moran)) %>%
              mutate(what = "local Moran's I")) %>%
  bind_rows(as_tibble(melt(I_moran_res)) %>%
              mutate(what = "local Moran's I (resistance)")) %>%
  group_by(sigma, kappa, what) %>%
  summarize(mean = mean(value),
            median = median(value),
            sd = sd(value)) -> data_Q

tibble(value = c(asym_moran, asym_moran_res, asym_raylegh),
       node = rep(1:n_nodes, 3),
       what = c(rep("asymptotic Moran's I", n_nodes),
                rep("asymptotic Moran's I (resistance)", n_nodes),
                rep("asymptotic Raylegh Laplacian", n_nodes))) -> data_asymptotic

save.image(file = here('sensitivity', 
                       'shiny', 
                       paste('sample_sensitivity_ISOCOV_',
                             dist,
                             ifelse(dist=='power',beta,''),
                             '.Rdata', sep = '')))



