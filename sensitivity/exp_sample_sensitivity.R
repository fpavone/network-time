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
graph$plot() -> graph_plot
# ggsave(graph_plot,
#        filename = here('sensitivity','graph_power1.png'),
#        width = 8,
#        height = 7)

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

# as_tibble(dist_res) %>% 
#   mutate(rows = paste('V', row_number(), sep = '')) %>%
#   pivot_longer(-rows, names_to = 'cols') %>%
#   mutate(what = "Resistant") %>%
#   bind_rows(as_tibble(dist_geo) %>%
#               mutate(rows = paste('V', row_number(), sep = '')) %>%
#               pivot_longer(-rows, names_to = 'cols') %>%
#               mutate(what = "Geodesic")) %>%
#   ggplot(aes(x = factor(rows, levels = paste('V',1:n_nodes,sep='')),
#              y = factor(cols, levels = paste('V',1:n_nodes,sep='')))) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = 'black',
#                       high = 'white') + 
#   scale_y_discrete(limits = rev(levels(factor(paste('V',1:n_nodes,sep=''),
#                                               levels = paste('V',1:n_nodes,sep=''))))) +
#   facet_grid(cols = vars(what)) +
#   theme(axis.title = element_blank()) 

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

# colnames(adj_res) <- colnames(adj)
# as_tibble(adj) %>%
#   mutate(rows = as.character(row_number())) %>%
#   pivot_longer(-rows, names_to = 'cols') %>%
#   mutate(what = 'Adj') %>%
#   bind_rows(as_tibble(adj_res) %>%
#               mutate(rows = as.character(row_number())) %>%
#               pivot_longer(-rows, names_to = 'cols') %>%
#               mutate(what = 'Adj Res')) %>%
#   ggplot(aes(x = factor(rows, levels = paste(1:n_nodes,sep='')),
#              y = factor(cols, levels = paste(1:n_nodes,sep='')))) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_gradient(low = 'white',
#                       high = 'black') + 
#   scale_y_discrete(limits = rev(levels(factor(paste(1:n_nodes,sep=''),
#                                               levels = paste(1:n_nodes,sep=''))))) +
#   facet_grid(cols = vars(what)) +
#   theme(axis.title = element_blank()) 


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

tau_seq <- c(0.5,1,2)
kappa_seq <- c(0.1,1,10)
alpha_seq <- c(1,2)
nsim <- 50

loc_inf_S <- array(NA, dim = c(nsim, 
                               n_nodes, 
                               length(tau_seq),
                               length(kappa_seq),
                               length(alpha_seq)),
                   dimnames = list(simulation = 1:nsim,
                                   node = 1:n_nodes,
                                   tau = tau_seq,
                                   kappa = kappa_seq,
                                   alpha = alpha_seq))
S <- array(NA, dim = c(nsim,
                       length(tau_seq),
                       length(kappa_seq),
                       length(alpha_seq)),
           dimnames = list(simulation = 1:nsim,
                           tau = tau_seq,
                           kappa = kappa_seq,
                           alpha = alpha_seq))

loc_inf_moran <- array(NA, dim = c(nsim, 
                                   n_nodes, 
                                   length(tau_seq),
                                   length(kappa_seq),
                                   length(alpha_seq)),
                       dimnames = list(simulation = 1:nsim,
                                       node = 1:n_nodes,
                                       tau = tau_seq,
                                       kappa = kappa_seq,
                                       alpha = alpha_seq))
I_moran <- array(NA, dim = c(nsim,
                             length(tau_seq),
                             length(kappa_seq),
                             length(alpha_seq)),
                 dimnames = list(simulation = 1:nsim,
                                 tau = tau_seq,
                                 kappa = kappa_seq,
                                 alpha = alpha_seq))

loc_inf_moran_res <- array(NA, dim = c(nsim, 
                                       n_nodes, 
                                       length(tau_seq),
                                       length(kappa_seq),
                                       length(alpha_seq)),
                           dimnames = list(simulation = 1:nsim,
                                           node = 1:n_nodes,
                                           tau = tau_seq,
                                           kappa = kappa_seq,
                                           alpha = alpha_seq))
I_moran_res <- array(NA, dim = c(nsim,
                                 length(tau_seq),
                                 length(kappa_seq),
                                 length(alpha_seq)),
                     dimnames = list(simulation = 1:nsim,
                                     tau = tau_seq,
                                     kappa = kappa_seq,
                                     alpha = alpha_seq))

loc_inf_raylegh <- array(NA, dim = c(nsim, 
                                     n_nodes, 
                                     length(tau_seq),
                                     length(kappa_seq),
                                     length(alpha_seq)),
                         dimnames = list(simulation = 1:nsim,
                                         node = 1:n_nodes,
                                         tau = tau_seq,
                                         kappa = kappa_seq,
                                         alpha = alpha_seq))
L_raylegh <- array(NA, dim = c(nsim,
                               length(tau_seq),
                               length(kappa_seq),
                               length(alpha_seq)),
                   dimnames = list(simulation = 1:nsim,
                                   tau = tau_seq,
                                   kappa = kappa_seq,
                                   alpha = alpha_seq))

# loc_inf_raylegh <- matrix(NA, nrow = nsim, ncol = n_nodes)
# colnames(loc_inf_raylegh) <- 1:n_nodes
# L_raylegh <- vector("numeric", nsim)

# tau <- 2
# kappa <- 10
# alpha <- 2


for(tau in tau_seq){
  for(kappa in kappa_seq){
    for(alpha in alpha_seq){
      for(i in 1:nsim){
        u_nodes <- sample_spde(tau = tau,
                               kappa = kappa,
                               graph = graph_exp,
                               alpha = alpha,
                               sigma_e = 0,
                               PtE = PtE_nodes)
        
        S[as.character(i),
          as.character(tau),
          as.character(kappa),
          as.character(alpha)] <- as.numeric(t(u_nodes) %*% L %*% u_nodes)
        loc_inf_S[as.character(i),
                  ,
                  as.character(tau),
                  as.character(kappa),
                  as.character(alpha)] <- apply(L,
                               MARGIN = 1,
                               FUN = function(x, u_nodes){as.numeric(2*t(x) %*% u_nodes)},
                               u_nodes = u_nodes,
                               simplify = TRUE)
        
        I_moran[as.character(i),
                as.character(tau),
                as.character(kappa),
                as.character(alpha)] <- as.numeric(t(u_nodes) %*% A_moran %*% u_nodes)/as.numeric(t(u_nodes) %*% B_moran %*% u_nodes)
        loc_inf_moran[as.character(i),
                      ,
                      as.character(tau),
                      as.character(kappa),
                      as.character(alpha)] <- mapply(Arow = as.list(data.frame(t(A_moran))),
                                                     Brow = as.list(data.frame(t(B_moran))),
                                                     FUN = loc_quadratic,
                                                     MoreArgs = list(z = u_nodes,
                                                                     Q = I_moran[i],
                                                                     B = B_moran),
                                                     USE.NAMES = TRUE)
        
        I_moran_res[as.character(i),
                    as.character(tau),
                    as.character(kappa),
                    as.character(alpha)] <- as.numeric(t(u_nodes) %*% A_moran_res %*% u_nodes)/as.numeric(t(u_nodes) %*% B_moran_res %*% u_nodes)
        loc_inf_moran_res[as.character(i),
                          ,
                          as.character(tau),
                          as.character(kappa),
                          as.character(alpha)] <- mapply(Arow = as.list(data.frame(t(A_moran_res))),
                                                         Brow = as.list(data.frame(t(B_moran_res))),
                                                         FUN = loc_quadratic,
                                                         MoreArgs = list(z = u_nodes,
                                                                         Q = I_moran_res[i],
                                                                         B = B_moran_res),
                                                         USE.NAMES = TRUE)
        
        
        L_raylegh[as.character(i),
                  as.character(tau),
                  as.character(kappa),
                  as.character(alpha)] <- as.numeric(t(u_nodes) %*% L %*% u_nodes)/as.numeric(t(u_nodes) %*% u_nodes)
        loc_inf_raylegh[as.character(i),
                        ,
                        as.character(tau),
                        as.character(kappa),
                        as.character(alpha)] <- mapply(Arow = as.list(data.frame(t(L))),
                                                       Brow = as.list(data.frame(t(diag(n_nodes)))),
                                                       FUN = loc_quadratic,
                                                       MoreArgs = list(z = u_nodes,
                                                                       Q = L_raylegh[i],
                                                                       B = diag(n_nodes)),
                                                       USE.NAMES = TRUE)
        
      }
    }
  }
}

asym_moran <- asym_quadratic(A_moran,B_moran)
names(asym_moran) <- 1:n_nodes
asym_moran_res <- asym_quadratic(A_moran_res,B_moran_res)
names(asym_moran_res) <- 1:n_nodes
asym_raylegh <- asym_quadratic(L,diag(n_nodes))
names(asym_raylegh) <- 1:n_nodes

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
                    Subgraph_centrality = subgraph_centrality(g_exp, diag = FALSE))


as_tibble(melt(loc_inf_moran)) %>%
  mutate(what = "local Moran's I") %>%
  bind_rows(melt(loc_inf_moran_res) %>%
              mutate(what = "local Moran's I (resistance)")) %>%
  bind_rows(melt(loc_inf_S) %>%
              mutate(what = "local S (Laplacian)")) %>%
  bind_rows(melt(loc_inf_raylegh) %>%
              mutate(what = "local Raylegh Laplacian")) %>%
  mutate_at(vars(node), as.character) -> data_plot
  # full_join(meta_data, by = c("node"))


as_tibble(melt(L_raylegh)) %>%
  mutate(what = 'local Raylegh Laplacian') %>%
  bind_rows(as_tibble(melt(I_moran)) %>%
              mutate(what = "local Moran's I")) %>%
  bind_rows(as_tibble(melt(I_moran_res)) %>%
              mutate(what = "local Moran's I (resistance)")) %>%
  group_by(tau, kappa, alpha, what) %>%
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
                       paste('sample_sensitivity_',
                             dist,
                             ifelse(dist=='power',beta,''),
                             '.Rdata', sep = '')))

# data_plot %>%
#   select(simulation, node, tau, kappa, alpha, value, what) %>%
#   group_by(node, tau, kappa, alpha, what) %>%
#   summarise(loc_inf = median(value)) %>%
#   mutate_at(vars(node), as.numeric) %>%
#   arrange(node) %>%
#   filter(tau == 1, kappa == 1, alpha == 1)
#   
# 
# data_plot %>%
#   filter(tau == 1,
#          kappa == 1,
#          alpha == 1) %>%
#   ggplot(aes(x = factor(node, levels = 1:n_nodes), 
#              y = abs(value))) +
#   geom_boxplot(aes(fill = Degree)) +
#   scale_fill_viridis_c() +
#   facet_grid(rows = vars(what), scale = 'free') +
#   labs(x = 'Node') +
#   theme_light() +
#   theme(axis.title.y = element_blank())


if(FALSE){
  as_tibble(melt(loc_inf_moran)) %>%
    mutate(what = "local Moran's I") %>%
    bind_rows(melt(loc_inf_moran_res) %>%
                mutate(what = "local Moran's I (resistance)")) %>%
    bind_rows(melt(loc_inf_S) %>%
                mutate(what = "local S (Laplacian)")) %>%
    bind_rows(melt(loc_inf_raylegh) %>%
                mutate(what = "local Raylegh Laplacian")) %>%
    mutate_at(vars(node), as.character) %>%
    full_join(meta_data, by = c("node")) %>%
    mutate_at(vars(alpha),
              \(x){paste('alpha =', x)}) %>%
    mutate_at(vars(tau),
              \(x){paste('tau =', x)}) %>%
    mutate_at(vars(kappa),
              \(x){paste('kappa =', x)}) %>%
    # filter(what == "local Raylegh Laplacian") %>%
    ggplot(aes(x = factor(node, levels = 1:n_nodes), 
               y = abs(value))) +
    geom_boxplot(aes(fill = Triangles)) +
    facet_nested(cols = vars(tau,kappa),
                 rows = vars(what,alpha),
                 scale = 'free') +
    labs(x = 'Node') +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) -> plot_mega
  
  ggsave(plot_mega,
         filename = here('sensitivity', 'mega.pdf'),
         width = 30,
         height = 20)
}

#### SINGLE COMBO PLOT
# Local influence plot
# as_tibble(loc_inf_moran) %>%
#   pivot_longer(everything(), names_to = 'node') %>%
#   mutate(what = "local Moran's I") %>%
#   bind_rows(as_tibble(loc_inf_moran_res) %>%
#               pivot_longer(everything(), names_to = 'node') %>%
#               mutate(what = "local Moran's I (resistance)")) %>%
#   bind_rows(as_tibble(loc_inf_S) %>%
#               pivot_longer(everything(), names_to = 'node') %>%
#               mutate(what = "local S (Laplacian)")) %>%
#   bind_rows(as_tibble(loc_inf_raylegh) %>%
#               pivot_longer(everything(), names_to = 'node') %>%
#               mutate(what = "local Raylegh Laplacian")) %>%
#   full_join(meta_data, by = c("node")) %>%
#   filter(what == "local Raylegh Laplacian") %>%
#   ggplot(aes(x = factor(node, levels = 1:n_nodes), 
#              y = abs(value))) +
#   geom_boxplot(aes(fill = Triangles)) + 
#   facet_grid(rows = vars(what), scale = 'free') +
#   labs(x = 'Node') +
#   theme(axis.title.y = element_blank())
# 
# # Asymptotic influence plot
# L_minmax <- range(eigen(L)$values)
# Adj_minmax <- range(eigen(adj)$values)
# Adj_res_minmax <- range(eigen(adj_res)$values)
# 
# data_minmax <- tibble(min = c(L_minmax[1], Adj_minmax[1], Adj_res_minmax[1]),
#                       max = c(L_minmax[2], Adj_minmax[2], Adj_res_minmax[2]),
#                       what = c("Asymptotic Raylegh Laplacian",
#                                "Asymptotic Moran's I",
#                                "Asymptotic Moran's I (resistance)"))
# 
# data_Q <- tibble(what = c("Asymptotic Raylegh Laplacian",
#                           "Asymptotic Moran's I",
#                           "Asymptotic Moran's I (resistance)"),
#                  Q = c(mean(L_raylegh),
#                        mean(I_moran),
#                        mean(I_moran_res)))
# 
# 
# 
# as_tibble(asym_moran, rownames = 'node') %>%
#   mutate(what = "Asymptotic Moran's I") %>%
#   bind_rows(as_tibble(asym_moran_res, rownames = 'node') %>%
#               mutate(what = "Asymptotic Moran's I (resistance)")) %>%
#   bind_rows(as_tibble(asym_raylegh, rownames = 'node') %>%
#               mutate(what = "Asymptotic Raylegh Laplacian")) %>%
#   full_join(meta_data, by = "node") %>%
#   ggplot(aes(x = factor(node, levels = 1:n_nodes), y = value)) +
#   geom_col(aes(fill = Degree_res)) +
#   # geom_hline(data = data_minmax, aes(yintercept = min),
#   #            color = 'darkred') +
#   # geom_hline(data = data_minmax, aes(yintercept = max),
#   #            color = 'darkred') +
#   facet_grid(rows = vars(what), scale = 'free') +
#   labs(x = 'Node') +
#   theme(axis.title.y = element_blank())




