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

n_nodes <- 10

g <- sample_graph_igraph(n_nodes = n_nodes,
                         degree_dist = \(x){x^{-1}})

powerlaw_graph <- igraph2metric(g)
powerlaw_graph_all1 <- powerlaw_graph$clone(deep = TRUE)
powerlaw_graph_all1$edge_lengths <- rep(1, length(powerlaw_graph$edge_lengths))

net_powerlaw_all1 <- graph_from_data_frame(d = tibble(source = powerlaw_graph_all1$E[,1],
                                                      target = powerlaw_graph_all1$E[,2],
                                                      weight = powerlaw_graph_all1$edge_lengths),
                                           vertices = 1:nrow(powerlaw_graph_all1$V),
                                           directed = FALSE)


graph_exp <- powerlaw_graph_all1
PtE_nodes <- PtE_allnodes(n_nodes, graph_exp$E)

ones <- rep(1,n_nodes)
C <- diag(n_nodes) - (1/n_nodes)*matrix(1, n_nodes, n_nodes)


dist_mat <- graph_exp$compute_resdist_PtE(PtE = PtE_nodes,
                                          normalized = TRUE,
                                          include_vertices = FALSE)
dist_mat <- as.matrix(dist_mat)

graph_exp$compute_laplacian(obs = FALSE)
L <- as.matrix(graph_exp$Laplacian$`__vertices`)

Linv <- ginv(L)

D2 <- -2*Linv + t(diag(Linv))%x%ones + t(ones)%x%diag(Linv)
dist_mat


### Sampling data

sigma <- 1
kappa <- 1
Sigma <- as.matrix(exp_covariance(dist_mat, c(sigma, kappa)))

m <- 200 # n_obs
X <- matrix(NA, nrow = m, ncol = n_nodes)
for(i in 1:m){
  X[i,] <- as.numeric(t(chol(Matrix::forceSymmetric(Sigma)))%*%rnorm(n_nodes))
}

X_mean <- colMeans(X)
X_c <- X - matrix(X_mean, m, n_nodes, byrow = TRUE)

Sigma_hat <- 1/(m-1)*t(X_c)%*%X_c
Corr_Sigma_hat <- diag(diag(Sigma_hat)^(-0.5)) %*% Sigma_hat %*% diag(diag(Sigma_hat)^(-0.5))

R_hat <- -(1/kappa)*log(Corr_Sigma_hat/sigma^2)
diag(R_hat) <- 0 # it's numeric 0, enforicing 0

R <- -(1/kappa)*log(Sigma/sigma^2)

Linv_hat <- -0.5*C%*% R_hat %*% C

L_hat <- ginv(-0.5*C%*% R_hat %*% C) # not a proper laplacian :(


# L = Q x Lambda x t(Q)
Q <- eigen(Linv)$vectors
Lambda <- diag(eigen(Linv)$values)
Lambda[n_nodes,n_nodes] <- 0
Y <- Q %*% sqrt(Lambda)


### TODO: try to start from true Sigma, and add perturbation
## Sigma_hat = Sigma + ones*eta, eta is the perturbation
## Do the inverse path to obtain L and check at each step some
## matrix norm. The idea is to check if at some point something
## explodes.




