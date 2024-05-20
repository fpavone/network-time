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

PtE_obs_edge <- function(nobs_edge, E)
{
  PtE <- NULL
  for(i in 1:nrow(E)){
    #add locations sampled at random to each edge
    PtE <- rbind(PtE, cbind(rep(i, nobs_edge), runif(nobs_edge)))
  }
  return(PtE)
}

PtE_obs_tot <- function(nobs, E, edge_length)
{
  # probability of observation on an edge is proportional to the length of the edge
  PtE <- matrix(NA, nrow = nobs, ncol = 2)
  PtE[,1] <- sample(1:nrow(E), 
                    size = nobs, 
                    prob = edge_length,
                    replace = TRUE)
  PtE[,2] <- runif(nobs)
  return(PtE)
}

PtE_obs_nodes_edge <- function(nobs, E, edge_length) # it generates observations on vertexes and edges
{
  # Total number of observations nobs = obs_on_nodes + obs_on_edges
  nv <- max(E)
  if(nobs <= nv) stop("too few observations")
  PtE <- rbind(PtE_allnodes(nv, E),
               PtE_obs_tot(nobs-nv, E, edge_length))
  return(PtE)
}


sample_graph_igraph <- function(n_nodes, degree_dist, max_it = 100)
{
  it <- 0
  check_simple <- check_connected <- FALSE
  while((it < max_it) & !(check_simple & check_connected)){
    degs <- sample(1:n_nodes, 
                   n_nodes, 
                   replace = TRUE, 
                   prob = degree_dist(1:n_nodes))
    if (sum(degs) %% 2 != 0) {
      degs[1] <- degs[1] + 1
    }
    g_ <- realize_degseq(degs, allowed.edge.types = "multi")
    g <- simplify(g_, remove.multiple = TRUE, remove.loops = TRUE)
    check_simple <- is_simple(g)
    check_connected <- is_connected(g)
    it <- it + 1
  }
  return(g)
}

igraph2metric <- function(g)
{
  V <- layout_nicely(g)
  E <- as_edgelist(g)
  graph_g <- metric_graph$new(V = V,
                              E = E)
  return(graph_g)
}


#################
## Porcu et al ##
#################

# gc_gc_cov <- nimbleFunction(     
#   run = function(dists = double(2),time = double(2), pars = double(1),
#                  sigma2 = double(0),tau2 = double(0)) {
#     returnType(double(2))
#     n <- dim(dists)[1]
#     result <- matrix(nrow = n, ncol = n, init = FALSE)
#     
#     for(i in 1:(n-1)){
#       for(j in (i+1):n){
#         temp1 <- (1 + (time[i,j]/pars[4])^(2*pars[1]))
#         temp <- sigma2/temp1^(pars[5] + 1) * (1 + (dists[i,j]/pars[3])/temp1)^(-pars[2])
#         
#         result[i, j] <- temp
#         result[j, i] <- temp
#       }
#     }
#     for(i in 1:(n)){
#       
#       temp1 <- (1 + (time[i,i]/pars[4])^(2*pars[1]))
#       temp <- sigma2/temp1^(pars[5] + 1) * (1 + (dists[i,i]/pars[3])/temp1)^(-pars[2])
#       
#       result[i, i] <- temp + tau2
#     }
#     return(result)
#   })
# 
# 
# fn_optim_gc <- function(par, y, dist_mat, time_mat, n){
#   cov_mat <- gc_gc_cov(dists = dist_mat,time = time_mat, 
#                        pars = c(1,2,exp(par[2]),exp(par[3]),1),
#                        sigma2 = exp(par[1]),
#                        tau2 = exp(par[4]))
#   dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
# }
