library(MetricGraph)

set.seed(35242)

## Complete Graph (10 total nodes)
tot_nodes <- 10
V_comp <- matrix(0, tot_nodes, 2)
for(i in 1:tot_nodes){
  V_comp[i,] <- c(cos(2*pi/tot_nodes*(i-1)), sin(2*pi/tot_nodes*(i-1)))
}
E_comp <- t(combn(tot_nodes,2))
comp_graph <- metric_graph$new(V = V_comp, E = E_comp)

## Simulate random locations
nobs_edge <- 10

PtE <- NULL
for (i in 1:nrow(comp_graph$E)) {
  PtE <- rbind(PtE, cbind(rep(i, nobs_edge), runif(nobs_edge)))
}

PtE_shuffle <- PtE[sample(nrow(PtE), size = nrow(PtE), replace = FALSE),]

tau <- 1
kappa <- 1
alpha <- 1

u <- sample_spde(tau = tau,
                 kappa = kappa,
                 graph = comp_graph,
                 alpha = alpha,
                 sigma_e = 0,
                 PtE = PtE)
data <- data.frame(edge_number = PtE[,1],
                   distance_on_edge =  PtE[,2],
                   u = u)

u_shuffle <- sample_spde(tau = tau,
                         kappa = kappa,
                         graph = comp_graph,
                         alpha = alpha,
                         sigma_e = 0,
                         PtE = PtE_shuffle)
data_shuffle <- data.frame(edge_number = PtE_shuffle[,1],
                           distance_on_edge =  PtE_shuffle[,2],
                           u = u_shuffle)

## Fit with ordered PtE
comp_graph$clear_observations()
comp_graph$add_observations(data = data, 
                            normalized = TRUE)

res_exp <- graph_lme(u ~ -1, 
                     graph = comp_graph, 
                     model = list(type = "isoCov"))
res_WM <- graph_lme(u ~ -1, 
                    graph = comp_graph, 
                    model = 'WM1')

## Fit with shuffled PtE
comp_graph$clear_observations()
comp_graph$add_observations(data = data_shuffle, 
                            normalized = TRUE)

res_exp_shuffle <- graph_lme(u ~ -1, 
                             graph = comp_graph, 
                             model = list(type = "isoCov"))
res_WM_shuffle <- graph_lme(u ~ -1, 
                            graph = comp_graph, 
                            model = 'WM1')

## Comparison
print(res_exp$coeff$random_effects)
print(res_exp_shuffle$coeff$random_effects)

print(res_WM$coeff$random_effects)
print(res_WM_shuffle$coeff$random_effects)

