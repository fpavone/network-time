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