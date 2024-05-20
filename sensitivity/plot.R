require(here)
require(tidyverse)
require(patchwork)
require(mvtnorm)
require(MetricGraph)
require(sp)
require(igraph)
theme_set(
  theme_light() +
    theme(
      strip.background = element_rect(color = 'gray', fill = 'white'),
      strip.text.x = element_text(color = 'black'),
      strip.text.y = element_text(color = 'black')
    )
)

################################
# Results of exp_sensitivity.R #
################################

load(here('sensitivity','sensitivity_obsperunit_270923.Rdata'))
lapply(X = PtE_list,
       FUN = function(x){
         data <- data.frame(edge_number = x[,1],
                            distance_on_edge =  x[,2],
                            m = rep(1,nrow(x)))
         return(data)
       }) -> data_graph

graph_list[[1]]$add_observations(data = data_graph[[1]], normalized = TRUE)
graph_list[[2]]$add_observations(data = data_graph[[2]], normalized = TRUE)
graph_list[[3]]$add_observations(data = data_graph[[3]], normalized = TRUE)
graph_list[[4]]$add_observations(data = data_graph[[4]], normalized = TRUE)

(graph_list[[1]]$plot(data = TRUE) + guides(color = "none") + labs(title = 'Complete')) |
  (graph_list[[2]]$plot(data = TRUE) + guides(color = "none") + labs(title = 'Complete (all 1)')) |
  (graph_list[[3]]$plot(data = TRUE) + guides(color = "none") + labs(title = 'Star')) |
  (graph_list[[4]]$plot(data = TRUE) + guides(color = "none") + labs(title = 'Random (all 1)')) -> graphs_plots

ggsave(graphs_plots,
       filename = here('sensitivity','graphs.jpg'),
       width = 15,
       height = 6)


results %>%
  map_dfr(.f = as_tibble,
          .id = 'graph') %>%
  pivot_longer(-c('graph','node','iter','h')) %>%
  mutate(model = str_extract(name,'^[:alpha:]*[:digit:]?'),
         parameter = str_extract(name,'[:alpha:]*$')) %>%
  select(-name) -> res_data

comp_meta %>%  
  mutate(graph = "Complete") %>%
  bind_rows(comp_all1_meta %>%
              mutate(graph = "Complete_all1")) %>%
  bind_rows(star_meta %>%
              mutate(graph = "Star")) %>%
  bind_rows(rand_all1_meta %>%
              mutate(graph = "Random_all1")) %>%
  mutate(closeness_trasf =  (closeness - median(closeness))/diff(range(closeness)),
         betweenness_trasf =  (betweenness - median(betweenness))/diff(range(betweenness)),
         degree_trasf = (degree - median(degree))/diff(range(degree)),
         .by = graph) %>%
  replace_na(replace = list(degree_trasf = 0,
                            betweenness_trasf = 0,
                            closeness_trasf = 0)) -> net_meta

model_name <- 'isoCov'
res_data %>%
  filter(model == model_name) %>%
  group_by(graph, node, h, model, parameter) %>%
  summarise(value = median(value)) %>%
  full_join(net_meta, by = c('node', 'graph')) %>%
  ggplot(aes(x = h, y = value, group = node)) +
  geom_line(aes(color = betweenness_trasf), alpha = 0.6) +
  scale_color_viridis_c() + 
  labs(title = model_name, color = 'Betweennes (rescaled)') +
  facet_grid(rows = vars(parameter), 
             cols = vars(graph),
             scale = 'free') -> plot_isocov

model_name <- 'WM1'
res_data %>%
  filter(model == model_name) %>%
  group_by(graph, node, h, model, parameter) %>%
  summarise(value = median(value)) %>%
  full_join(net_meta, by = c('node', 'graph')) %>%
  ggplot(aes(x = h, y = value, group = node)) +
  geom_line(aes(color = betweenness_trasf), alpha = 0.6) +
  scale_color_viridis_c() + 
  labs(title = model_name, color = 'Betweennes (rescaled)') +
  facet_grid(rows = vars(parameter), 
             cols = vars(graph),
             scale = 'free') -> plot_wm1

plot_isocov / plot_wm1 + plot_layout(guides = "collect") -> plot_isocov_wm1
ggsave(plot_isocov_wm1,
       filename = 'isocov_wm1.jpg',
       width = 15,
       height = 10)

### Averaged locations
rm(list=ls())
load(here('sensitivity','sensitivity_obsperunit_avgloc_061023.Rdata'))

results %>%
  map_dfr(.f = as_tibble,
          .id = 'graph') %>%
  pivot_longer(-c('graph','node','iter','h')) %>%
  mutate(model = str_extract(name,'^[:alpha:]*[:digit:]?'),
         parameter = str_extract(name,'[:alpha:]*$')) %>%
  select(-name) -> res_data

comp_meta %>%  
  mutate(graph = "Complete_all1") %>%
  # bind_rows(comp_all1_meta %>%
  #             mutate(graph = "Complete_all1")) %>%
  bind_rows(star_meta %>%
              mutate(graph = "Star")) %>%
  bind_rows(rand_all1_meta %>%
              mutate(graph = "Random_all1")) %>%
  bind_rows(powerlaw_all1_meta %>%
              mutate(graph = "Powerlaw_all1")) %>%
  mutate(closeness_trasf =  (closeness - median(closeness))/diff(range(closeness)),
         betweenness_trasf =  (betweenness - median(betweenness))/diff(range(betweenness)),
         degree_trasf = (degree - median(degree))/diff(range(degree)),
         .by = graph) %>%
  replace_na(replace = list(degree_trasf = 0,
                            betweenness_trasf = 0,
                            closeness_trasf = 0)) -> net_meta

model_name <- 'isoCov'
res_data %>%
  filter(model == model_name) %>%
  group_by(graph, node, h, model, parameter) %>%
  summarise(value = median(value)) %>%
  full_join(net_meta, by = c('node', 'graph')) %>%
  ggplot(aes(x = h, y = value, group = node)) +
  geom_line(aes(color = degree_trasf), alpha = 0.6) +
  scale_color_viridis_c() + 
  labs(title = model_name, color = 'Betweennes (rescaled)') +
  facet_grid(rows = vars(parameter), 
             cols = vars(graph),
             scale = 'free') #-> plot_isocov_avgloc

model_name <- 'WM'
res_data %>%
  filter(model == model_name) %>%
  group_by(graph, node, h, model, parameter) %>%
  summarise(value = median(value, na.rm = TRUE)) %>%
  full_join(net_meta, by = c('node', 'graph')) %>%
  # filter(value > 1.5, graph == 'Powerlaw_all1', parameter == 'alpha') %>% view()
  ggplot(aes(x = h, y = value, group = node)) +
  geom_line(aes(color = degree_trasf), alpha = 0.6) +
  scale_color_viridis_c() + 
  # geom_label(aes(label = node)) +
  labs(title = model_name, color = 'Betweennes (rescaled)') +
  facet_grid(rows = vars(parameter), 
             cols = vars(graph),
             scale = 'free') #-> plot_wm1_avgloc

plot_isocov_avgloc / plot_wm1_avgloc + plot_layout(guides = "collect") -> plot_isocov_wm1_avgloc
ggsave(plot_isocov_wm1_avgloc,
       filename = 'isocov_wm1_avgloc.jpg',
       width = 15,
       height = 10)

##OTHER

res_data %>%
  filter(node %in% c(4,7),
         graph == 'Complete_all1',
         value < 10) %>%
  ggplot(aes(x = factor(h, levels = seq(-2, 2, by = 0.2)),
             y = value)) +
  geom_boxplot() +
  facet_grid(rows = vars(parameter),
             cols = vars(model),
             scale = 'free')

PtE <- PtE_list[[1]]
data_tmp <- data.frame(edge_number = PtE[,1],
                       distance_on_edge =  PtE[,2],
                       u = rep(1,nrow(PtE)))
graph_list[["Complete"]]$add_observations(data = data_tmp, 
                                             normalized = TRUE)
graph_list[["Complete"]]$plot(data = "u")


comp_graph_all1$compute_laplacian(obs = FALSE)
L <- comp_graph_all1$Laplacian$`__vertices`

aa <- rep(0,10)
aa[1] <- 1
t(aa)%*%L%*%aa

################################
# Results of exp_stability.R   #
################################

load(here('sensitivity','stability_test_PtE_obs_nodes_edge_sorted270923.Rdata'))

results %>%
  map_dfr(.f = as_tibble,
          .id = 'graph') %>%
  pivot_longer(-graph) %>%
  mutate(model = str_extract(name,'^[:alpha:]*[:digit:]?'),
         parameter = str_extract(name,'[:alpha:]*$')) %>%
  select(-name) -> res_data

res_data %>%
  mutate(value_trunc = ifelse(value > 5, 5, value)) %>%
  ggplot(aes(x = graph, y = value_trunc)) +
  geom_boxplot(aes(fill = model), position = "dodge2") +
  facet_grid(rows = vars(parameter), scale = 'free')




  
  