source("additional_simulations_for_revisions.R")

library(akima)

mycols <- c("#2e4052", "#d95d39", "#754042", "#157a6e")
zscore <- function(vec) (vec - mean(vec, na.rm=TRUE)) / sd(vec, na.rm=TRUE)
theme_set(theme_bw())
time_max <- 200
num_nets <- 1024
num_sims <- 20

# size and connectance
minsize <- 100; maxsize <- 500
minconn <- 0.0005; maxconn <- 0.15
# sizes <- runif(num_nets, minsize, maxsize) %>% round()
# powers <- runif(num_nets, 0.1, 3)
# growths <- runif(num_nets, 1, 10) %>% round()
# connectances <- runif(num_nets, minconn, maxconn)

sizes <- rep(seq(minsize, maxsize, length.out=ceiling(sqrt(num_nets))), each=ceiling(sqrt(num_nets)))
connectances <- rep(seq(minconn, maxconn, length.out=ceiling(sqrt(num_nets))), ceiling(sqrt(num_nets)))

networks <- lapply(1:num_nets, function(ii) play_erdos_renyi(sizes[ii], connectances[ii], directed=FALSE))
# networks <- lapply(1:num_nets, function(ii) play_barabasi_albert(sizes[ii], power=powers[ii], growth=growths[ii], directed=FALSE))

# connectances <- lapply(networks, . %>%  with_graph(graph_size() / graph_order() / (graph_order() - 1))) %>% unlist()

# network_structure <- lapply(networks, measure_network_structure) %>%
#   tibble(data=.) %>%
#   mutate(net_size=sizes, connectance=connectances) %>%
#   unnest(data) %>%
#   pivot_longer(c(-component, -net_size, -connectance), names_to="metric", values_to="value") %>%
#   mutate(metric = fct_inorder(metric))
#
# network_structure_plot <- network_structure %>%
#   na.omit() %>%
#   group_by(metric) %>%
#   group_modify(function(data, ...) {
#     a <- interp(x=data$net_size, y=data$connectance, z=zscore(data$value),
#                 xo=seq(minsize, maxsize, length.out=400),
#                 yo=seq(minconn, maxconn, length.out=400),
#                 duplicate="mean")
#     res <- a$z %>%
#       set_colnames(a$y) %>%
#       as_tibble() %>%
#       mutate(net_size=a$x) %>%
#       pivot_longer(-net_size, names_to="connectance", values_to="z-score",
#                    values_drop_na=TRUE,
#                    names_transform=list(connectance=as.numeric))
#     return(res)
#   }) %>%
#   {ggplot(.) +
#       aes(x=net_size, y=connectance, fill=`z-score`) +
#       geom_tile() +
#       facet_wrap(~metric) +
#       scale_fill_viridis_c(option="E")}

# distances <- lapply(1:length(networks), function(ii) {
#   networks[[ii]] %>% get_distances() %>% enframe(name=NULL) %>%
#     mutate(network=ii, size=sizes[ii], connectance=connectances[ii])
# }) %>%
#   bind_rows() %>%
#   group_by(network, size, connectance) %>%
#   group_modify(~tibble(metric=c("mean", "median", "maximum"),
#                        value=c(mean(.$value), median(.$value), max(.$value))))
#
# distances_plot <- distances %>%
#   group_by(metric) %>%
#   group_modify(function(data, ...) {
#     a <- interp(x=data$size, y=data$connectance, z=zscore(data$value),
#                 xo=seq(minsize, maxsize, length.out=400),
#                 yo=seq(minconn, maxconn, length.out=400),
#                 duplicate="mean")
#     res <- a$z %>%
#       set_colnames(a$y) %>%
#       as_tibble() %>%
#       mutate(size=a$x) %>%
#       pivot_longer(-size, names_to="connectance", values_to="z-score",
#                    values_drop_na=TRUE,
#                    names_transform=list(connectance=as.numeric))
#     return(res)
#   }) %>%
#   {ggplot(.) +
#       aes(x=size, y=connectance, fill=`z-score`) +
#       geom_tile() +
#       facet_wrap(~metric) +
#       scale_fill_viridis_c(option="E")}

disease_simulations <- mclapply(networks, mc.cores=7, epidemic_summary, num_sims=num_sims, cores=1) %>%
  bind_rows() %>%
  mutate(size=rep(sizes, each=num_sims), connectance=rep(connectances, each=num_sims)) %>%
  pivot_longer(c(-rep, -size, -connectance), names_to="metric", values_to="value")

disease_simulation_plot <- disease_simulations %>%
  filter(metric != "dynamics_unfinished") %>%
  group_by(metric) %>%
  mutate(zscore=zscore(value)) %>%
  group_by(metric, size, connectance) %>%
  summarise(zscore=mean(zscore)) %>%
  {ggplot(.) +
      aes(x=size, y=connectance) +
      geom_raster(aes(fill=zscore), interpolate=TRUE) +
      facet_wrap(~metric, ncol=1) +
      scale_x_continuous(expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0,0)) +
      scale_fill_viridis_c(option="E")}

TTP_variance_plot <- disease_simulations %>%
  filter(metric == "Time to Peak") %>%
  mutate(zscore=zscore(value)) %>%
  group_by(metric, size, connectance) %>%
  summarise(zscore=var(zscore)) %>%
  {ggplot(.) +
      aes(x=size, y=connectance) +
      geom_raster(aes(fill=zscore), interpolate=TRUE) +
      facet_wrap(~metric, ncol=1) +
      scale_x_continuous(expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0,0)) +
      scale_fill_viridis_c(option="E")}

save.image("../../results/additional_sims_ER_size_2.RData")
# save.image("../../results/additional_sims_BA_size.RData")

# ggsave(network_structure_plot,  filename="../../figures/increasing_size_simulations/network_structure.png", width=9, height=3)
# ggsave(distances_plot,          filename="../../figures/increasing_size_simulations/distances.png",         width=9, height=3)
ggsave(disease_simulation_plot, filename="../../figures/increasing_size_simulations/disease_dynamics.png",  width=6, height=8)
ggsave(TTP_variance_plot, filename="../../figures/increasing_size_simulations/TTP_variance.png",  width=6, height=3)
