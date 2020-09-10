source("additional_simulations_for_revisions.R")
source("../generate_EEB_networks.R")

library(janitor)

mycols <- c("#2e4052", "#d95d39", "#754042", "#157a6e")
theme_set(theme_bw())
time_max <- 200
num_nets <- 200
num_sims <- 50

EEB_nets <- generate_EEB_networks()

make_islands_graph <- function(subnet_dist) {
  lapply(1:sum(subnet_dist$n), function(n) {
    play_erdos_renyi(sample(subnet_dist$size,
                            subnet_dist$n %>% sum(),
                            prob=subnet_dist$percent,
                            replace=TRUE),
                     1, directed=FALSE)
  }) %>%
    bind_graphs()
}

empirical_subnet_dist <- components(EEB_nets$lab)$csize %>%
  sort() %>% enframe(value="size") %>% tabyl(size)

# N.b. this does not create the few cases of two labs joined
#      by one individual that we see in the empirical data
networks <- mclapply(1:num_nets, mc.cores=7, function(x) make_islands_graph(empirical_subnet_dist))

network_structure <- mclapply(networks, mc.cores=7, measure_network_structure) %>%
  bind_rows() %>%
  pivot_longer(-component, names_to="metric", values_to="value") %>%
  mutate(metric = fct_inorder(metric))

network_structure_plot <- ggplot(network_structure %>% na.omit()) +
  aes(x=metric, y=value) +
  geom_jitter(width=0.25, height=0) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200),
                     minor_breaks=c(6, 7, 8, 9, 60, 70, 80, 90),
                     limits=c(NA, 160)) +
  scale_colour_manual(values=mycols)+ scale_fill_manual(values=mycols) +
  ggtitle("Component-wise structural metrics",
          str_glue("for {num_nets} rewired block model randomizations of the shared lab-space network")) +
  theme(axis.title=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(b = 5.5)))

# # All shortest path lengths will be 1 in this model
# distances <- lapply(1:length(networks), function(ii) {
#   networks[[ii]] %>% get_distances() %>% enframe(name=NULL) %>% mutate(network=ii)
# }) %>% bind_rows() %>%
#   count(network, value)
#
# distances_plot <- ggplot(distances) +
#   aes(x=as.factor(value), y=n) +
#   geom_jitter(width=0.25, height=0) +
#   scale_x_discrete(breaks=1:8) +
#   xlab("Shortest Path Length") + ylab("Number of Paths") +
#   ggtitle("Distribution of shortest paths",
#           str_glue("for {num_nets} rewired block model randomizations of the shared lab-space network")) +
#   theme(plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
#                                   vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
#         plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))

disease_simulations <- lapply(networks, epidemic_summary, num_sims=num_sims, cores=7) %>%
  bind_rows() %>%
  pivot_longer(-rep, names_to="metric", values_to="value")

disease_simulation_plot <- ggplot(disease_simulations %>% filter(metric != "dynamics_unfinished")) +
  aes(x=metric, y=value) +
  geom_jitter(width=0.25, height=0) +
  scale_colour_manual(values=mycols) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200),
                     minor_breaks=c(6, 7, 8, 9, 60, 70, 80, 90),
                     limits=c(NA, 160)) +
  ggtitle("Epidemic Simulation Outcomes",
          str_glue("for {num_sims} simulations on each of {num_nets} randomizations of the shared lab-space network")) +
  theme(axis.title=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))

save.image("../../results/additional_sims_lab_block.RData")

# ggsave(network_structure_plot,  filename=str_glue("../../figures/configuration_model_simulations/network_structure_{str_remove(net, '\\\\s.*')}.png"), width=5, height=5)
# ggsave(distances_plot,          filename=str_glue("../../figures/configuration_model_simulations/distances_{str_remove(net, '\\\\s.*')}.png"),         width=5, height=5)
# ggsave(disease_simulation_plot, filename=str_glue("../../figures/configuration_model_simulations/disease_dynamics_{str_remove(net, '\\\\s.*')}.png"),  width=5, height=5)
