library(magrittr)
library(tidyverse)
library(tidygraph)
library(parallel)
library(igraph)

theme_set(theme_bw())

set.seed(0)

load("additional_sims_lab_block.RData")

lab_disease_simulations <- disease_simulations
lab_distances <- distances
lab_network_structure <- network_structure

EEB_nets <- lapply(EEB_nets, as_tbl_graph)

network_structure <- bind_rows(
  lab_network_structure %>% mutate(net = "Just Shared Lab Space",
                                   rand = "Block Model Randomization"),
  measure_network_structure(EEB_nets$lab) %>%
    mutate(net = "Just Shared Lab Space",
           rand = "Empirical Network") %>%
    bind_rows() %>%
    pivot_longer(c(-net, -rand, -component), names_to="metric", values_to="value") %>%
    mutate(metric = fct_inorder(metric))
)

network_structure_plot <- ggplot(network_structure %>% na.omit()) +
  aes(x=metric, y=value, shape=rand,
      group=str_c(net, rand) %>% factor(levels=c("Just Shared Lab SpaceBlock Model Randomization",
                                                 "Just Shared Lab SpaceEmpirical Network"))) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.75)) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20),
                     minor_breaks=c(6, 7, 8, 9),
                     limits=c(NA, 20)) +
  scale_colour_manual(values=mycols, aesthetics=c("colour", "fill")) +
  scale_shape_manual(values=c(4, 16)) +
  ggtitle("Component-wise structural metrics",
          str_glue("for {num_nets} block model randomizations of each empirical network")) +
  theme(axis.title=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(b = 5.5)),
        legend.box="vertical", legend.margin=margin(0, 0, 0, 0), legend.spacing.y=unit(0, "pt"))
ggsave(network_structure_plot,  filename="figures/block_model_simulations/network_structure_compare_rand.png", width=6, height=5)

disease_simulations <- bind_rows(
  lab_disease_simulations %>% mutate(net = "Just Shared Lab Space",
                                     rand = "Configuration Model Randomization"),
  epidemic_summary(EEB_nets$lab, cores=7) %>%
    mutate(net = "Just Shared Lab Space", rand = "Empirical Network") %>%
    pivot_longer(c(-rep, -net, -rand), names_to="metric", values_to="value")
)

disease_simulation_plot <- ggplot(disease_simulations %>% filter(metric != "dynamics_unfinished")) +
  aes(x=metric, y=value, shape=rand,
      group=str_c(net, rand) %>% factor(levels=c("Just Shared Lab SpaceBlock Model Randomization",
                                                 "Just Shared Lab SpaceEmpirical Network"))) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.75)) +
  scale_colour_manual(values=mycols) +
  scale_shape_manual(values=c(4, 16)) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20),
                     minor_breaks=c(6, 7, 8, 9),
                     limits=c(NA, 160)) +
  ggtitle("Epidemic Simulation Outcomes",
          str_glue("for {num_sims} simulations on each of {num_nets} randomizations of each empirical network")) +
  theme(axis.title=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)),
        legend.box="vertical", legend.margin=margin(0, 0, 0, 0), legend.spacing.y=unit(0, "pt"))
ggsave(disease_simulation_plot, filename="figures/block_model_simulations/disease_dynamics_compare_rand.png", width=6, height=5)
