library(magrittr)
library(tidygraph)
library(tidyverse)
library(parallel)

source("simulate_disease_on_network.R")
source("generate_EEB_networks.R")

set.seed(0)

mycols <- c("#b26e63",
            "#525b76",
            "#606d5d",
            "#766c7f",
            "#869d96")
theme_set(theme_bw())
time_max <- 200

measure_network_structure <- function(network) {
  as_tbl_graph(network) %>%
    morph(to_components) %>%
    crystallise() %>%
    rowwise() %>%
    mutate(order = with_graph(graph, graph_order()),
           diameter = with_graph(graph, graph_diameter(directed=FALSE)),
           mean_path_length = with_graph(graph, graph_mean_dist(directed=FALSE))) %>%
    ungroup() %>%
    select(-graph) %>%
    rename(component = name)
}

summarise_disease_output <- function(disease_output) {
  timeseries <- disease_output %>%
      select(-number_infectious_neighbors, -prob_change_status) %>%
      as_tibble() %>%
      pivot_longer(names_to="time", values_to="status", matches("\\d+")) %>%
      mutate(time=as.integer(time),
             status=factor(status, levels=c("S", "I", "R"), labels=c("Susceptible", "Infectious", "Removed"))) %>%
    count(status, time) %>%
    complete(status, nesting(time), fill=list(n=0))

  tibble(max_infected = timeseries %>% filter(status == "Infectious") %>% use_series(n) %>% max(),
         total_infected = timeseries %>% filter(status == "Removed", time == max(time)) %>% use_series(n),
         # dynamics_unfinished = timeseries %>% filter(status == "Infectious", time == max(time)) %>% use_series(n) %>% sum(),
         time_to_max_infected = timeseries %>% filter(status == "Infectious", n == max_infected) %>% use_series(time) %>% head(1))
}

epidemic_summary <- function(network, num_sims=100, cores=1) {
  mclapply(1:num_sims, mc.cores=cores, function(xx) {
    run_disease_network_simulation(time_max, network, "SIR", beta=0.5, gamma=0.25) %>%
      summarise_disease_output() %>%
      mutate(rep = xx)
  }) %>% bind_rows()
}

EEB_nets <- generate_EEB_networks()
g<-as_tbl_graph(EEB_nets$office)
h<-as_tbl_graph(EEB_nets$lab)
full_graph <- graph_join(g, h, by="name") %>% to_undirected()

network_structure <- list(measure_network_structure(g) %>% mutate(network = str_c("office (", n(), "components)")),
                          measure_network_structure(h) %>% mutate(network = str_c("lab (", n(), "components)")),
                          measure_network_structure(full_graph) %>% mutate(network = str_c("full (", n(), "components)"))) %>%
  bind_rows() %>%
  pivot_longer(c(-network, -component), names_to="metric", values_to="value")

ggplot(network_structure %>% na.omit()) +
  aes(x=metric, y=value, colour=network) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.75)) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200)) +
  scale_colour_manual(values=mycols) +
  theme(axis.title.x=element_blank())
ggsave(filename="figures/network_structure.pdf", width=7, height=5)

get_distances <- . %>% distances() %>% as.vector() %>% .[is.finite(.)]
bind_rows(tibble(distance=get_distances(g), network="office"),
          tibble(distance=get_distances(h), network="lab"),
          tibble(distance=get_distances(full_graph), network="full")) %>%
          {ggplot(.) +
              aes(x=distance, fill=network) +
              geom_histogram(position="dodge", binwidth=0.5) +
              scale_fill_manual(values=mycols) +
              xlab("Shortest Path Length")}
ggsave(filename="../figures/shortest_path_length_distribution.pdf", width=7, height=5)

disease_simulations <- bind_rows(epidemic_summary(g, cores=7) %>% mutate(network="office"),
          epidemic_summary(h, cores=7) %>% mutate(network="lab"),
          epidemic_summary(full_graph, cores=7) %>% mutate(network="full")) %>%
  pivot_longer(c(-rep, -network), names_to="metric", values_to="value")
ggplot(disease_simulations) +
  aes(x=metric, y=value, colour=network) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.75)) +
  scale_colour_manual(values=mycols) +
  theme(axis.title.x=element_blank())
ggsave(filename="../figures/disease_outcomes.pdf", width=7, height=5)
