# The work is interesting and technically sound. I have only one request to
# corroborate their conclusions and to fully appreciate the proposed study. In
# finite networks randomness can drive the system to different equilibria. I
# suggest the authors to perform the following additional experiments in order
# to have a sensitivity analysis of dynamics on contact networks.
# [a] Run simulations over random networks with the same node statistics. Do
#     simulations converge to the same equilibria?
# [b] Perform the same experiments by generating networks with increasing size.

library(tidygraph)
library(tidyverse)
library(parallel)

source("simulate_disease_on_network.R")

measure_network_structure <- function(network) {
  as_tbl_graph(network) %>%
    morph(to_components) %>%
    crystallise() %>%
    rowwise() %>%
    mutate(`Size` = with_graph(graph, graph_order()),
           `Diameter` = with_graph(graph, graph_diameter(directed=FALSE)),
           `Mean Path Length` = with_graph(graph, graph_mean_dist(directed=FALSE))) %>%
    ungroup() %>%
    select(-graph) %>%
    rename(component = name)
}

get_distances <- . %>% distances() %>% .[upper.tri(.)] %>% as.vector() %>% .[is.finite(.)]

summarise_disease_output <- function(disease_output) {
  timeseries <- disease_output %>%
    select(-number_infectious_neighbors, -prob_change_status) %>%
    as_tibble() %>%
    pivot_longer(names_to="time", values_to="status", matches("\\d+")) %>%
    mutate(time=as.integer(time),
           status=factor(status, levels=c("S", "E", "I", "R"), labels=c("Susceptible", "Exposed", "Infectious", "Removed"))) %>%
    count(status, time) %>%
    complete(status, nesting(time), fill=list(n=0))

  tibble(`Epidemic Peak Size` = timeseries %>% filter(status == "Infectious") %>% use_series(n) %>% max(),
         `Final Epidemic Size` = timeseries %>% filter(status == "Removed", time == max(time)) %>% use_series(n),
         dynamics_unfinished = timeseries %>% filter(status == "Infectious", time == max(time)) %>% use_series(n) %>% sum(),
         `Time to Peak` = timeseries %>% filter(status == "Infectious", n == `Epidemic Peak Size`) %>% use_series(time) %>% head(1))
}

epidemic_summary <- function(network, num_sims=100, cores=1) {
  mclapply(1:num_sims, mc.cores=cores, function(xx) {
    run_disease_network_simulation(time_max, network, "SEIR", beta=0.175, sigma=1/3.69, gamma=1/9.5) %>%
      summarise_disease_output() %>%
      mutate(rep = xx)
  }) %>% bind_rows()
}
