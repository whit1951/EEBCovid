library(magrittr)
library(tidygraph)
library(tidyverse)
library(parallel)

source("../../simulate_disease_on_network.R")
source("../../generate_EEB_networks.R")

set.seed(0)


mycols <- c("#2e4052", "#d95d39")
theme_set(theme_bw())
time_max <- 150

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

EEB_nets <- generate_EEB_networks("../../")
g<-as_tbl_graph(EEB_nets$office)
h<-as_tbl_graph(EEB_nets$lab)
full_graph <- graph_join(g, h, by="name") %>% to_undirected() %>% simplify() %>% as_tbl_graph()

disease_simulations <- bind_rows(epidemic_summary(h, cores=7) %>% mutate(network="Just Shared Lab Space"),
                                 epidemic_summary(full_graph, cores=7) %>% mutate(network="Combined Lab and Office")) %>%
  pivot_longer(c(-rep, -network), names_to="metric", values_to="value")

disease_dynamics <- ggplot(disease_simulations %>% filter(metric != "dynamics_unfinished")) +
  aes(x=metric, y=value, colour=network) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.5)) +
  scale_colour_manual(values=mycols) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200),
                     minor_breaks=c(6, 7, 8, 9, 60, 70, 80, 90),
                     limits=c(NA, 160)) +
  theme(axis.title=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))
ggsave(disease_dynamics, filename="disease_dynamics.tiff", width=4.5, height=4.5)
ggsave(disease_dynamics, filename="disease_dynamics.png", width=4.5, height=4.5)
