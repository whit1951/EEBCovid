library(patchwork)
library(magrittr)
library(tidygraph)
library(ggraph)
library(tidyverse)
library(parallel)

source("../../generate_EEB_networks.R")

mycols <- c("#2e4052", "#d95d39", "#754042")
theme_set(theme_bw())
time_max <- 200

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

EEB_nets <- generate_EEB_networks("../../")
g<-as_tbl_graph(EEB_nets$office)
h<-as_tbl_graph(EEB_nets$lab)
full_graph <- graph_join(g, h, by="name") %>% to_undirected()

network_structure <- list(measure_network_structure(h) %>% mutate(network = str_c("Just Shared Lab Space (", n(), " components)")),
                          measure_network_structure(full_graph) %>% mutate(network = str_c("Combined Lab and Office (", n(), " components)"))) %>%
  bind_rows() %>%
  pivot_longer(c(-network, -component), names_to="metric", values_to="value") %>%
  mutate(metric = fct_inorder(metric))

set.seed(0)
full_net <- ggraph(full_graph, layout=igraph::layout_nicely(full_graph)) +
  geom_edge_link(edge_width=0.66, colour="#635E5B") +
  geom_node_point(aes(colour=I(mycols[1])), size=2) +
  ggtitle("A) Combined Lab and Office",
          str_c(with_graph(full_graph, graph_component_count()), " distinct components")) +
  theme_graph(base_family="") +
  theme(legend.position="none",
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))

set.seed(0)
lab_net <- ggraph(h, layout=igraph::layout_nicely(h)) +
  geom_edge_link(edge_width=0.66, colour="#635E5B") +
  geom_node_point(aes(colour=I(mycols[2])), size=2) +
  ggtitle("B) Just Shared Labspace",
          str_c(with_graph(h, graph_component_count()), " distinct components")) +
  theme_graph() +
  theme(legend.position="none",
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))

path_lengths <- ggplot(bind_rows(tibble(distance=get_distances(h), network="lab"),
                                 tibble(distance=get_distances(full_graph), network="full"))) +
  aes(x=distance, fill=network) +
  geom_histogram(position="dodge", binwidth=0.5) +
  scale_x_continuous(breaks=1:8) +
  scale_fill_manual(values=mycols) +
  xlab("Shortest Path Length") + ylab("Number of Paths") +
  ggtitle("C) Distribution of shortest paths") +
  theme(legend.position="none",
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))

full_plot <- full_net + lab_net + path_lengths + plot_layout(widths=c(1,1,1.5))
ggsave(full_plot, filename="network_structure.tiff", width=13, height=5)
ggsave(full_plot, filename="network_structure.png", width=13, height=5)

component_structure <- ggplot(network_structure %>% na.omit()) +
  aes(x=metric, y=value, colour=network) +
  geom_point(position=position_jitterdodge(jitter.width=0.1, dodge.width=0.5)) +
  coord_trans(y="log1p") +
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200),
                     minor_breaks=c(6, 7, 8, 9, 60, 70, 80, 90),
                     limits=c(NA, 160)) +
  scale_colour_manual(values=mycols)+ scale_fill_manual(values=mycols) +
  # ggtitle("Component-wise structural metrics") +
  theme(axis.title=element_blank(),
        legend.position="none",
        axis.ticks.x=element_blank(),
        plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                  vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
        plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(b = 5.5)))
ggsave(component_structure, filename="component_structure.tiff", width=5, height=4)
ggsave(component_structure, filename="component_structure.png", width=5, height=4)
