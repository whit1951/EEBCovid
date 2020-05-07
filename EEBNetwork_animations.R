library(magrittr)
library(tidygraph)
library(tidyverse)
library(ggraph)
library(magick)
library(gganimate)
library(igraph)
library(Matrix)
library(tools)
library(patchwork)

source("simulate_disease_on_network.R")
source("generate_EEB_networks.R")

mycols <- c("Susceptible"="#2e294e", "Infectious"="#d00000")
time_max <- 150

combine_gifs <- function(plot1, plot2) {
  # read the plots and store them
  plot1 <- image_read(plot1)
  plot2 <- image_read(plot2)
  # sync the number of frames in each plot
  n1=length(plot1)
  n2=length(plot2)
  # match number of frames of the two plots
  if (!(n1 == n2)) plot1 <- rep(plot1, n2)
  if (!(n1 == n2)) plot2 <- rep(plot2, n1)
  # initialize the combined plot
  p <- image_append(c(plot1[1], plot2[1]))
  # grow the combined plot frame by frame
  n=ifelse(n1 == n2, n1, n1 * n2)
  n=min(1000, n)  # set max to 1000
  for (i in 2:(n-1)) {
    tmp <- image_append(c(plot1[i], plot2[i]))
    p <- c(p, tmp)
  }
  return(p)
}

make_composite_disease_gif <- function(network, title, filename) {
  set.seed(0)

  network <- as_tbl_graph(network)

  output <- run_disease_network_simulation(time_max, network, "SI", beta=0.025)

  tmp <- igraph::layout_nicely(output)

  tidy_output <- output %>%
    select(-number_infectious_neighbors, -prob_change_status) %>%
    as_tibble() %>%
    pivot_longer(names_to="time", values_to="status", matches("\\d+")) %>%
    mutate(time=as.integer(time),
           status=factor(status, levels=c("S", "I"), labels=c("Susceptible", "Infectious")))

  layout_matrix <- tmp[rep(1:nrow(tmp), times=time_max+1), ]

  net_plot <- ggraph(tbl_graph(nodes=tidy_output %>% arrange(time),
                               edges=output %>% activate(edges) %>% as_tibble()),
                     layout=layout_matrix) +
    geom_edge_link(edge_width=0.66) +
    geom_node_point(aes(colour=status), size=3) +
    scale_colour_manual(values=mycols) +
    transition_manual(time) +
    ggtitle(title) +
    theme_graph() +
    theme(legend.position="none")

  p <- ggraph(output %>% activate(nodes) %>% rename(status = time_max) %>%
                mutate(status = factor(status, levels=c("S", "I"),
                                       labels=c("Susceptible", "Infectious"))),
              layout=tmp) +
    geom_edge_link(edge_width=0.66, colour="#635E5B") +
    geom_node_point(aes(colour=status), size=3) +
    scale_colour_manual(values=mycols) +
    ggtitle(title) +
    theme_graph() +
    theme(legend.position="none")
  ggsave(p, filename=str_c(file_path_sans_ext(filename), "_final.tiff"), width=6, height=6)

  net_gif <- animate(net_plot, end_pause=30, width=450, height=450)

  timeseries <- tidy_output %>%
    count(status, time) %>%
    complete(status, nesting(time), fill=list(n=0))

  line_plot <- ggplot(timeseries) +
    aes(x=time, y=n, fill=status) +
    geom_area() +
    geom_segment(aes(xend=time_max, yend=n), linetype=2, colour="grey",
                 data=filter(timeseries, status == "Infectious")) +
    geom_label(aes(x=time_max + 0.25, label=status), nudge_y=-5, vjust="middle",
               hjust="left", size=6, label.padding=unit(0.5, "lines"), colour="white",
               data=filter(timeseries, status == "Infectious")) +
    geom_label(aes(x=time_max + 0.25, label="Suceptible"), nudge_y=5, vjust="middle",
               hjust="left", size=6, label.padding=unit(0.5, "lines"), colour="white",
               fill=mycols["Susceptible"],
               data=filter(timeseries, status == "Infectious")) +
    transition_reveal(time) +
    scale_fill_manual(values=mycols) +
    coord_cartesian(clip="off", expand=FALSE) +
    xlab("Time") + ylab("Number of Individuals") +
    theme_minimal() +
    theme(legend.position="none") +
    theme(plot.margin=margin(5.5, 100, 5.5, 5.5),
          axis.text.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.text.y=element_text(size=14),
          axis.title=element_text(size=16))
  line_gif <- animate(line_plot, end_pause=30, width=600, height=450)

  new_gif <- combine_gifs(net_gif, line_gif)

  image_write(new_gif, filename)

  # p1 <- ggplot(timeseries) +
  #   aes(x=time, y=n, fill=status) +
  #   geom_area() +
  #   scale_fill_manual(values=mycols) +
  #   scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
  #   coord_cartesian(clip="off", expand=FALSE) +
  #   xlab("Time") + ylab("Number of Individuals") +
  #   theme_minimal() +
  #   theme(legend.position="bottom",
  #         legend.title=element_blank(),
  #         plot.margin=margin(5.5, 100, 5.5, 5.5),
  #         axis.text.x=element_blank(),
  #         panel.grid.major.y=element_blank(),
  #         panel.grid.minor.y=element_blank(),
  #         axis.text.y=element_text(size=14),
  #         axis.title=element_text(size=16))
  # p + p1
  # ggsave(filename=str_c(file_path_sans_ext(filename), "_final_combined.tiff"), width=14, height=7)
}

EEB_nets <- generate_EEB_networks()

g<-as_tbl_graph(EEB_nets$office)
h<-as_tbl_graph(EEB_nets$lab)

make_composite_disease_gif(g, "Office Network", "figures/Office_Network.gif")
make_composite_disease_gif(h, "Lab Network", "figures/Lab_Network.gif")

full_graph<-graph_join(g, h, by="name") %>% to_undirected()
make_composite_disease_gif(full_graph, "Combined Office and Lab Network", "figures/Full_Network.gif")
