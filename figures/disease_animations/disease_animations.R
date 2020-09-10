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

source("../../code/simulate_disease_on_network.R")
source("../../code/generate_EEB_networks.R")

mycols <- c("Susceptible" = "#2e4052",
            "Exposed"     = "#157a6e",
            "Infectious"  = "#d95d39",
            "Removed"     = "#754042")
time_max <- 75


combine_gifs <- function(plot1, plot2, vertical=FALSE) {
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
  p <- image_append(c(plot1[1], plot2[1]), stack=vertical)
  # grow the combined plot frame by frame
  n=ifelse(n1 == n2, n1, n1 * n2)
  n=min(1000, n)  # set max to 1000
  for (i in 2:(n-1)) {
    tmp <- image_append(c(plot1[i], plot2[i]), stack=vertical)
    p <- c(p, tmp)
  }
  return(p)
}

make_composite_disease_gif <- function(network, title, filename) {
  set.seed(0)

  network <- as_tbl_graph(network)
  tmp <- igraph::layout_nicely(network)

  output <- run_disease_network_simulation(time_max, network, "SEIR", beta=0.175, sigma=1/3.69, gamma=1/9.5)

  tidy_output <- output %>%
    select(-number_infectious_neighbors, -prob_change_status) %>%
    as_tibble() %>%
    pivot_longer(names_to="time", values_to="status", matches("\\d+")) %>%
    mutate(time=as.integer(time),
           status=factor(status, levels=c("S", "E", "I", "R"), labels=c("Susceptible", "Exposed", "Infectious", "Removed")))

  layout_matrix <- tmp[rep(1:nrow(tmp), times=time_max+1), ]

  net_plot <- ggraph(tbl_graph(nodes=tidy_output %>% arrange(time),
                               edges=output %E>% as_tibble()),
                     layout=layout_matrix) +
    geom_edge_link(edge_width=0.66) +
    geom_node_point(aes(colour=status), size=3.5) +
    scale_colour_manual(values=mycols) +
    transition_manual(time) +
    ggtitle(title) +
    theme_graph() +
    theme(legend.position="none",
          plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                    vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
          plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 5.5)))

  net_gif <- animate(net_plot, end_pause=30, width=450, height=450)

  timeseries <- tidy_output %>%
    count(status, time) %>%
    complete(status, nesting(time), fill=list(n=0))

  label_df <- tibble(time   = rep(0:time_max, 4),
                     x      = time_max + 0.25,
                     status = rep(c("Susceptible", "Exposed", "Infectious", "Removed"),
                                   each=time_max + 1) %>% fct_inorder()) %>%
    left_join(timeseries, by=c("time", "status")) %>%
    pivot_wider(names_from=status, values_from=n) %>%
    mutate(ynudge_Susceptible = 0,
           ynudge_Exposed = ifelse(Exposed < 5, 4, 0) + ifelse(Infectious < 5, 8, 0),
           ynudge_Infectious = ifelse(Infectious < 5, 4, 0),
           ynudge_Removed = ifelse(Removed < 5, -4, 0),
           Susceptible = Removed + Infectious + Exposed + Susceptible / 2,
           Exposed     = Removed + Infectious + Exposed / 2,
           Infectious  = Removed + Infectious / 2,
           Removed     = Removed / 2) %>%
           {full_join(
             select(., -starts_with("ynudge")) %>%
               pivot_longer(c(Susceptible, Exposed, Infectious, Removed),
                            names_to="status", values_to="n"),
             select(., -Susceptible, -Exposed, -Infectious, -Removed) %>%
               rename_at(vars(starts_with("ynudge")), ~str_remove(., "ynudge_")) %>%
               pivot_longer(c(Susceptible, Exposed, Infectious, Removed),
                            names_to="status", values_to="ynudge"))} %>%
    group_by(status) %>%
    mutate(xnudge = ifelse(n == max(n) & (status == "Exposed" | status == "Infectious"), 1000, 0))

  segment_df <- timeseries %>%
    pivot_wider(names_from=status, values_from=n) %>%
    mutate(Exposed    = ifelse(cumsum(Exposed) == sum(Exposed) & Exposed == 0,
                               n_max + 10, Removed + Infectious + Exposed),
           Infectious = ifelse(cumsum(Infectious) == sum(Infectious) & Infectious == 0,
                               n_max + 10, Removed + Infectious)) %>%
    select(-Susceptible) %>%
    pivot_longer(c(Exposed, Infectious, Removed), names_to="status", values_to="n") %>%
    group_by(status) %>%
    mutate(n = ifelse(cumsum(n) == sum(n), n_max + 10, n))

  line_plot <- ggplot(timeseries) +
    aes(x=time, y=n, fill=status) +
    geom_area() +
    geom_segment(aes(xend=time_max, yend=n), linetype=2, colour="grey", data=segment_df) +
      geom_label(aes(x=x + xnudge, y=n + ynudge, label=status),
                 vjust="middle", hjust="left", size=6, label.padding=unit(0.5, "lines"), colour="white",
                 data=label_df) +
    transition_reveal(time) +
    scale_fill_manual(values=mycols) +
    coord_cartesian(clip="off", expand=FALSE, xlim=c(0, time_max), ylim=c(0, n_max)) +
    xlab("Time") + ylab("Proportion of Individuals") +
    theme_minimal() +
    theme(legend.position="none") +
    theme(plot.margin=margin(5.5, 115, 5.5, 5.5),
          axis.text=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title=element_text(size=16),
          plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                    vjust = 1, margin = margin(b = 5.5)), plot.title.position = "panel",
          plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(b = 5.5)))
  line_gif <- animate(line_plot, end_pause=30, width=600, height=450)

  new_gif <- combine_gifs(net_gif, line_gif)

  image_write(new_gif, filename)
}

EEB_nets <- generate_EEB_networks("../../code/")

g<-as_tbl_graph(EEB_nets$office)
h<-as_tbl_graph(EEB_nets$lab)
full_graph<-graph_join(g, h, by="name") %>% to_undirected()

n_max <- with_graph(full_graph, graph_order()) # note this is different for office net

lab <- make_composite_disease_gif(h, "Just Shared Lab Space", "Lab_Network.gif")
full <- make_composite_disease_gif(full_graph, "Combined Lab and Office", "Full_Network.gif")
# fig <- combine_gifs(full, lab, TRUE)
# image_write(fig, "Full_Figure.gif")
