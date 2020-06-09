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

source("../../simulate_disease_on_network.R")
source("../../generate_EEB_networks.R")

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
  p <- image_append(c(plot1[1], plot2[1]))
  # grow the combined plot frame by frame
  n=ifelse(n1 == n2, n1, n1 * n2)
  n=min(1000, n)  # set max to 1000
  for (i in 2:(n-1)) {
    tmp <- image_append(c(plot1[i], plot2[i]), stack=vertical)
    p <- c(p, tmp)
  }
  return(p)
}

make_row <- function(network, title1, subtitle="", title2) {
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

  net <- ggraph(output %N>% rename(status = !!sym(str_c(time_max))) %>%
                  mutate(           status=factor(status, levels=c("S", "E", "I", "R"),
                                                  labels=c("Susceptible", "Exposed", "Infectious", "Removed"))),
                layout=tmp) +
    geom_edge_link(edge_width=0.66, colour="#635E5B") +
    geom_node_point(aes(colour=I(mycols[status])), size=3.5, show.legend=FALSE) +
    ggtitle(title1, subtitle) +
    theme_graph() +
    theme(plot.margin=margin(0,0,0,0),
          plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                    vjust = 1, margin = margin(b = 0)), plot.title.position = "panel",
          plot.subtitle = element_text(size=rel(1), hjust = 0, vjust = 1, margin = margin(b = 0)))

  timeseries <- tidy_output %>%
    count(status, time) %>%
    complete(status, nesting(time), fill=list(n=0))

  area <- ggplot(timeseries) +
    aes(x=time, y=n, fill=status) +
    geom_area() +
    scale_fill_manual(values=mycols) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125)) +
    coord_cartesian(clip="off", expand=FALSE) +
    xlab("Time") + ylab("Proportion of Individuals") +
    ggtitle(title2) +
    theme_minimal() +
    theme(plot.margin=margin(0,0,0,0),
          axis.text=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          axis.title=element_text(size=16),
          plot.title = element_text(size = rel(1.2), hjust = 0, face="bold", family="Computer Modern",
                                    vjust = 1, margin = margin(b = 0)), plot.title.position = "panel",
          plot.subtitle = element_text(hjust = 0, vjust = 1, margin = margin(b = 0)))

  net + area +
    plot_layout(widths = c(2, 3))
}

EEB_nets <- generate_EEB_networks("../../")

g<-as_tbl_graph(EEB_nets$office)
h<-as_tbl_graph(EEB_nets$lab)
full_graph<-graph_join(g, h, by="name") %>% to_undirected()

make_row(full_graph, title1="A)", title2="B)") / make_row(h, title1="C)", title2="D)") +
  plot_layout(guides="collect") & theme(legend.position='bottom',
                                        legend.text=element_text(size=14),
                                        legend.title=element_blank())
ggsave("static_figure.png", width=12, height=10)
ggsave("static_figure.tiff", width=12, height=10)
