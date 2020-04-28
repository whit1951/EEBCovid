library(magrittr)
library(tidygraph)
library(tidyverse)
library(ggraph)
library(magick)
library(gganimate)
library(igraph)
library(Matrix)
library(tools)

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

run_disease_network_simulation <- function(timesteps, contact_network, model_structure,
                                           beta, sigma=0, gamma=0, xi=0) {
  # parse model structure
  model_states <- str_split(model_structure, "")[[1]]
  if (head(model_states, 1) == tail(model_states, 1)) model_states %<>% head(-1)

  # initialize with one infected node (chosen at random)
  if (is.matrix(contact_network)) contact_network %<>% as_tbl_graph()
  contact_network %<>%
    activate(nodes) %>%
    mutate("0"=sample(c("I", rep("S", nrow(as_tibble(.)) - 1))) %>%
             factor(levels=model_states))

  # progress through the discrete time epidemic simulation
  for (timestep in 1:timesteps) {
    contact_network %<>%
      activate(edges) %>%
      # identify which links have the potential to produce new infectious individuals
      # i.e. which links connect an infectious individual with a susceptible one?
      # note that if useing a directed graph, only use the first clause of this "or"
      mutate(pot_inf_link=or(and(.N()[.E()$from, str_c(timestep - 1)] == "I",
                                   .N()[.E()$to, str_c(timestep - 1)]   == "S"),
                               and(.N()[.E()$from, str_c(timestep - 1)] == "S",
                                   .N()[.E()$to, str_c(timestep - 1)]   == "I")) %>%
               as.integer()) %>%
      activate(nodes) %>%
      # see how many potentially infecting connections each individual has
      mutate(number_infectious_neighbors =
               centrality_degree(weights=pot_inf_link, mode="in"),
             # the probability of a state-change this timestep is determined by
             # the model parameters:
             prob_change_status=case_when(
               !!sym(str_c(timestep - 1)) == "S" ~
                 1 - (1 - beta)^(number_infectious_neighbors), # infection rate
               !!sym(str_c(timestep - 1)) == "E" ~ sigma,      # rate of disease progression
               !!sym(str_c(timestep - 1)) == "I" ~ gamma,      # recovery rate
               !!sym(str_c(timestep - 1)) == "R" ~ xi,         # rate of immunity loss
               # if there is not a rate parameter for this state, then don't change
               TRUE ~ 0),
             # create a new column for this timestep, updating statuses when neccessary
             !!str_c(timestep) :=
               model_states[!!sym(str_c(timestep - 1)) %>%
                              as.integer() %>%
                              add(rbernoulli(n(), prob_change_status)) %>%
                              # in some models, we loop around to the beginning again
                              {ifelse(. > length(model_states), 1, .)}] %>%
               factor(levels=model_states))
  }
  return(contact_network)
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
}

EEBdata<-read.csv("EEB-individ-data-cleaned.csv")
EEBdata$OfficeComb<- with(EEBdata, paste(O.bdg, Office))
EEBdata$LabComb<- with(EEBdata, paste(L.bdg, Lab))

## Office Layer
office_el<-cbind(ID=EEBdata$UNIQUE, Office=EEBdata$OfficeComb)
office_el<-as.data.frame(office_el)
office_el<-distinct(office_el) #remove duplicate entries
office_el<-office_el[!(office_el$Office == " "),] #remove IDs with no office affiliation

#Create bipartite network using sparse matrices
Office_bi <- spMatrix(nrow=length(unique(office_el$ID)),
                      ncol=length(unique(office_el$Office)),
                      i=as.numeric(factor(office_el$ID)),
                      j=as.numeric(factor(office_el$Office)),
                      x=rep(1, length(as.numeric(office_el$ID))) )
row.names(Office_bi) <- levels(factor(office_el$ID))
colnames(Office_bi) <- levels(factor(office_el$Office))

#Create a unipartite network where individuals have edges if they are assigned the same office room
Office_uni <- tcrossprod(Office_bi)

#plot the network
g<-graph_from_adjacency_matrix(Office_uni, mode="undirected", weighted=NULL, diag=FALSE,  add.colnames=NULL, add.rownames=NA)
make_composite_disease_gif(g, "Office Network", "figures/Office_Network.gif")

## Lab Layer
lab_el<-EEBdata %>% select("ID"=UNIQUE, "Lab"=GROUP)
lab_el<-distinct(lab_el) #remove duplicate entries
lab_el<-lab_el[!(lab_el$Lab == " "),] #remove IDs with no office affiliation

#Create bipartite network using sparse matrices
Lab_bi <- spMatrix(nrow=length(unique(lab_el$ID)),
                   ncol=length(unique(lab_el$Lab)),
                   i=as.numeric(factor(lab_el$ID)),
                   j=as.numeric(factor(lab_el$Lab)),
                   x=rep(1, length(as.numeric(lab_el$ID))) )
row.names(Lab_bi) <- levels(factor(lab_el$ID))
colnames(Lab_bi) <- levels(factor(lab_el$Lab))


#Create a unipartite network where individuals have edges if they are assigned the same office room
Lab_uni <- tcrossprod(Lab_bi)

#plot the network
h<-graph_from_adjacency_matrix(Lab_uni, mode="undirected", weighted=NULL, diag=FALSE,  add.colnames=NULL, add.rownames=NA)
make_composite_disease_gif(h, "Lab Network", "figures/Lab_Network.gif")

## Combine Layers
g<-as_tbl_graph(g)
h<-as_tbl_graph(h)
full_graph<-graph_join(g, h, by="name") %>% to_undirected()
make_composite_disease_gif(full_graph, "Combined Office and Lab Network", "figures/Full_Network.gif")
