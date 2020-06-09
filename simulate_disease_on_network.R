library(magrittr)
library(tidygraph)
library(tidyverse)

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
      mutate(pot_inf_link=edge_is_between(.N()[str_c(timestep - 1)] == "S",
                                          .N()[str_c(timestep - 1)] == "I")) %>%
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
