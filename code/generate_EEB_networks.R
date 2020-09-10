library(tidygraph)
library(igraph)
library(Matrix)

generate_EEB_networks <- function(alternative_path="") {
  EEBdata<-read.csv(paste0(alternative_path, "../data/EEB-individ-data-cleaned.csv"))
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

  ## Lab Layer
  lab_el<-EEBdata %>% select("ID"=UNIQUE, "Lab"=GROUP)
  lab_el<-distinct(lab_el) #remove duplicate entries
  lab_el<-lab_el[!(lab_el$Lab == " "),] #remove IDs with no lab affiliation

  #Create bipartite network using sparse matrices
  Lab_bi <- spMatrix(nrow=length(unique(lab_el$ID)),
                     ncol=length(unique(lab_el$Lab)),
                     i=as.numeric(factor(lab_el$ID)),
                     j=as.numeric(factor(lab_el$Lab)),
                     x=rep(1, length(as.numeric(lab_el$ID))) )
  row.names(Lab_bi) <- levels(factor(lab_el$ID))
  colnames(Lab_bi) <- levels(factor(lab_el$Lab))


  #Create a unipartite network where individuals have edges if they are assigned the same lab room
  Lab_uni <- tcrossprod(Lab_bi)

  #plot the network
  h<-graph_from_adjacency_matrix(Lab_uni, mode="undirected", weighted=NULL, diag=FALSE,  add.colnames=NULL, add.rownames=NA)

  return(list(office=g, lab=h))
}
