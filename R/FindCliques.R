
#' Clique Partition
#'
#' @param clique.set Input set of cliques which may overlap
#' @param randomize Whether we randomize the input list
#'
#' @return A list of non-overlapping cliques
#' @export
#'
#'
CliquePartition <- function(clique.set, randomize = F){
  K <- length(clique.set)
  idx <- 1:K
  if(randomize){
    idx <- sample(idx,K)
  }

  included.cliques <- c()
  included.idx <- c()
  for(k in idx){
    clq.tmp <- clique.set[[k]]
    if(!(any(clq.tmp %in% included.cliques))){
      included.cliques <- c(included.cliques,clq.tmp)
      included.idx <- c(included.idx, k)
    }
  }
  reduced.cliques <- list()
  k.sub <- 1
  for(k in included.idx){
    reduced.cliques[[k.sub]] <- clique.set[[k]]
    k.sub <- k.sub + 1
  }

  return(reduced.cliques)
}







#' Search For Cliques
#'
#' @param G Adjacency Matrix
#' @param min_clique_size Minimum clique size
#' @param res Tuning parameter for clustering
#' @param verbose Whether to include additional messages
#'
#' @return A list of non-overlapping cliques
#' @export
#'
#'
ClusterCliqueSearch <- function(G, min_clique_size = 8, res = 1, verbose = F){
  res.tmp = res
  n = nrow(G)
  clique.set <- list()
  iter = 1

  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  comm <- igraph::cluster_leading_eigen(g)
  table(comm$membership)
  K <- max(comm$membership)
  clique.set <- list()
  for(k in seq(K)){
    #print(paste("Block:", k,"/",K))
    idx.sub <- which(as.numeric(comm$membership) == k )
    G.sub <- G[idx.sub,idx.sub]
    if(length(idx.sub) > 1){
      n.sub <- nrow(G.sub)
      g.sub <- igraph::graph_from_adjacency_matrix(G.sub, mode= "undirected")
      clique.set.sub <- igraph::maximal.cliques(g.sub, min = min_clique_size)
      if(length(clique.set.sub) > 0){
        clique.set.sub <- CliquePartition(clique.set.sub)
        cliques.sub.re.idx <- list()
        R <- length(clique.set.sub)
        for(r in seq(R)){
          cliques.sub.re.idx[[r]] <- idx.sub[clique.set.sub[[r]]]
        }
        clique.set <- append(clique.set,cliques.sub.re.idx)
      }

    }
  }
  if(verbose){
    print(paste("Number of Cliques of size,",min_clique_size,":", length(clique.set)))
  }
  return(clique.set)
}


#' Exhaustive Search For Cliques
#'
#' @param G Adjacency Matrix
#' @param min_clique_size Minimum clique size
#' @param verbose Whether to include additional messages
#'
#' @return A list of non-overlapping cliques
#' @export
#'
#'
CliqueSearch <- function(G, min_clique_size = 8, verbose = F){
  g <- igraph::graph_from_adjacency_matrix(G, mode= "undirected")
  clique.set = igraph::max_cliques(g,min = min_clique_size)
  clique.set = CliquePartition(clique.set)
  if(verbose){
    print(paste("Number of Cliques of size,",min_clique_size,":", length(clique.set)))
  }
  return(clique.set)
}




