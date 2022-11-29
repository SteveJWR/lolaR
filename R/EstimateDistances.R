

#' Estimate Random Effects By Degree Ratios
#'
#' @param G Full Data Adjacency Matrix
#' @param cliques List of indices corresponding to cliques
#'
#' @return Estimated Set of Random Effects
#' @export
#'
#'
estimateRandeff <- function(G,cliques){
  K <- length(cliques)
  nu.set <- list()
  for(k in 1:K){
    p.hats <- Matrix::colMeans(G[,as.vector(cliques[[k]])])
    nu.hats <- log(p.hats) - log(max(p.hats))
    nu.set[[k]] <- nu.hats
  }
  return(nu.set)
}



#' Floyd-Warshall Algorithm
#'
#' @param D0.init Initialization for Floyd-Warshall Algorithm
#'
#' @return Distance Matrix which satisfied triangle inequality
#' @export
#'
# TODO: Add @examples
FWA <- function(D0.init){
  # first pass D0 may not be a metric
  min.d0 <- min(D0.init)
  D0 = D0.init - min.d0
  diag(D0) = 0
  K <- nrow(D0)
  for(k in seq(K)){
    for(j in seq(K)){
      for(i in seq(j)){
        if(D0[i,j] > D0[i,k] + D0[k,j]){
          D0[i,j] = D0[i,k] + D0[k,j]
          D0[j,i] = D0[i,j]
        }
      }
    }
  }
  max.d = max(D0[!is.infinite(D0)])
  # trimming unconnected cliques
  D0[D0 > max.d] = max.d
  return(D0)
}


#' Count of Cliques
#'
#' @param cliques A list of clique indices
#'
#' @return List of indices corresponding to cliques and corresponding vector of labels
#' @export
#'
# TODO: Add @examples
CliquesCount <- function(cliques){
  clique.vec <- c()
  cliques.idx <- c()
  K = length(cliques)
  for(k in seq(K)){
    clique.vec <- c(clique.vec, rep(k, length(cliques[[k]])))
    cliques.idx <- c(cliques.idx, as.numeric(cliques[[k]]))
  }

  return(list("idx" = cliques.idx, "labels" = clique.vec))
}


#' Good Initialization For Distance Matrix
#'
#' @param G Full Data Adjacency Matrix
#' @param cliques List of indices corresponding to cliques
#' @param fixed.effect.vec Vector of corresponding estimated random effects
#' @param global.randeff Whether to approximate and group the random effects into a single one
#' @param add.offset Add a minimum distance for better optimization
#'
#' @return Initialization of Estimated Distance Matrix which obeys triangle inequalities
#' @export
#'
# TODO: Add @examples
InitD0 <- function(G, cliques, fixed.effect.vec, global.randeff = T, add.offset = 0.001){

  cliques.subset <- CliquesCount(cliques)

  G.subset = G[cliques.subset$idx, cliques.subset$idx]
  if(global.randeff){
    K = length(cliques)
    D0 = matrix(0, K,K)
    glb.rndf <- log(mean(exp(fixed.effect.vec)))
    for(k1 in seq(K)){
      for(k2 in seq(K)){
        idx1 = which(cliques.subset$labels == k1)
        idx2 = which(cliques.subset$labels == k2)
        p.xy = Matrix::mean(G.subset[idx1,idx2])
        D0[k1,k2] = -log(p.xy) + 2*glb.rndf
      }
    }
  } else {
    K = length(cliques)
    D0 = matrix(0, K,K)
    for(k1 in seq(K)){
      for(k2 in seq(K)){
        idx1 = which(cliques.subset$labels == k1)
        idx2 = which(cliques.subset$labels == k2)
        p.xy = Matrix::mean(G.subset[idx1,idx2])
        D0[k1,k2] = -log(p.xy) + log(mean(exp(fixed.effect.vec[idx1]))) + log(mean(exp(fixed.effect.vec[idx2])))
      }
    }
  }

  D0[D0 < 0 ] = 0
  D0 = FWA(D0)

  D0 = D0 + add.offset
  diag(D0) = 0

  return(D0)
}





#' Value of Clique Subset Likelihood
#' @export
likelihood_value <- function(G, cliques, fixed.effect.vec, D0){

  K = length(cliques)

  cliques.subset = CliquesCount(cliques)
  cliques.subset.labels = cliques.subset$labels
  cliques.subset.idx = cliques.subset$idx
  G.subset <- G[cliques.subset.idx, cliques.subset.idx]

  nu.big <- outer(fixed.effect.vec,fixed.effect.vec, "+")
  obj.val <- 0

  K.large <- dim(G.subset)[1]
  D.big <- matrix(0,K.large,K.large)

  #full.idx <- seq(length(cliques))
  for(k1 in seq(K)){
    for(k2 in seq(K)){
      idx1 = which(cliques.subset.labels == k1)
      idx2 = which(cliques.subset.labels == k2)
      D.big[idx1,idx2] = D0[k1,k2]
      # nu.block <- nu.big[idx1,idx2]
      # A.block <- G[idx1,idx2]
    }
  }


  lik.mat = G.subset*(nu.big - D.big) + (1 - G.subset)*log(1 - exp(nu.big - D.big))
  lik.mat[is.infinite(lik.mat)] = 0 # not counting self-clique connections
  lik.mat[is.na(lik.mat)] = 0 # not counting self-clique connections
  obj.val <- sum(lik.mat, na.rm = T) # removing the erroneous terms
  return(obj.val)
}



#' Estimate Distance Matrix From Cliques
#'
#' @param G Full Data Adjacency Matrix
#' @param cliques List of indices corresponding to cliques
#' @param fixed.effect.vec List of estimated fixed Effects
#' @param D0 Initial Point To Take Taylor Series of Objective
#' @param thresh Default Threshold for Increasing Likelihood
#' @param max.iter Maximum number of iterations of the second order approximation
#' @param solver Which solver should be passed to CVX
#' @param verbose Include additional messages
#'
#' @return Estimated Distance Matrix
#' @export
#'
EstimateD <- function(G, cliques, D0,
                      thresh = 10**(-6), max.iter = 50, solver = "MOSEK",
                      verbose = F, rand.eff.0 = F){

  # numerical smoothing for some gradient terms.
  eps = 10**(-9) # precision for terms

  if(rand.eff.0) {
    fixed.effect.vec <- rep(0, length(unlist(cliques)))
  } else {
    fixed.effect.vec <- unlist(estimateRandeff(G, cliques))
  }


  if(solver %in% c("MOSEK","GUROBI")){
    if(! solver %in% CVXR::installed_solvers()){
      cvx_solver = "OSQP"
    } else {
      cvx_solver = solver
    }
  } else {
    cvx_solver = solver
  }

  K = length(cliques)
  if(missing(D0)){
    D0 <- InitD0(G, cliques, fixed.effect.vec) # provide an initial value
  }


  D <- CVXR::Variable(K,K, name = "Distance Matrix")

  #D.big <- Variable(n.subg,n.subg, name = "Distance Matrix")
  nu.big <- outer(fixed.effect.vec,fixed.effect.vec, "+")
  #d.vec <- Variable(K^2, name = "Distance Vector")

  d.vec = CVXR::vec(D) #define a vectorization of the distance matrix
  constraints <- list(
    CVXR::diag(D) == 0,
    D == t(D),
    D >= 0
  )

  # the triangle inequalities

  index.block <- expand.grid(seq(K),seq(K),seq(K))
  ib.1 <- index.block[,1] < index.block[,2]
  ib.2 <- index.block[,2] < index.block[,3]
  index.block <- index.block[ib.1 & ib.2, ]


  ref.block <- index.block
  ref.block[,1] <- (index.block[,1] - 1)*K + index.block[,2]
  ref.block[,2] <- (index.block[,2] - 1)*K + index.block[,3]
  ref.block[,3] <- (index.block[,1] - 1)*K + index.block[,3]

  E.block.1 <- matrix(0,nrow(index.block),K^2)
  E.block.1 <- as(E.block.1, "sparseMatrix")
  n.ref <- nrow(ref.block)

  E.block.1 <- Matrix::sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
                                    j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
                                    x = c(rep(1,n.ref),rep(1,n.ref),rep(-1,n.ref)),
                                    dims = c(n.ref,K^2))
  E.block.2 <- Matrix::sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
                                    j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
                                    x = c(rep(1,n.ref),rep(-1,n.ref),rep(1,n.ref)),
                                    dims = c(n.ref,K^2))
  E.block.3 <- Matrix::sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
                                    j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
                                    x = c(rep(-1,n.ref),rep(1,n.ref),rep(1,n.ref)),
                                    dims = c(n.ref,K^2))

  E.mat <- rbind(E.block.1,
                 E.block.2,
                 E.block.3)
  # These are the biggest memory objects, so we don't want to duplicate them
  rm(E.block.1)
  rm(E.block.2)
  rm(E.block.3)

  # turn the many restrictions into a single one.
  constraints <- append(constraints, E.mat %*% d.vec >= 0)


  cliques.subset = CliquesCount(cliques)
  cliques.subset.idx = cliques.subset$idx

  G.subset <- G[cliques.subset.idx, cliques.subset.idx]
  ######
  # for fast computation of M and B matrices
  clique.sizes = table(cliques.subset$labels)
  l.max = max(clique.sizes)
  Id.tens = array(0,c(K,K,l.max**2))
  G.tens = array(0,c(K,K,l.max**2)) #tensor version of the clique adjacency
  nu.tens = array(0,c(K,K,l.max**2))

  # each pair has a length
  for(k1 in seq(K)){
    for(k2 in seq(K)){
      if(k1 != k2){
        n.potential.connections = clique.sizes[k1]*clique.sizes[k2]
        idx1 = which(cliques.subset$labels == k1)
        idx2 = which(cliques.subset$labels == k2)
        nu.block <- nu.big[idx1,idx2]

        Id.tens[k1,k2,1:n.potential.connections] = 1
        nu.tens[k1,k2,1:n.potential.connections] = as.vector(nu.block)
        G.tens[k1,k2,1:n.potential.connections] = as.vector(G.subset[idx1,idx2])
        #D.prev.tens[k1,k2,] = D.prev[k1,k2]
      }
    }
  }

  ######
  # successive second order approximation
  # allows for a fast implementation of QP solutions
  D.prev <- D0
  D.prev.prev <- D0
  diff = Inf
  lik.rat = 1
  iter = 1
  lik.prev = -Inf
  while(lik.rat > thresh & iter < max.iter){
    time.1 <- Sys.time()
    D.prev.tens = array(0,c(K,K,l.max**2))
    for(k1 in seq(K)){
      for(k2 in seq(K)){
        if(k1 != k2){
          D.prev.tens[k1,k2,] = D.prev[k1,k2]
        }
      }
    }

    M.tens <- ((1 - G.tens)*(-1/2)*(exp(nu.tens - D.prev.tens - eps))/((1 - exp(nu.tens - D.prev.tens - eps))**2))*Id.tens

    B.tens <- (-G.tens + (1 - G.tens)*(exp(nu.tens - D.prev.tens - eps))/(1 - exp(nu.tens - D.prev.tens - eps)) - 2*M.tens*D.prev.tens)*Id.tens


    M = rowSums(M.tens, dims = 2)
    B = rowSums(B.tens, dims = 2)

    diag(M) = 0
    diag(B) = 0

    time.2 <- Sys.time()
    # print(time.2 - time.1)


    b = as.vector(B)
    w = sqrt(as.vector(-M))
    obj.arg = -CVXR::sum_squares(w*d.vec) + sum(b * d.vec)

    obj <- CVXR::Maximize(obj.arg)
    prob <- CVXR::Problem(obj, constraints)

    time.1 <- Sys.time()
    # can get weirdly stuck in numerical error

    result <- tryCatch({
      CVXR::solve(prob, solver = cvx_solver)
    }, error = function(e) {

      return(NULL)})

    time.2 <- Sys.time()
    #print(time.2 - time.1)
    # unstuck the problem sometimes

    if(!is.null(result)){
      if(result$status == "solver_error"){
        D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
        lik.prev = -Inf
      } else {
        D.next <- result$getValue(D)
      }

    } else {
      D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
      lik.prev <- -Inf
    }

    D.next[D.next < eps] = 0
    # large fraction at zero, we should reset
    if(mean(D.next == 0) > 0.5){
      D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
      lik.prev <- -Inf
    }

    D.next[D.next < eps] = 0
    # TODO: Update the likelihood_value function
    lik.next = likelihood_value(G, cliques, fixed.effect.vec, D.next)
    diff = lik.next - lik.prev
    lik.rat = abs(1 - lik.next/lik.prev)
    if(is.nan(diff) | lik.prev == -Inf | is.nan(lik.rat)){
      diff = Inf
      lik.rat = 1
    }

    lik.prev = lik.next
    D.prev.prev = D.prev
    D.prev = D.next

    if(verbose){
      cat(paste("Num Steps:", iter, "Likelihood Stopping Criteria:", round(lik.rat,7)), end = "\r")
    }
    iter = iter + 1
  }
  D.hat = D.next
  return(D.hat)
}
# EstimateD <- function(G, cliques, D0,
#                       thresh = 10**(-3), max.iter = 50, solver = "MOSEK",
#                       verbose = F, rand.eff.0 = F){
#
#   if(rand.eff.0) {
#     fixed.effect.vec <- rep(0, length(unlist(cliques)))
#   } else {
#     fixed.effect.vec <- unlist(estimateRandeff(G, cliques))
#   }
#
#
#
#
#   # numerical smoothing for some gradient terms.
#   eps = 10**(-8) # precision for terms
#
#   if(solver %in% c("MOSEK","GUROBI")){
#     if(! solver %in% installed_solvers()){
#       cvx_solver = "OSQP"
#     } else {
#       cvx_solver = solver
#     }
#   } else {
#     cvx_solver = solver
#   }
#
#   K = max(cliques)
#   if(missing(D0)){
#     D0 <- InitD0(G, cliques, fixed.effect.vec) # provide an initial value
#   }
#
#
#   D <- Variable(K,K, name = "Distance Matrix")
#
#   #D.big <- Variable(n.subg,n.subg, name = "Distance Matrix")
#   nu.big <- outer(fixed.effect.vec,fixed.effect.vec, "+")
#   #d.vec <- Variable(K^2, name = "Distance Vector")
#
#   d.vec = vec(D) #define a vectorization of the distance matrix
#   constraints <- list(
#     diag(D) == 0,
#     D == t(D),
#     D >= 0
#   )
#
#   # the triangle inequalities
#
#   index.block <- expand.grid(seq(K),seq(K),seq(K))
#   ib.1 <- index.block[,1] < index.block[,2]
#   ib.2 <- index.block[,2] < index.block[,3]
#   index.block <- index.block[ib.1 & ib.2, ]
#
#
#   ref.block <- index.block
#   ref.block[,1] <- (index.block[,1] - 1)*K + index.block[,2]
#   ref.block[,2] <- (index.block[,2] - 1)*K + index.block[,3]
#   ref.block[,3] <- (index.block[,1] - 1)*K + index.block[,3]
#
#   E.block.1 <- matrix(0,nrow(index.block),K^2)
#   E.block.1 <- as(E.block.1, "sparseMatrix")
#   n.ref <- nrow(ref.block)
#   E.block.1 <- sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
#                             j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
#                             x = c(rep(1,n.ref),rep(1,n.ref),rep(-1,n.ref)),
#                             dims = c(n.ref,K^2))
#   E.block.2 <- sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
#                             j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
#                             x = c(rep(1,n.ref),rep(-1,n.ref),rep(1,n.ref)),
#                             dims = c(n.ref,K^2))
#   E.block.3 <- sparseMatrix(i = c(seq(n.ref),seq(n.ref),seq(n.ref)),
#                             j = c(ref.block[,1],ref.block[,2],ref.block[,3]),
#                             x = c(rep(-1,n.ref),rep(1,n.ref),rep(1,n.ref)),
#                             dims = c(n.ref,K^2))
#
#   E.mat <- rbind(E.block.1,
#                  E.block.2,
#                  E.block.3)
#
#   rm(E.block.1)
#   rm(E.block.2)
#   rm(E.block.3)
#
#   # turn the many restrictions into a single one.
#   constraints <- append(constraints, E.mat %*% d.vec >= 0)
#
#
#
#   ######
#   # for fast computation of M and B matrices
#   clique.sizes = table(cliques)
#   l.max = max(clique.sizes)
#   Id.tens = array(0,c(K,K,l.max**2))
#   G.tens = array(0,c(K,K,l.max**2)) #tensor version of the clique adjacency
#   nu.tens = array(0,c(K,K,l.max**2))
#
#   # each pair has a length
#   for(k1 in seq(K)){
#     for(k2 in seq(K)){
#       if(k1 != k2){
#         n.potential.connections = clique.sizes[k1]*clique.sizes[k2]
#         idx1 = which(cliques == k1)
#         idx2 = which(cliques == k2)
#         nu.block <- nu.big[idx1,idx2]
#
#         Id.tens[k1,k2,1:n.potential.connections] = 1
#         nu.tens[k1,k2,1:n.potential.connections] = as.vector(nu.block)
#         G.tens[k1,k2,1:n.potential.connections] = as.vector(G[idx1,idx2])
#         #D.prev.tens[k1,k2,] = D.prev[k1,k2]
#       }
#     }
#   }
#
#   ######
#   # successive second order approximation
#   # allows for a fast implementation of QP solutions
#   D.prev <- D0
#   D.prev.prev <- D0
#   diff = Inf
#   iter = 1
#   lik.prev = -Inf
#   while(diff > thresh & iter < max.iter){
#     time.1 <- Sys.time()
#     D.prev.tens = array(0,c(K,K,l.max**2))
#     for(k1 in seq(K)){
#       for(k2 in seq(K)){
#         if(k1 != k2){
#           D.prev.tens[k1,k2,] = D.prev[k1,k2]
#         }
#       }
#     }
#
#     M.tens <- ((1 - G.tens)*(-1/2)*(exp(nu.tens - D.prev.tens - eps))/((1 - exp(nu.tens - D.prev.tens - eps))**2))*Id.tens
#
#     B.tens <- (-G.tens + (1 - G.tens)*(exp(nu.tens - D.prev.tens - eps))/(1 - exp(nu.tens - D.prev.tens - eps)) - 2*M.tens*D.prev.tens)*Id.tens
#
#
#     M = rowSums(M.tens, dims = 2)
#     B = rowSums(B.tens, dims = 2)
#
#     diag(M) = 0
#     diag(B) = 0
#
#     time.2 <- Sys.time()
#     # print(time.2 - time.1)
#
#
#     b = as.vector(B)
#     w = sqrt(as.vector(-M))
#     obj.arg = -sum_squares(w*d.vec) + sum(b * d.vec)
#
#     obj <- Maximize(obj.arg)
#     prob <- Problem(obj, constraints)
#
#     time.1 <- Sys.time()
#     # can get weirdly stuck in numerical error
#
#     result <- tryCatch({
#       solve(prob, solver = cvx_solver)
#     }, error = function(e) {
#
#       return(NULL)})
#
#     time.2 <- Sys.time()
#     #print(time.2 - time.1)
#     # unstuck the problem sometimes
#
#     if(!is.null(result)){
#       if(result$status == "solver_error"){
#         D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
#         lik.prev = -Inf
#       } else {
#         D.next <- result$getValue(D)
#       }
#
#     } else {
#       D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
#       lik.prev <- -Inf
#     }
#
#     D.next[D.next < eps] = 0
#     # large fraction at zero, we should reset
#     if(mean(D.next == 0) > 0.5){
#       D.next <- D.prev + (0.2)*(D.prev - D.prev.prev)
#       lik.prev <- -Inf
#     }
#
#     D.next[D.next < eps] = 0
#     lik.next = likelihood_value(G, cliques, fixed.effect.vec, D.next)
#     diff = lik.next - lik.prev
#     if(is.nan(diff) | lik.prev == -Inf){
#       diff = Inf
#     }
#     lik.prev = lik.next
#     D.prev.prev = D.prev
#     D.prev = D.next
#
#     if(verbose){
#       cat(paste("Num Steps:", iter, "Diff Likelihood:", round(diff,4)), end = "\r")
#     }
#     iter = iter + 1
#   }
#
#   D.hat = D.next
#
#   return(D.hat)
# }


