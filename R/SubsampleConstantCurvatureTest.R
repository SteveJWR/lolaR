
#Upper bound on curvature Estimate
# TODO: hide these functions from user interface
#' Estimating function for the upper bound on curvature
#'
#'
#' @param kappa curvature
#' @param dxy distance between points x and y
#' @param dxz distance between points x and z
#' @param dyz distance between points y and z
#' @param dxm distance between points x and m
#' @param dym distance between points y and m
#' @param dzm distance between points z and m
#'
#' @return curvature estimating function bound
#' @export
#'
# @examples
g_u <- function(kappa, dxy, dxz, dyz, dxm, dym, dzm){
  out <- gEF(kappa,dxy,dxz,dyz,dxm + MidDist(kappa,dym,dzm,dyz))
}

# TODO: hide these functions from user interface
#' Estimating function for the lower bound on curvature
#'
#' @param kappa curvature
#' @param dxy distance between points x and y
#' @param dxz distance between points x and z
#' @param dyz distance between points y and z
#' @param dxm distance between points x and m
#' @param dym distance between points y and m
#' @param dzm distance between points z and m
#'
#' @return  curvature estimating function bound
#' @export
#'
# @examples
g_l <- function(kappa, dxy, dxz, dyz, dxm, dym, dzm){
  out <- gEF(kappa,dxy,dxz,dyz,dxm - MidDist(kappa,dym,dzm,dyz))
}



#' Estimator of upper bound on the curvature
#'
#' @param dxy distance between points x and y
#' @param dxz distance between points x and z
#' @param dyz distance between points y and z
#' @param dxm distance between points x and m
#' @param dym distance between points y and m
#' @param dzm distance between points z and m
#' @param kappa.prec granularity of kappa grid
#' @param min.curvature minimum curvature value searched
#'
#' @return  Upper Bound Estimate
#' @export
#'
# @examples
kappa_u <- function(dxy, dxz, dyz, dxm, dym, dzm,
                    kappa.prec = 10**(-5),
                    min.curvature = -1000){

  # first check for whether the midpoint estimate is already too far for any curvature:

  if(dxm < min(dxy,dxz) - (1/2)*dyz){
    warning("Triangle Inequality is not satisfied")
    return(-Inf)
  }


  # Picking good initialization for the grid search.
  max.curvature = (pi/max(c(dxy,dxz,dyz,dxm,dym,dzm)))**2

  kappa.upper <- max.curvature
  kappa.lower <- min.curvature

  g.u.upper <- g_u(kappa.upper, dxy, dxz, dyz, dxm, dym, dzm)
  if(is.nan(g.u.upper)){
    g.u.upper = 1
  }
  g.u.lower <- g_u(kappa.lower, dxy, dxz, dyz, dxm, dym, dzm)


  if(g.u.upper < 0 ) {
    return(max.curvature)
  } else if(g.u.lower > 0) {
    return(min.curvature)
  } else {
    kappa.gap <- kappa.upper - kappa.lower
    while(kappa.gap > kappa.prec){
      kappa.mid <-  mean(c(kappa.upper, kappa.lower))
      g.u.mid <- g_u(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
      #print(kappa.mid)
      # Handles a few numerical errors
      tries <- 0
      while(is.nan(g.u.mid) & tries < 10){
        #print(tries)

        kappa.mid <- runif(1,kappa.lower, kappa.upper)
        #print(kappa.mid)
        g.u.mid <- g_u(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
        tries <- tries + 1
      }
      if(tries >= 10){
        return(max.curvature)
      }

      if(g.u.mid > 0){
        kappa.upper <- kappa.mid
      } else if(g.u.mid == 0) {
        break
      } else {
        kappa.lower <- kappa.mid
      }
      kappa.gap <- kappa.upper - kappa.lower
    }
    return(kappa.mid)
  }
}

#' Estimator of lower bound on the curvature
#'
#' @param dxy distance between points x and y
#' @param dxz distance between points x and z
#' @param dyz distance between points y and z
#' @param dxm distance between points x and m
#' @param dym distance between points y and m
#' @param dzm distance between points z and m
#' @param kappa.prec granularity of kappa grid
#' @param min.curvature minimum curvature value searched
#'
#' @return Lower bound estimate
#' @export
#'
# @examples
kappa_l <- function(dxy, dxz, dyz, dxm, dym, dzm,
                    kappa.prec = 10**(-5),
                    min.curvature = -1000){

  # first check for whether the midpoint estimate is already too far for any curvature:

  if(dxm < min(dxy,dxz) - (1/2)*dyz){
    warning("Triangle Inequality is not satisfied")
    return(-Inf)
  }


  # Picking good initialization for the grid search.
  max.curvature = (pi/max(c(dxy,dxz,dyz,dxm,dym,dzm)))**2

  kappa.upper <- max.curvature
  kappa.lower <- min.curvature

  g.l.upper <- g_l(kappa.upper, dxy, dxz, dyz, dxm, dym, dzm)
  if(is.nan(g.l.upper)) {
    g.l.upper = 1
  }
  g.l.lower <- g_l(kappa.lower, dxy, dxz, dyz, dxm, dym, dzm)


  if(g.l.upper < 0 ) {
    return(max.curvature)
  } else if(g.l.lower > 0) {
    return(min.curvature)
  } else {
    kappa.gap <- kappa.upper - kappa.lower
    while(kappa.gap > kappa.prec){
      kappa.mid <- mean(c(kappa.upper, kappa.lower))

      g.l.mid <- g_l(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)

      # Handles a few numerical errors
      tries <- 0
      while(is.nan(g.l.mid) & tries < 10){
        kappa.mid <- runif(1,kappa.lower, kappa.upper)
        g.l.mid <- g_l(kappa.mid, dxy, dxz, dyz, dxm, dym, dzm)
        tries <- tries + 1
      }
      if(tries >= 10){
        return(min.curvature)
      }

      if(g.l.mid > 0){
        kappa.upper <- kappa.mid
      } else if(g.l.mid == 0) {
        break
      } else {
        kappa.lower <- kappa.mid
      }
      kappa.gap <- kappa.upper - kappa.lower
    }
    return(kappa.mid)
  }
}


#' estimate curvature bounds using a reference set x
#'
#' @param D Distance Matrix
#' @param y reference point y
#' @param z reference point z
#' @param m surrogate midpoint of y and z
#' @param x.set set of reference x coordinates
#'
#' @return Upper And Lower Bounds for estimate
#' @export
#'
# @examples
estimateBounds <- function(D,y,z,m,x.set){
  kappa.us <- c()
  kappa.ls <- c()

  for(x in x.set){

    dxy <- D[x,y]
    dxz <- D[x,z]
    dyz <- D[y,z]
    dxm <- D[x,m]
    dym <- D[y,m]
    dzm <- D[z,m]

    kap.u <- kappa_u(dxy, dxz, dyz, dxm, dym, dzm,
                     kappa.prec = 10**(-5),
                     min.curvature = -1000)

    kap.l <- kappa_l(dxy, dxz, dyz, dxm, dym, dzm,
                     kappa.prec = 10**(-5),
                     min.curvature = -1000)

    kappa.us <- c(kappa.us,kap.u)
    kappa.ls <- c(kappa.ls,kap.l)
  }
  return(list("upper.bounds" = kappa.us,
              "lower.bounds" = kappa.ls))
}



# Can Change the default settings for the midpoint search.
#' Selection of good midpoints and reference points. Simple default values are included for ease of use
#'
#' @param D Distance Matrix
#' @param J Number of midpoints
#' @param tri.const filter triangle constant
#' @param d.yz.min minimum distance for yz
#' @param d.yz.max maximum distance for yz
#' @param sub.idx subset of indices searched
#'
#' @return Reference points
#' @export
#'
# @examples
selectReference <- function(D,J,tri.const = 1.4, d.yz.min = 1, d.yz.max = 5, sub.idx){
  K = nrow(D)
  if(missing(sub.idx)){
    sub.idx = 1:K
    midpoints <- optimal_midpoint_search(D, top.k = J, d.yz.min = d.yz.min,
                                         d.yz.max = d.yz.max, subset = sub.idx)
    out.list <- list()
    for(j in seq(J)){
      y <- midpoints[j,1]
      z <- midpoints[j,2]
      m <- midpoints[j,3]
      x.set <- filter_indices(D,y,z,m, tri.const = tri.const)
      out.list[[j]] <- list("y" = y, "z" = z, "m" = m, "xset" = x.set)
    }
  } else {
    midpoints <- optimal_midpoint_search(D, top.k = J, d.yz.min = d.yz.min,
                                         d.yz.max = d.yz.max, subset = sub.idx)
    out.list <- list()
    for(j in seq(J)){
      y <- midpoints[j,1]
      z <- midpoints[j,2]
      m <- midpoints[j,3]
      x.set <- filter_indices(D,y,z,m, tri.const = tri.const)
      out.list[[j]] <- list("y" = y, "z" = z, "m" = m, "xset" = x.set)
    }
  }

  return(out.list)
}


#' Constant Curvature Test Based on Subsampling approaches.
#'
#' @param A Adjacency matrix
#' @param clique.set indices of where the cliques are located in A
#' @param reference.set reference set for where curvature is measured
#' @param subsample.rate drop out number for the subsampling distributions
#' @param B Number of subsample iterations
#' @param store.D.samples stores subsamples of the matrix D for later use
#' @param max.iter number of iterations performed for each maximum likelihood step.  Setting this to 1 uses the one-step estimator
#'
#' @return Test of constant curvature based on subsampling
#' @export
#'
# @examples
SubSampleConstantCurvatureTest <- function(A,clique.set,reference.set,
                                           subsample.rate = 1,B = 100, store.D.samples = T,
                                           max.iter = 1){
  # Max-iter set to default at 1 in order to use the one-step estimator
  K = length(clique.set)
  J = length(reference.set)
  upper.bound.sub <- rep(NA, B)
  lower.bound.sub <- rep(NA, B)
  clique.subsample <- list()
  D.hat <- EstimateD(A,clique.set,max.iter = max.iter)
  D.subsample.list <- list()
  for(b in seq(B)){
    for(k in seq(K)){
      ell = length(clique.set[[k]])
      clique.subsample[[k]] <- sample(clique.set[[k]],
                                      size = ell - subsample.rate,
                                      replace = F)
    }
    sub.idx <- CliquesCount(clique.subsample)

    D.hat.sub <- EstimateD(A,clique.subsample, D0 = D.hat,
                           max.iter = max.iter, verbose = F)
    upper.min <- Inf
    lower.max <- -Inf
    if(store.D.samples){
      D.subsample.list[[b]] <- D.hat.sub
    }
    for(j in seq(J)){
      y <- reference.set[[j]][["y"]]
      z <- reference.set[[j]][["z"]]
      m <- reference.set[[j]][["m"]]
      x.set <- reference.set[[j]][["xset"]]
      bounds <- estimateBounds(D.hat.sub, y,z,m,x.set)
      upper.min <- min(c(upper.min, median(bounds$upper.bounds)))
      lower.max <- max(c(lower.max, median(bounds$lower.bounds)))
      upper.bound.sub[b] <- upper.min
      lower.bound.sub[b] <- lower.max
    }
    cat("Subsample", b, "/", B, end = "\r")
  }
  p.value <- ComputePvalue(upper.bound.sub, lower.bound.sub)
  return(list("p.value" = p.value, "upper.bounds" = upper.bound.sub,
              "lower.bounds" = lower.bound.sub, "D.subs" = D.subsample.list))
}

# Recompute the test values using cached distance matrix set.
# This is a function purely for the simulations in the paper.
# This can be only used when A and clique.set were the same that generated the constant curvature test.
# TODO: We probably don't need this in the simulation in the end.
#
#' Compute the Constant Curvature Test using Multiple Thresholds. Used in replicating the simulations.
#'
#' @param D.hat initial estimated distance matrix
#' @param D.subsample list of subsampled distance matrices
#' @param reference.set reference set of selected midpoints and filtered reference points
#' @param tri.const.seq sequence of filtered x.values
#' @param verbose print details
#'
#' @return Test of constant curvature using multiple thresholds
#' @export
#'
# @examples
SubSampleConstantCurvatureTestMultipleThresholds <- function(D.hat,D.subsample, reference.set,
                                                             tri.const.seq, verbose = F){
  output.data <- matrix(rep(NA,2*length(tri.const.seq)), ncol = 2)
  output.data <- as.data.frame(output.data)
  colnames(output.data) <- c("threshold.constant", "p.value")
  B <- length(D.subsample)
  J <- length(reference.set)

  for(k in seq(length(tri.const.seq))){
    if(verbose){
      cat(paste0("Constant: ", k, "/", length(tri.const.seq), end = "\r"))
    }
    tri.const <- tri.const.seq[k]
    reference.set.tmp = reference.set
    for(j in seq(J)){
      y <- reference.set.tmp[[j]]$y
      z <- reference.set.tmp[[j]]$z
      m <- reference.set.tmp[[j]]$m
      x.set.tmp <- filter_indices(D.hat,y,z,m, tri.const = tri.const)
      reference.set.tmp[[j]]$xset = x.set.tmp
    }
    upper.bound.sub <- rep(NA, B)
    lower.bound.sub <- rep(NA, B)
    for(b in seq(B)){
      D.hat.sub <- D.subsample[[b]]
      upper.min <- Inf
      lower.max <- -Inf
      for(j in seq(J)){
        y <- reference.set.tmp[[j]][["y"]]
        z <- reference.set.tmp[[j]][["z"]]
        m <- reference.set.tmp[[j]][["m"]]
        x.set <- reference.set.tmp[[j]][["xset"]]
        bounds <- estimateBounds(D.hat.sub, y,z,m,x.set)

        upper.min <- min(c(upper.min, median(bounds$upper.bounds)))
        lower.max <- max(c(lower.max, median(bounds$lower.bounds)))
        upper.bound.sub[b] <- upper.min
        lower.bound.sub[b] <- lower.max
      }
    }

    p.value <- ComputePvalue(upper.bound.sub, lower.bound.sub)
    output.data[k,] <- c(tri.const,p.value)
  }
  return(output.data)
}



# Overlapping Bounds p-value
#' Compute p-value from overlap of the subsampled upper and lower bounds.
#'
#' @param upper.bounds upper bound subsample set
#' @param lower.bounds lower bound subsample set
#'
#' @return p-value from subsampled upper and lower bounds.
#' @export
#'
# @examples
ComputePvalue <- function(upper.bounds, lower.bounds){
  B = length(upper.bounds)
  sort.u <- sort(upper.bounds, decreasing = F)
  sort.l <- sort(lower.bounds, decreasing = T)

  for(b in seq(B)){
    u = sort.u[b]
    l = sort.l[b]
    if(u > l){
      break
    }
  }
  p.value = min(c(1,2*(1 - b/B)))
  return(p.value)
}




