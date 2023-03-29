




# estimating function for kappa
#' Estimating Function Value
#'
#' @param kappa Curvature of Embedding Space
#' @param dxy Triangle Side Length
#' @param dxz Triangle Side Length
#' @param dyz Triangle Side Length
#' @param dxm Length from point x to the midpoint of y and z
#'
#' @return Value of estimating function. This will equal 0 when triangle can be embedded in the space of curvature kappa
#' @export
#'
#' @examples gEF(0,1,1,1,sqrt(3/4))
gEF <- function(kappa,dxy,dxz,dyz,dxm){
  # using a single value to construct the estimating function use 0i to compute using complex numbers.
  thresh = 10**(-9)
  threshi = thresh*sqrt(-1 + 0i)
  # Terms added for vectorization and numerical stability
  g = (abs(kappa)  > 2*thresh)*((2*cos(dxm*sqrt(kappa + thresh + threshi)) - sec(dyz/2*sqrt(kappa + thresh + threshi))*(cos(dxy*sqrt(kappa + thresh + threshi)) + cos(dxz*sqrt(kappa + thresh + threshi))))/(kappa + thresh + threshi))+ (abs(kappa)  <= 2*thresh)*(((dxy**2 + dxz**2)/2) - dyz**2/4 - dxm**2)
  g = Re(g)
  return(g)
}


# true distance to midpoint as a function of kappa, and triangle side lengths
#' Value of dxm as a function of the triangle side lengths and the curvature
#'
#' @param kappa Curvature of Embedding Space
#' @param dxy Triangle Side Length
#' @param dxz Triangle Side Length
#' @param dyz Triangle Side Length
#'
#' @return Distance of the length xm under the specified embedding space
#' @export
#'
#' @examples MidDist(0,1,1,1)
MidDist <- function(kappa,dxy,dxz,dyz){
  # numerical stability term
  thresh = 10^(-9)
  threshi = thresh*sqrt(-1 + 0i)
  out.com <- (kappa != 0)*((1/(sqrt(kappa + thresh*sign(kappa) + threshi)))*acos((1/2)*sec(dyz/2*sqrt(kappa + thresh*sign(kappa) + threshi))*(cos(dxy*sqrt(kappa + thresh*sign(kappa) + threshi)) + cos(dxz*sqrt(kappa + thresh*sign(kappa) + threshi))))) + (kappa == 0)*(sqrt((dxy**2 + dxz**2)/2 - dyz**2/4))
  out <- Re(out.com)
  return(out)
}


# derivative of the estimating function.
#' Derivative of implicit estimating function
#'
#' @param kappa Curvature of Embedding Space
#' @param dxy Triangle Side Length
#' @param dxz Triangle Side Length
#' @param dyz Triangle Side Length
#' @param dxm Length from point x to the midpoint of y and z
#'
#' @return Gradient of the estimating function
#' @export
#'
#' @examples gGradKappa(0,1,1,1,1/2)
gGradKappa <- function(kappa,dxy,dxz,dyz,dxm){
  thresh = 10**(-9)
  threshi = thresh*sqrt(-1 + 0i)
  # log-scale computing for small kappa
  if(kappa != 0){
    gp = (-1)*(1/(4*(kappa + thresh*sign(kappa) + threshi)^2))*(8*cos(dxm*sqrt(kappa + thresh*sign(kappa) + threshi)) - 4*cos(dxz*sqrt(kappa + thresh*sign(kappa) + threshi))*sec((dyz*sqrt(kappa + thresh*sign(kappa) + threshi))/2) +
                                 4*dxm*sqrt(kappa + thresh*sign(kappa) + threshi)*sin(dxm*sqrt(kappa + thresh*sign(kappa) + threshi)) - 2*dxy*sqrt(kappa + thresh*sign(kappa) + threshi)*sec((dyz*sqrt(kappa + thresh*sign(kappa) + threshi))/2)*sin(dxy*sqrt(kappa + thresh*sign(kappa) + threshi)) -
                                 2*dxz*sqrt(kappa + thresh*sign(kappa) + threshi)*sec((dyz*sqrt(kappa + thresh*sign(kappa) + threshi))/2)*sin(dxz*sqrt(kappa + thresh*sign(kappa) + threshi)) +
                                 dyz*sqrt(kappa + thresh*sign(kappa) + threshi)*cos(dxz*sqrt(kappa + thresh*sign(kappa) + threshi))*sec((dyz*sqrt(kappa + thresh*sign(kappa) + threshi))/2)*tan((dyz*sqrt(kappa + thresh*sign(kappa) + threshi))/2) +
                                 cos(dxy*sqrt(kappa + thresh*sign(kappa) + threshi))*sec(dyz*sqrt(kappa + thresh*sign(kappa) + threshi)/2)*(-4 + dyz*sqrt(kappa + thresh*sign(kappa) + threshi)*tan(dyz*sqrt(kappa + thresh*sign(kappa) + threshi)/2)))
    gp = Re(gp)
  } else if(kappa == 0){
    gp = (-1)*(1/192)*(-16*dxm^4 + 8*dxy^4 + 8*dxz^4 - 12*dxy^2*dyz^2 - 12*dxz^2*dyz^2 + 5*dyz^4)
  }
  return(gp)
}

# gradient with respect to d
# gGradD <- function(kappa, dxy, dxz, dyz, dxm){
#   g1 <- (-1/(sqrt(kappa + 0i)))*sec(dyz/2*sqrt(kappa + 0i))*sin(dxy/2*sqrt(kappa + 0i))
#   g2 <- (-1/(sqrt(kappa + 0i)))*sec(dyz/2*sqrt(kappa + 0i))*sin(dxz/2*sqrt(kappa + 0i))
#   g3 <- (1/(2*sqrt(kappa + 0i)))*sec(dyz/2*sqrt(kappa + 0i))*tan(dyz/2*sqrt(kappa + 0i))*(cos(dxy*sqrt(kappa + 0i)) + cos(dxz*sqrt(kappa + 0i)))
#   g4 <- (2/(sqrt(kappa + 0i)))*sin(dxm*sqrt(kappa + 0i))
#   g.vec <- c(g1,g2,g3,g4)
#   return(g.vec)
# }



#' Find the solution to the implicit estimating function to estimate curvature
#'
#' @param dxy Triangle Side Length
#' @param dxz Triangle Side Length
#' @param dyz Triangle Side Length
#' @param dxm Length from point x to the midpoint of y and z
#' @param max.iter Maximum number of iterations
#' @param ee.thresh Threshold for the ratio of the solution of the estimating function
#' @param abs.ee.thresh Threshold for the absolute value of the solution of the estimating function
#' @param kappa.thresh Threshold for the gridsearch precision
#' @param min.curvature Minimum possible value of curvature
#' @param init.gridsize Gridsize to pick for good initialization
#' @param gridsearch Whether to use grid-search for estimating curvature as opposed to a Newton Method.
#'
#' @return Estimated Curvature Value
#' @export
#'
#' @examples EstimateKappa(1,1,1,sqrt(3/4))
EstimateKappa <- function(dxy,dxz,dyz,dxm,
                          max.iter = 10, ee.thresh = 0.001, abs.ee.thresh = 10**(-9), kappa.thresh = 0.0001,
                          min.curvature = -5000, init.gridsize = 10000, gridsearch = T){

  # first check for whether the midpoint estimate is already too far for any curvature:

  if(dxm < min(dxy,dxz) - (1/2)*dyz){
    warning("Triangle Inequality is not satisfied")
    return(-Inf)
  }


  # Picking good initialization for Newton Method
  max.curvature = (pi/max(c(dxy,dxz,dyz)))**2
  kappa.upper = max.curvature
  kappa.lower = min.curvature

  kap.seq <- c(seq(min.curvature,0, length.out = init.gridsize),
               seq(0,max.curvature - ee.thresh, length.out = init.gridsize))
  kap.seq <- kap.seq[-init.gridsize]
  suppressWarnings(d.xm.vec <- MidDist(kap.seq,dxy,dxz,dyz))

  min.dxm = min(d.xm.vec, na.rm = T)
  max.dxm = max(d.xm.vec, na.rm = T)
  if(any(is.na(c(dxy,dxz,dyz)))){
    warning("Input Distance Was Undefined")
    return(NA)
  } else if(dxm > max.dxm){
    warning("Midpoint Distance is too large for embedding geometries")
    return(max.curvature)
  } else if(dxm < min.dxm){
    warning("Midpoint Distance is too small for embedding geometries")
    return(min.curvature)
  }

  if(gridsearch){
    while(abs(kappa.upper - kappa.lower)  > kappa.thresh){
      best.idx = which.min(abs(gEF(kap.seq,dxy,dxz,dyz,dxm)))
      kappa.est = kap.seq[best.idx]
      kappa.upper = kap.seq[min(c(best.idx + 1, length(kap.seq)))]
      kappa.lower = kap.seq[max(c(best.idx - 1, 1))]
      kap.seq = seq(kappa.lower, kappa.upper, length.out = init.gridsize)
    }
    return(kappa.est)
  } else {
    g.vec <- gEF(kap.seq,dxy,dxz,dyz,dxm)

    #Find the first endpoint
    kappa.endpoint = which.max(g.vec > 0 & is.finite(g.vec))
    kappa.max <- kap.seq[kappa.endpoint]
    grid.size <- kap.seq[kappa.endpoint] - kap.seq[kappa.endpoint - 1]
    kappa.max <- kappa.max + 2*grid.size
    kappa.min <- kappa.max - 5*grid.size

    kap.seq <- seq(kappa.min,min(kappa.max, max.curvature - ee.thresh), length.out = init.gridsize)

    g.vec <- gEF(kap.seq,dxy,dxz,dyz,dxm)
    kappa.init <- kap.seq[which.min(abs(g.vec))]

    kappa.prev = kappa.init
    rel.change = Inf
    iter = 0
    kappa.seq <- c(kappa.prev)
    rel.change.seq <- c()

    while((abs(gEF(kappa.prev,dxy,dxz,dyz,dxm)) > abs.ee.thresh ) & (rel.change > ee.thresh & iter < max.iter & kappa.prev < max.curvature & kappa.prev > min.curvature)){
      kappa.next = kappa.prev - gEF(kappa.prev,dxy,dxz,dyz,dxm)/gGradKappa(kappa.prev,dxy,dxz,dyz,dxm)
      rel.change = abs(1 - abs(gEF(kappa.next,dxy,dxz,dyz,dxm))/abs(gEF(kappa.prev,dxy,dxz,dyz,dxm)))
      kappa.seq <- c(kappa.seq,kappa.next)
      rel.change.seq <- c(rel.change.seq, rel.change)

      kappa.old = kappa.prev
      kappa.prev = kappa.next

      iter = iter + 1
    }
    # cases for impossible midpoint distances

    kappa.est = kappa.prev
    return(kappa.est)
  }
}


#
#' Midpoint objective function
#'
#' @export
#'
midpoint_objective <- function(dym,dzm,dyz){
  out = (dym**2 + dzm**2)/(dyz**2)
  return(out)
}


#
#' Midpoint Search Function
#'
#' @export
#'
optimal_midpoint_search <- function(D,top.k = 1, d.yz.min = 1.0,d.yz.max = 2.5, subset){

  K <- nrow(D)
  if (missing(subset)) {
    subset = 1:K
  }
  obj.upper.bound = 4*max(D)^2/min(D[D > 0])^2

  idx.set <- expand.grid(subset, subset, subset)
  idx.set <- idx.set[idx.set[, 1] < idx.set[, 2] & idx.set[,2] != idx.set[, 3] & idx.set[, 1] != idx.set[, 3], ]
  n.idx <- nrow(idx.set)
  obj.val <- sapply(1:n.idx, function(i) {
    y = idx.set[i, 1]
    z = idx.set[i, 2]
    m = idx.set[i, 3]
    dyz <- D[y, z]
    dym <- D[y, m]
    dzm <- D[z, m]
    dist.filter <- ifelse(dzm * dym * dyz == Inf, NaN, 1)
    obj = midpoint_objective(dzm, dym, dyz) * dist.filter
    out <- abs(obj - 1/2) + abs(dym - dzm)/dyz + obj.upper.bound *
      (dyz > d.yz.max) + obj.upper.bound * (dyz < d.yz.min)
    return(out)
  })

  sort.obj <- obj.val[order(obj.val)]
  ranked.idx <- idx.set[order(obj.val), ]
  selected.set <- c()
  best.pairs <- matrix(NA, nrow = 0, ncol = 6)
  i = 1
  while (nrow(best.pairs) < top.k & i <= nrow(ranked.idx)) {
    if (!(any(as.vector(ranked.idx[i, ]) %in% selected.set))) {
      D.tmp <- D[as.numeric(ranked.idx[i, ]), as.numeric(ranked.idx[i,
      ])]
      best.pairs <- rbind(best.pairs, c(as.numeric(ranked.idx[i,
      ]), D.tmp[1, 3], D.tmp[2, 3], D.tmp[1, 2]))
      selected.set <- c(selected.set, c(ranked.idx[i, 1],
                                        ranked.idx[i, 2], ranked.idx[i, 3]))
    }
    i = i + 1
  }
  colnames(best.pairs) <- c("y", "z", "m", "dym", "dzm", "dyz")
  rownames(best.pairs) <- NULL
  return(best.pairs)
}





#' Triangle Inequality Filter For Distance Matrices
#'
#'
#' @param D Distance matrix
#' @param y Index of midpoint set
#' @param z Index of midpoint set
#' @param m Index of midpoint set
#' @param tri.const Scaling constant for triangle inequality
#'
#' @return Indices which satisfy the stronger triangle inequality
#' @export
#'
filter_indices <- function(D,y,z,m, tri.const = sqrt(2)){
  K = nrow(D)
  x.set <- 1:K
  x.set <- x.set[-c(y,z,m)]

  #indicator of whether to include an x
  x.id <- sapply(x.set, function(x){
    # ensuring that the distance is not infinite to each of x.y.m
    if(is.infinite(D[x,y]) | is.infinite(D[x,z]) | is.infinite(D[x,m])){
      return(F)
    } else {
      i1 = D[y,z] + D[x,z] >= tri.const*D[x,y]
      i2 = D[y,z] + D[x,y] >= tri.const*D[x,z]
      i3 = D[x,y] + D[x,z] >= tri.const*D[y,z]

      return(ifelse(i1*i2*i3 == 1, T, F))
    }

  })
  x.filtered <- x.set[x.id]
  return(x.filtered)
}



#' Estimate a set of curvatures using a set of indices by x
#'
#' @param D Distance matrix
#' @param y Index of midpoint set
#' @param z Index of midpoint set
#' @param m Index of midpoint set
#' @param x.set Vector of indices for which we want to estimate
#'
#' @return Vector of kappa estimates corresponding to x.set indices
#' @export
#'
EstimateKappaSet <- function(D,y,z,m,x.set){
  d1 = nrow(D)
  d2 = ncol(D)
  if(d1 > d2){
    D = t(D)
  }

  kappa.set <- rep(NA,length(x.set))

  for(i in seq(length(x.set))){
    x = x.set[i]
    dxy = D[y,x]
    dxz = D[z,x]
    dyz = D[y,z]
    dxm = D[m,x]
    d.vec = c(dxy,dxz,dyz,dxm)
    if(any(d.vec == Inf)){
      kappa.hat.x <- NA
    } else{
      kappa.hat.x <- EstimateKappa(dxy,dxz,dyz,dxm)
    }
    kappa.set[i] = kappa.hat.x
  }

  return(kappa.set)
}



#' Estimate Curvature for an input Distance Matrix
#'
#' @param D Estimated Distance matrix
#' @param tri.const Triangle Constant to use for selecting good x points
#' @param num.midpoints Number of midpoint sets to use
#' @param d.yz.min Minimum distance of y z pairs to consider
#' @param d.yz.max Maximum distance of y z pairs to consider
#' @param verbose Additional plotting details
#'
#' @return Return estimated values of curvature
#' @export
EstimateCurvature <- function(D,
                              tri.const = 1.4,
                              num.midpoints = 3,
                              d.yz.min = 1.5,
                              d.yz.max = Inf,
                              verbose = F){


  mid.search <- optimal_midpoint_search(D,top.k = 10,
                                        d.yz.min = min(d.yz.min, max(D)/2),
                                        d.yz.max = d.yz.max)
  if(verbose){
    print("Midpoints: ")
    max.row <- min(3,nrow(mid.search))
    print(mid.search[1:max.row,])
  }

  y.opt = mid.search[1,1]
  z.opt = mid.search[1,2]
  m.opt = mid.search[1,3]

  opt.vec <- c(y.opt,z.opt,m.opt)
  if(!any(is.na(opt.vec))){
    x.set <- filter_indices(D,
                              y.opt,
                              z.opt,
                              m.opt,
                              tri.const = tri.const)

    kappa.set <- EstimateKappaSet(D,
                                  y.opt,
                                  z.opt,
                                  m.opt,
                                  x.set)

    out.set <- list("kappas" = kappa.set,
                    "kappa.med" = median(kappa.set, na.rm = T),
                    "midpoints" = mid.search)

  } else {
    out.set <- list("kappas" = NULL,
                    "kappa.med" = NULL,
                    "midpoints" = mid.search)
  }
  return(out.set)
}




#' Test for Constant curvature of a distance matrix
#'
#' @param D Estimated Distance matrix
#' @param tri.const Triangle Constant to use for selecting good x points
#' @param num.midpoints Number of midpoint sets to use
#' @param d.yz.min Minimum distance of y z pairs to consider
#' @param d.yz.max Maximum distance of y z pairs to consider
#' @param curve.scale Threshold to use to minimize effect of outliers
#' @param verbose Additional plotting details
#'
#' @return Test for constant curvature of a given distance matrix
#' @export
#'
ConstantCurvatureTest <- function(D, tri.const = 1.4,
                                  num.midpoints = 3,
                                  d.yz.min = 1.5,
                                  d.yz.max = Inf,
                                  curve.scale = 10,
                                  verbose = F){

  mid.search <- optimal_midpoint_search(D,top.k = num.midpoints,
                                        d.yz.min = min(d.yz.min, max(D)/2),
                                        d.yz.max = d.yz.max)

  location.vec <- c()
  curvature.vec <- c()
  rescaled.curvature.vec <- c()
  K = num.midpoints
  for(k in seq(K)){
    y.opt = mid.search[k,1]
    z.opt = mid.search[k,2]
    m.opt = mid.search[k,3]
    opt.vec <- c(y.opt,z.opt,m.opt)
    if(!any(is.na(opt.vec))){
      x.set <- filter_indices(D,
                              y.opt,
                              z.opt,
                              m.opt,
                              tri.const = tri.const)

      kappa.set <- EstimateKappaSet(D,
                                    y.opt,
                                    z.opt,
                                    m.opt,
                                    x.set)
      location.vec <- c(location.vec, rep(k,length(kappa.set)))
      curvature.vec <- c(curvature.vec,kappa.set)
      y = SoftThreshold(kappa.set, curve.scale)
      # transformation scale
      rescaled.kappa <- (y - median(y))/mad(y) + median(y)
      rescaled.curvature.vec <- c(rescaled.curvature.vec, rescaled.kappa)
    }
  }

  if(length(unique(location.vec)) > 1 ){
    test.dat <- data.frame("loc" = location.vec, "est" = curvature.vec)
    trim.test.dat <- data.frame("loc" = location.vec, "est" = rescaled.curvature.vec)
    norm.test <- kruskal.test(est ~ loc, data = trim.test.dat)
    test <-  kruskal.test(est ~ loc, data = test.dat)
    out.list <- list("p.value" =  test$p.value,"norm.p.value" =norm.test$p.value, "loc" = location.vec, "estimates" = curvature.vec, "transformed_estimates"=rescaled.curvature.vec)
  } else {
    out.list <- list("p.value" =  NULL, "norm.p.value" =  NULL, "loc" = NULL,  "estimates" = NULL, "transformed_estimates"= NULL)
  }

  return(out.list)
}





