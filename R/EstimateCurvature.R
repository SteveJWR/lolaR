




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
#' @examples gEF(0,1,1,1)
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
#' @param min.curvature Minimum possible value of curvature
#' @param init.gridsize Gridsize to pick for good initialization
#'
#' @return Estimated Curvature Value
#' @export
#'
#' @examples EstimateKappa(1,1,1,sqrt(3/4))
EstimateKappa <- function(dxy,dxz,dyz,dxm,
                          max.iter = 10, ee.thresh = 0.1, abs.ee.thresh = 10**(-9),
                          min.curvature = -5000, init.gridsize = 10000){

  # first check for whether the midpoint estimate is already too far for any curvature:

  if(dxm < min(dxy,dxz) - (1/2)*dyz){
    warning("Triangle Inequality is not satisfied")
    return(-Inf)
  }


  # Picking good initialization for Newton Method
  max.curvature = (pi/max(c(dxy,dxz,dyz)))**2

  kap.seq <- seq(min.curvature,max.curvature, length.out = init.gridsize)
  g.vec <- gEF(kap.seq,dxy,dxz,dyz,dxm)
  kappa.init <- kap.seq[which.min(abs(g.vec))]


  suppressWarnings(d.xm.vec <- MidDist(kap.seq,dxy,dxz,dyz))

  kap.seq <- seq(min.curvature,0, 1)
  d.xm.vec <- c()
  for(kap in kap.seq){
    suppressWarnings(d.xm.vec <- c(d.xm.vec, MidDist(kap,dxy,dxz,dyz)))
  }
  min.dxm = min(d.xm.vec, na.rm = T)
  max.dxm = max(d.xm.vec, na.rm = T)

  # cases for impossible midpoint distances
  if(any(is.na(c(dxy,dxz,dyz)))){
    return(NA)
  } else if(dxm > max.dxm){
    return(Inf)
  } else if(dxm < min.dxm){
    return(-min.curvature)
  } else {

    kappa.prev = kappa.init
    rel.change = Inf
    iter = 0
    while(gEF(kappa.prev,dxy,dxz,dyz,dxm) < abs.ee.thresh|(rel.change > ee.thresh & iter < max.iter & kappa.prev < max.curvature & kappa.prev > min.curvature)){
      kappa.next = kappa.prev - gEF(kappa.prev,dxy,dxz,dyz,dxm)/gGradKappa(kappa.prev,dxy,dxz,dyz,dxm)
      rel.change = abs(1 - gEF(kappa.next,dxy,dxz,dyz,dxm)/(gEF(kappa.prev,dxy,dxz,dyz,dxm)))

      kappa.old = kappa.prev
      kappa.prev = kappa.next

      iter = iter + 1
    }
  }
  kappa.est = kappa.prev
  return(kappa.est)
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
optimal_midpoint_search <- function(D,top.k = 1, d.yz.min = 1.0,d.yz.max = 2.5){
  K <- nrow(D)
  idx.set <- expand.grid(1:K,1:K,1:K)
  idx.set <- idx.set[idx.set[,1] < idx.set[,2] & idx.set[,2] != idx.set[,3] & idx.set[,1] != idx.set[,3],]
  n.idx <- nrow(idx.set)
  obj.val <- sapply(1:n.idx, function(i){
    y = idx.set[i,1]
    z = idx.set[i,2]
    m = idx.set[i,3]
    dyz <- D[y,z]
    dym <- D[y,m]
    dzm <- D[z,m]

    dist.filter <- ifelse(dzm*dym*dyz == Inf, NaN, 1)

    obj.tmp = midpoint_objective(dzm,dym,dyz)*dist.filter
    out <- abs(obj.tmp - 1/2) + abs(dym - dzm)/dyz + 10000*(dyz > d.yz.max) +  10000*(dyz < d.yz.min)
    return(out)
  })

  sort.obj <- obj.val[order(obj.val)]
  ranked.idx <- idx.set[order(obj.val),]
  selected.set <- c()
  best.pairs <- matrix(NA, nrow = 0, ncol = 6)

  i = 1
  # stops the loop if we do not have a full top k which don't overlap
  while(nrow(best.pairs) < top.k & i <= nrow(ranked.idx)){
    if(!(any(as.vector(ranked.idx[i,]) %in% selected.set))){
      D.tmp <- D[as.numeric(ranked.idx[i,]),as.numeric(ranked.idx[i,])]
      best.pairs <- rbind(best.pairs, c(as.numeric(ranked.idx[i,]), D.tmp[1,3],D.tmp[2,3], D.tmp[1,2]))
      selected.set <- c(selected.set, c(ranked.idx[i,1],ranked.idx[i,2],ranked.idx[i,3]))
    }

    i = i+1
  }
  colnames(best.pairs) <- c("y","z","m","dym", "dzm", "dyz")
  rownames(best.pairs) <- NULL
  return(best.pairs)
}



#' Estimate Curvature for an input Distance Matrix
#'
#' @param D.hat Estimated Distance Matrix
#' @param tri.const Triangle Constant to use for selecting good x points
#' @param num.midpoints Number of midpoint sets to use
#' @param d.yz.min Minimum distance of y z pairs to consider
#' @param d.yz.max Maximum distance of y z pairs to consider
#' @param verbose Additional plotting details
#'
#' @return Return estimated values of curvature
#' @export
#'
#'
EstimateCurvature <- function(D.hat,
                              tri.const = 1.4,
                              num.midpoints = 3,
                              d.yz.min = 1.5,
                              d.yz.max = Inf,
                              verbose = T){


  mid.search <- optimal_midpoint_search(D.hat,top.k = 10,
                                        d.yz.min = min(d.yz.min, max(D.hat)/2),
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
    x.set <- filter_indices_2(D.hat,
                              y.opt,
                              z.opt,
                              m.opt,
                              tri.const = tri.const)

    kappa.set <- EstimateKappa_set(D.hat,
                                    y.opt,
                                    z.opt,
                                    m.opt,
                                    x.set)

    out.set <- list("kappas" = kappa.set,
                    "kappa.med" = median(kappa.set, na.rm = T),
                    "D" = D.hat,
                    "midpoints" = mid.search)

  } else {
    out.set <- list("kappas" = NULL,
                    "kappa.med" = NULL,
                    "D" = D.hat,
                    "midpoints" = mid.search)
  }

  return(out.set)
}



ConstantCurvatureTest <- function(D.hat, num.points = 3, tri.const = 1.4){

  mid.search <- optimal_midpoint_search(D.hat,top.k = num.points,
                                        d.yz.min = min(d.yz.min, max(D.hat)/2),
                                        d.yz.max = d.yz.max)


}





