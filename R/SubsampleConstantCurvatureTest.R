

#Upper bound on curvature Estimate
g_u <- function(kappa, dxy, dxz, dyz, dxm, dym, dzm){
  out <- gEF(kappa,dxy,dxz,dyz,dxm + MidDist(kappa,dym,dzm,dyz))
}

g_l <- function(kappa, dxy, dxz, dyz, dxm, dym, dzm){
  out <- gEF(kappa,dxy,dxz,dyz,dxm - MidDist(kappa,dym,dzm,dyz))
}




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

#D = D.hat.sub # , y,z,m,x.set
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
selectReference <- function(D,J,tri.const = 1.4){
  midpoints <- optimal_midpoint_search(D, top.k = J, d.yz.min = 1, d.yz.max = 5)
  out.list <- list()
  for(j in seq(J)){
    y <- midpoints[j,1]
    z <- midpoints[j,2]
    m <- midpoints[j,3]
    x.set <- filter_indices(D,y,z,m, tri.const = tri.const)
    out.list[[j]] <- list("y" = y, "z" = z, "m" = m, "xset" = x.set)

  }
  return(out.list)
}

### Subsample

SubSampleConstantCurvatureTest <- function(A,clique.set,reference.set,
                                           subsample.rate = 1,B = 100){
  K = length(clique.set)
  J = length(reference.set)
  upper.bound.sub <- rep(NA, B)
  lower.bound.sub <- rep(NA, B)
  clique.subsample <- list()
  D.hat <- EstimateD(A,clique.set)
  for(b in seq(B)){
    for(k in seq(K)){
      ell = length(clique.set[[k]])
      clique.subsample[[k]] <- sample(clique.set[[k]],
                                      size = ell - subsample.rate,
                                      replace = F)
    }
    sub.idx <- CliquesCount(clique.subsample)

    D.hat.sub <- EstimateD(A,clique.subsample, D0 = D.hat, verbose = F)
    upper.min <- Inf
    lower.max <- -Inf
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
  return(list("p.value" = p.value, "upper.bounds" = upper.bound.boot, "lower.bounds" = lower.bound.boot))
}


# Overlapping Bounds p-value
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

