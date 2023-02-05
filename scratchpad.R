EstimateKappa(1,1,1,sqrt(3/4))


kappa = seq(-0.000001,0.000001,length.out = 10000)
kappa = seq(-0.001,0.001,length.out = 10000)
kappa = seq(-0.1,0.1,length.out = 10000)



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







y = y.opt
z = z.opt
m = m.opt


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
  print(i)
}



