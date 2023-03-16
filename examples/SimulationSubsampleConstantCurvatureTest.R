library(lolaR)
library(Matrix)
library(CVXR)
library(Rmosek)
library(MASS)

#setwd("/Users/Owner/Documents/PhD/Networks/local_geometry_R/local_network_geometry")



### A Null example

source("00_functions.R")
source("clique_finder.R")
library(lolaR)
rm(filter_indices)
##### Simulated Datasets


scale.idx <- 2
kappa.idx <- 2


n.sims = 10
sim.idx <- 1:n.sims


kappa.set <- c(-2,-1,-0.5,0,0.5,1)
scale.set <- c(1/sqrt(2),1,2,4)


mu = -3
sd = 3

# midpoint objective constant > 0.5
#Cm <- 0.75

sim.avg.variance <- 0.25**2
p = 3
num.midpoints = 3

res = 1

kappa = kappa.set[kappa.idx]
if(kappa < 0){
  centers.radius = 2.5 # 2
} else {
  centers.radius = 2.5
}

centers.variance = 0.5**2

kappa.ests.results <- matrix(NA,nrow = n.sims*length(tri.const.seq),
                             ncol = length(scale.set) + 1)
sl.kappa.est.results <- matrix(NA,nrow = n.sims,
                               ncol = length(scale.set) + 1)
p.val.results <- matrix(NA, nrow = n.sims*length(tri.const.seq),
                        ncol = length(scale.set) + 1)
normalized.p.val.results <- p.val.results




scale <- scale.set[scale.idx]
n <- round(5000*scale)

n.centers <- round(100*sqrt(scale))
ell = round(8 + 4*log2(scale)) # min clique-size
approximate.variance <- sim.avg.variance

d.yz.min = 1.5
if(ell < 8){
  d.yz.min = 1
}

PI <- as.numeric(rdirichlet(1, rep(2,n.centers)))

cluster.model.variance = rgamma(n.centers, shape = approximate.variance)

lpcm <- latent_position_cluster_model(n,n.centers, p,
                                      centers.radius,
                                      kappa,
                                      cluster.model.variance,
                                      PI = PI)
# lpcm <- latent_position_cluster_model_2(n,n.centers, p, kappa,
#                                         centers.variance =centers.variance,
#                                         cluster.variance = approximate.variance,
#                                         PI = PI)
Z <- lpcm$Z

nu.vec <- rnorm(n,mean = mu*scale, sd = sd*scale)
nu.vec <- nu.vec*(nu.vec < 0)
nu.vec <- (exp(nu.vec) < 2/sqrt(n))*log(2/sqrt(n)) + (exp(nu.vec) >= 2/sqrt(n))*nu.vec


A.sim <- sim_ls_network_fast_2(nu.vec, Z, kappa)
d.sim <- colSums(A.sim)

# if(sim == 1 & plot.graph = T){
#   A.noiso <- A.sim[d.sim > 0, d.sim > 0]
#   #A.noiso <- A
#   A.noiso <- A.noiso[1:2000,1:2000]
#   d.noiso <- colSums(A.noiso)
#   A.noiso <- A.noiso[d.noiso > 0, d.noiso > 0 ]
#   g <- igraph::graph_from_adjacency_matrix(A.noiso, mode = "undirected")
#   V(g)$labels = NA
#   plot(g, vertex.size= 2,vertex.label=NA)
#
# }


# print(paste("Max cliques size:",max.clique.size))


# making sure at least ~ 50 cliques are found.
# when scale > 8 we have to use an approximate clique search
clique.set <- guided_clique_set(A.sim,lpcm$cluster_labels,
                                min_clique_size = ell)
clique.set <- clique_split(clique.set, min_clique_size = ell)

print(paste("Number of Cliques of size,",ell,":", length(clique.set)))

if(length(clique.set) > 40 ){
  clique.set <- clique.set[1:60]
}

J = 3

D.hat <- EstimateD(A.sim, clique.set, verbose = T)
reference.set <- selectReference(D.hat, J = J, tri.const = 1.4)
cc.test <- SubSampleConstantCurvatureTest(A.sim, clique.set, reference.set)


# Number of midpoints

tri.const <- 1.4
midpoints <- estimates$midpoints


lower.bounds <- c()
upper.bounds <- c()

lower.bounds.med <- c()
upper.bounds.med <- c()

locations <- c()
for(j in seq(J)){
  y = midpoints[j,1]
  z = midpoints[j,2]
  m = midpoints[j,3]

  x.set <- filter_indices_2(D.hat, y,
                            z,m,
                            tri.const = tri.const)
  bounds <- estimateBounds(D.hat,y,z,m,x.set)

  lower.bounds <- c(lower.bounds, bounds$lower.bounds)
  upper.bounds <- c(upper.bounds, bounds$upper.bounds)

  lower.bounds.med <- c(lower.bounds.med, median(bounds$lower.bounds))
  upper.bounds.med <- c(upper.bounds.med, median(bounds$upper.bounds))

  locations <- c(locations, rep(j, length(bounds$upper.bounds)))
}

idx <- seq(length(locations))

dat <- data.frame("idx" = idx,
                  "upper" = upper.bounds,
                  "lower" = lower.bounds,
                  "locations" = locations)
dat$locations <- as.factor(dat$locations)

dat.med <- data.frame("idx" = seq(J),
                      "upper" = upper.bounds.med,
                      "lower" = lower.bounds.med,
                      "locations" = seq(J))
dat.med$locations <- as.factor(dat.med$locations)

ggplot(data = dat, aes(x = idx, y = upper, shape = locations, col = "upper")) +
  geom_errorbar(data = dat, aes(x = idx, ymin = lower, ymax = upper, col = locations)) +
  ylim(-10,3) +
  geom_hline(yintercept = kappa) +
  ggtitle("Simulated Graph Constant Curvature K = 0.5, l = 12") +
  ylab("Curvature Estimate")

ggplot(data = dat.med, aes(x = idx, y = upper, shape = locations)) +
  geom_errorbar(data = dat.med, aes(x = idx, ymin = lower, ymax = upper, col = locations)) +
  ylim(-20,3) +
  geom_hline(yintercept = kappa) +
  ggtitle("Simulated Graph Constant Curvature K = 0.5, l = 6") +
  ylab("Curvature Estimate")

ggplot(data = dat.med, aes(x = idx, y = upper - lower)) +
  geom_line() +
  ylim(0,3) +
  ggtitle("Bias of Estimates Plot") +
  ylab("Curvature Estimate")









