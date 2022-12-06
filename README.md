
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lolaR

<!-- badges: start -->
<!-- badges: end -->

lolaR is a package built to estimate the curvature in latent distance
network models as described in Hoff, Raftery, and Handcock
([2002](#ref-Hoff2002latent)) as well as test whether a network could
have been generated by a latent space of constant curvature. We estimate
the latent distances on a subset of the network using cliques. We then
exploit the fact that triangles allow for identification of the latent
space provided we have a midpoint. See Wilkins-Reeves and McCormick
([2022](#ref-wilkins2022asymptotically)) for additional details

## Installation

You can install the development version of lolaR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SteveJWR/lolaR")
```

## Example

We highlight this application with two collaboration networks in Physics
article co-authorship from Leskovec, Kleinberg, and Faloutsos
([2007](#ref-Leskovec2007GraphEvolution)) and available at Leskovec and
Krevl ([2014](#ref-snapnets)). We consider the problem of testing
whether a latent space model has constant curvature. To do so we first
consider two collaboration networks, one for Astrophysics and one for
Condensed Matter Physics.

``` r
library(lolaR)
library(ggplot2)
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(Matrix)

load("data/Astro_adjacency_matrix.rda")
load("data/CondMat_adjacency_matrix.rda")

# Adjacency Matrix Names
#G.astro
#G.cm 
```

We next plot the largest connected components from the first 2000
indices.

``` r

g.sub.astro <- igraph::graph_from_adjacency_matrix(G.astro[1:2000,1:2000], mode = "undirected")
g.sub.cm <- igraph::graph_from_adjacency_matrix(G.cm[1:2000,1:2000], mode = "undirected")


connected.idx.astro  <- which(components(g.sub.astro, mode = "strong")$membership == 4)
connected.idx.cm  <- which(components(g.sub.cm, mode = "strong")$membership == 35)


g.sub.astro <- igraph::graph_from_adjacency_matrix(G.astro[connected.idx.astro,connected.idx.astro], mode = "undirected")
g.sub.cm <- igraph::graph_from_adjacency_matrix(G.cm[connected.idx.cm,connected.idx.cm], mode = "undirected")

V(g.sub.astro)$name = NA
V(g.sub.cm)$name = NA
```

``` r
plot(g.sub.astro, vertex.size=4)
```

<img src="man/figures/README-Example Plot 1-1.png" width="100%" />

``` r
plot(g.sub.cm, vertex.size=4)
```

<img src="man/figures/README-Example Plot 2-1.png" width="100%" />

We first find a list of cliques of size $19$ and $12$ respectively. This
is in order to have a moderately large ($40-60$) number of rows/columns
in our distance matrix.

We can either conduct an approximate search for cliques at least of size
($\ell = \{19,12\}$ respectively)

``` r
cliques.sample.astro = ClusterCliqueSearch(G.astro, min_clique_size = 19)
cliques.sample.cm = ClusterCliqueSearch(G.cm, min_clique_size = 12)
```

But since these networks are relatively sparse, we can compute the exact
set of maximal cliques and partition them into non-overlapping sets.

``` r
g.astro <- igraph::graph_from_adjacency_matrix(G.astro, mode = "undirected")
g.cm <- igraph::graph_from_adjacency_matrix(G.cm, mode = "undirected")


cliques.astro = igraph::maximal.cliques(g.astro, min = 19)
cliques.cm = igraph::maximal.cliques(g.cm, min = 12)

cliques.astro <- CliquePartition(cliques.astro)
cliques.cm <- CliquePartition(cliques.cm)
```

Next using the set of cliques, we can estimate a distance matrix under
the latent distance model without first choosing an embedding space.

``` r

D.astro <- EstimateD(G.astro, cliques.astro, verbose = T)
#> Num Steps: 1 Likelihood Stopping Criteria: 1 Num Steps: 2 Likelihood Stopping Criteria: 0.0125474 Num Steps: 3 Likelihood Stopping Criteria: 0.0024412 Num Steps: 4 Likelihood Stopping Criteria: 0.0002361 Num Steps: 5 Likelihood Stopping Criteria: 3.9e-05 Num Steps: 6 Likelihood Stopping Criteria: 4.63e-05 Num Steps: 7 Likelihood Stopping Criteria: 2.32e-05 Num Steps: 8 Likelihood Stopping Criteria: 1.1e-06 Num Steps: 9 Likelihood Stopping Criteria: 0 
D.cliques <- EstimateD(G.cm, cliques.cm, verbose = T)
#> Num Steps: 1 Likelihood Stopping Criteria: 1 Num Steps: 2 Likelihood Stopping Criteria: 0.0002213 Num Steps: 3 Likelihood Stopping Criteria: 3.2e-06 Num Steps: 4 Likelihood Stopping Criteria: 0 
```

We first can search for the best midpoint and estimate the latent
curvature for each model.

``` r
kappa.astro <- EstimateCurvature(D.astro, verbose = T, d.yz.min = 1, d.yz.max = 4.5)
#> [1] "Midpoints: "
#>       y  z  m      dym      dzm      dyz
#> [1,] 25 38 35 2.038540 2.038771 4.077311
#> [2,] 21 55 39 1.716070 1.718500 3.434570
#> [3,]  6 16 29 1.560494 1.580823 3.141317
kappa.cm <- EstimateCurvature(D.cm, verbose = T, d.yz.min = 1, d.yz.max = 4.5)
#> [1] "Midpoints: "
#>       y  z  m      dym      dzm      dyz
#> [1,] 22 33 36 2.041491 2.056636 4.098127
#> [2,] 24 30 20 2.182764 2.160485 4.343248
#> [3,]  8 27 34 1.778742 1.725676 3.504418
```

From these distance matrices we are able to search for indices which
approximately form a set \$ y,z,m \$ where \$ m \$ is the midpoint of \$
y \$ and \$ z $, i.e. ($ d\_{yz} = 2d\_{ym} = 2d\_{zm} \$). After
finding a set of \$ R \$ non-overlapping indices such approximately
satisfying this midpoint equation \$
{y<sup>{(r)},z</sup>{(r)},m^{(r)}}\_{r = 1}^R \$.

We would like to test the hypothesis
$$ H_0: \kappa(r) = \kappa \text{ for all } r \in \{1,2,\dots, R\} $$
Where \$ (r) \$ corresponds to the curvature at the corresponding point
\$ r \$.

``` r

test.astro <- ConstantCurvatureTest(D.astro, num.midpoints = 3, d.yz.min = 1, d.yz.max = 4.5)
test.cm <- ConstantCurvatureTest(D.cm, num.midpoints = 3, d.yz.min = 1, d.yz.max = 4.5)
```

``` r
med.vec.astro <- rep(0,3)
med.vec.cm <- rep(0,3)

for(k in seq(3)){
  med.vec.astro[k] = SoftThreshold(median(test.astro$estimates[test.astro$loc == k]), 10)
  med.vec.cm[k] = SoftThreshold(median(test.cm$estimates[test.cm$loc == k]),10)
}
```

``` r
library(ggplot2)


test.astro$scaled_estimates = SoftThreshold(test.astro$estimates, 10)
dat.astro <- data.frame(matrix(c(test.astro$loc, test.astro$scaled_estimates), ncol = 2))

names(dat.astro) = c("loc", "scaled_estimates")
ggplot(dat.astro, aes(y = scaled_estimates, x = loc)) +
  geom_jitter() +
  geom_vline(xintercept = 0.5, col = "red") +
  geom_vline(xintercept = 1.5, col = "red") +
  geom_vline(xintercept = 2.5, col = "red") +
  geom_vline(xintercept = 3.5, col = "red") +
  geom_segment(aes(x=0.5,xend=1.5,y=med.vec.astro[1],yend=med.vec.astro[1]), col = "blue") +
  geom_segment(aes(x=1.5,xend=2.5,y=med.vec.astro[2],yend=med.vec.astro[2]), col = "blue") +
  geom_segment(aes(x=2.5,xend=3.5,y=med.vec.astro[3],yend=med.vec.astro[3]), col = "blue") +
  ggtitle("Curvature Estimates Within Astrophysics Network") +
  ylab("Trimmed Curvature") +
  xlab("Midpoint Set")
```

<img src="man/figures/README-Plotting Curvature Estimated In Each Region Astrophysics-1.png" width="100%" />

``` r



test.cm$scaled_estimates = SoftThreshold(test.cm$estimates, 10)
dat.cm <- data.frame(matrix(c(test.cm$loc, test.cm$scaled_estimates), ncol = 2))

names(dat.cm) = c("loc", "scaled_estimates")
ggplot(dat.cm, aes(y = scaled_estimates, x = loc)) +
  geom_jitter() +
  geom_vline(xintercept = 0.5, col = "red") +
  geom_vline(xintercept = 1.5, col = "red") +
  geom_vline(xintercept = 2.5, col = "red") +
  geom_vline(xintercept = 3.5, col = "red") +
  geom_segment(aes(x=0.5,xend=1.5,y=med.vec.cm[1],yend=med.vec.cm[1]), col = "blue") +
  geom_segment(aes(x=1.5,xend=2.5,y=med.vec.cm[2],yend=med.vec.cm[2]), col = "blue") +
  geom_segment(aes(x=2.5,xend=3.5,y=med.vec.cm[3],yend=med.vec.cm[3]), col = "blue") +
  ggtitle("Curvature Estimates Within CMP Network") +
  ylab("Trimmed Curvature") +
  xlab("Midpoint Set")
```

<img src="man/figures/README-Plotting Curvature Estimated In Each Region CMP-1.png" width="100%" />

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Hoff2002latent" class="csl-entry">

Hoff, Peter D, Adrian E Raftery, and Mark S Handcock. 2002. “Latent
Space Approaches to Social Network Analysis.” *Journal of the American
Statistical Association* 97 (460): 1090–98.

</div>

<div id="ref-Leskovec2007GraphEvolution" class="csl-entry">

Leskovec, Jure, Jon Kleinberg, and Christos Faloutsos. 2007. “<span
class="nocase">Graph evolution</span>.” *ACM Transactions on Knowledge
Discovery from Data (TKDD)* 1 (1).
<https://doi.org/10.1145/1217299.1217301>.

</div>

<div id="ref-snapnets" class="csl-entry">

Leskovec, Jure, and Andrej Krevl. 2014. “SNAP Datasets: Stanford Large
Network Dataset Collection.” <http://snap.stanford.edu/data>.

</div>

<div id="ref-wilkins2022asymptotically" class="csl-entry">

Wilkins-Reeves, Steven, and Tyler McCormick. 2022. “Asymptotically
Normal Estimation of Local Latent Network Curvature.” *arXiv Preprint
arXiv:2211.11673*.

</div>

</div>
