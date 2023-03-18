library(devtools)

use_github()
create_package("/Users/Owner/Documents/PhD/Rpackages/lolaR")

use_r("BasicMath")
use_r("FindCliques")
use_r("EstimateDistances")
use_r("EstimateCurvature")
use_r("SubsampleConstantCurvatureTest.Rurvature")

#list of required packages
# igraph, CVXR, ggplot2, Matrix

use_package("igraph")
use_package("CVXR")
use_package("Rmosek", "Suggests")
use_package("ggplot2")
use_package("Matrix")


use_testthat()
use_test("BasicMath")
use_test("FindCliques")
use_test("EstimateDistances")
use_test("EstimateCurvature")

test()

export(BasicMath)
export(FindCliques)
export(EstimateDistances)
export(EstimateCurvature)

document()

check()


### Datasets

my_pkg_data <- sample(1000)
usethis::use_data(my_pkg_data)

usethis::use_data()


load_all()

install()

# Add a Readme
use_readme_rmd()

## Add a Vignette
usethis::use_vignette("example-curvature-estimation")





