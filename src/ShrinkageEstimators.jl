# This package is based on the shrinkage estimator from
# https://cran.r-project.org/web/packages/entropy/

module ShrinkageEstimators

include("ExportedFunctions.jl")
include("Discretization.jl")
include("Shrinkage.jl")
include("Formulae.jl")

end # module
