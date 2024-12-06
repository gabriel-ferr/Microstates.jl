#
#           Microstates.jl DEV Version
#           https://github.com/gabriel-ferr/Microstates.jl
#
#       By Gabriel Vinicius Ferreira, Sérgio Roberto Lopes and Thiago de Lima Prado
#
#       References:
#   [1] G. Corso, T. De Lima Prado, G. Z. D. S. Lima, J. Kurths, and S. R. Lopes, 
#   Quantifying Entropy Using Recurrence Matrix Microstates, Chaos an Interdisciplinary Journal of Nonlinear Science 28, (2018).
#
#
#       Contact: gabriel.vferreira@icloud.com
#
# ===========================================================================================
#           To use this dev version:
#       julia> Pkg.add("https://github.com/gabriel-ferr/Microstates.jl")
#       ...
#       julia> using Microstates
# ===========================================================================================
module Microstates
# -------------------------------------------------------------------------------------------
using Distances
using StatsBase
using Statistics
# -------------------------------------------------------------------------------------------
include("cpu/compute.jl")
include("utils/entropy.jl")
include("utils/std_recurrence.jl")
include("utils/crd_recurrence.jl")
# -------------------------------------------------------------------------------------------
export power_vector
export microstates
# ===========================================================================================
end