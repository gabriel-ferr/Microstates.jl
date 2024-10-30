#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#
#       GitHub: https://github.com/gabriel-ferr/Microstates.jl
#
#   References:
#       

module Microstates
#
#       Import the libraries that we need...
using Distances
using StatsBase
#       Define the float type that Microstates use.
__FLOAT_TYPE = Float64
#       Change the float type between Float32 or Float64.
function change_float_type(type::DataType)
    if (type != Float64 && type != Float32)
        println(string("The data type `", type, "` is not valid. Please use `Float64` or `Float32`."))
        return
    end
    __FLOAT_TYPE = type
end
#
#       Include the project files...
include("cpu/compute.jl")
include("utils/std_recurrence.jl")
#       Export the functions...
export microstates
#
end
