#
#               By Gabriel Ferreira
#                   Orientation: Prof. Dr. Thiago de Lima Prado
#                                Prof. Dr. Sergio Roberto Lopes
#
# =================================================================================================
#           I don't know how to compile a library with Julia and add it to the PKG path, then I
#   did't do it... so, to use this library you need to include it in your project and
#   say that you want to use it, like:
#       julia> include("microstates.jl")
#       julia> using .Microstates
#       
#           You will also need some libraries, so remember to import them into your PKG
#       repository before using this code =3
#           Libraries you may have a need for:
#           ∙ StatsBase.jl
#           ∙ Distances.jl
#
#       Obs: I think that you could try "import Pkg; Pkg.add("Microstates")" to install it in Pkg path ...
#
#       
module Microstates

#       This variable defines the float type used in the system. 
#       You can use it to select the data type returned by the functions.
#           !! Only Float64 and Float32 !!
__FLOAT_TYPE = Float64

using StatsBase
using Distances
using CUDA

include("src/utils.jl")
include("src/entropy.jl")
include("src/max_entroy.jl")
include("src/compute_cpu.jl")
include("src/compute_gpu.jl")

#       Changes the float format =D
function change_float_type(type::DataType)
    if (type != Float64 && type != Float32)
        println("You can only change the data format to Float32 or Float64!")
        return
    end
    global __FLOAT_TYPE = type
    println(string("Done: ", __FLOAT_TYPE))
end

export change_float_type
export entropy
export gpu_microstates
export microstates
export power_vector

end