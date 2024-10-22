#
#       -- Microstates.jl: read more in https://github.com/gabriel-ferr/Microstates.jl
#

#       I'm going to make a function here to find the maximum entropy for standard recurrence...
function max_entropy(serie::AbstractArray{__FLOAT_TYPE,3}, ε_range::Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}, n::Int; sampling::Int=floor(Int, size(serie, 2) * 0.1))
    #   IN DEV
end