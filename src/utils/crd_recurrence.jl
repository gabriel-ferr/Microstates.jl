#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#
"""
    crd_recurrence(x::AbstractVector{__FLOAT_TYPE}, y::AbstractVector{__FLOAT_TYPE}, ε::Tuple{__FLOAT_TYPE, __FLOAT_TYPE})

Use the corridor recurrence to calculate the recurrence between `x` and `y`, i.e:
    `R = Θ(ε_max - |x-y|) * Θ(|x-y| - ε_min)`
"""
function crd_recurrence(x::AbstractVector{__FLOAT_TYPE}, y::AbstractVector{__FLOAT_TYPE}, ε::Tuple{__FLOAT_TYPE,__FLOAT_TYPE})
    return ((ε[2] - euclidean(x, y)) * (euclidean(x, y) - ε[1])) >= 0 ? Int8(1) : Int8(0)
end