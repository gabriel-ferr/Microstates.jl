#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#
"""
    std_recurrence(x::AbstractVector{__FLOAT_TYPE}, y::AbstractVector{__FLOAT_TYPE}, ε::__FLOAT_TYPE)

Use the standard recurrence to calculate the recurrence between `x` and `y`, i.e:
    `R = Θ(ε - |x-y|)`
"""
function std_recurrence(x::AbstractVector{__FLOAT_TYPE}, y::AbstractVector{__FLOAT_TYPE}, ε::__FLOAT_TYPE)
    return (ε - euclidean(x, y)) >= 0 ? Int8(1) : Int8(0)
end