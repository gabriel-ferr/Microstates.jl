#
#       Microstates.jl - utils
#
#           To compute the microstates, I need to define the recurrence function. Here, I created a function
#   to apply corridor recurrence (which is the same as recurrence using a corridor threshold).
#
#           R = Θ(ε_max - ||x-y||) * Θ(||x-y|| - ε_min)
function crd_recurrence(x::AbstractVector{Float64}, y::AbstractVector{Float64}, threshold::Tuple{Float64, Float64})
    return ((threshold[2] - euclidean(x, y)) >= 0 ? 1 : 0) * ((euclidean(x,y) - threshold[1]) >= 0 ? 1 : 0)
end