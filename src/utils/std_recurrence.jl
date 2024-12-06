#
#       Microstates.jl - utils
#
#           To compute the microstates, I need to define the recurrence function. 
#   Here, I created a function to apply the default form of recurrence.
function std_recurrence(x::AbstractVector{Float64}, y::AbstractVector{Float64}, threshold::Float64)
    return (threshold - euclidean(x, y)) >= 0 ? 1 : 0
end