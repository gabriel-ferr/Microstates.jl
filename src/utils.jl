#
#       -- Microstates.jl: read more in https://github.com/gabriel-ferr/Microstates.jl
#

#           Power vector used to convert the microstates from a binary format to 
#   decimal format, then we can use this decimal value as a key to identify the microstate =D
function power_vector(n::Int)
    vec = zeros(Int64, (n * n))
    for i in eachindex(vec)
        vec[i] = Int64(2^(i - 1))
    end
    return vec
end

#           Function to Standard recurrence...
function standard_recurrence(x::AbstractVector{__FLOAT_TYPE}, y::AbstractVector{__FLOAT_TYPE}, ε::__FLOAT_TYPE)
    return (ε - euclidean(x, y)) >= 0 ? 1 : 0
end