#
#       Microstates.jl - utils
#
#           Okay..., here, we have a function that calculates the entropy of our microstates of recurrence.
function entropy(probs::Vector{Float64})
    s = 0.0
    for p in probs
        if (p >= 0)
            s += (-1) * p * log(p)
        end
    end
    return s
end