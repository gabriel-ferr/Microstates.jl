#
#       -- Microstates.jl: read more in https://github.com/gabriel-ferr/Microstates.jl
#

#       Shenon entropy for microstates =D
function entropy(probs::Dict{Int64,__FLOAT_TYPE})
    s = 0.0
    for k in keys(probs)
        if (probs[k] > 0)
            s += (-1) * probs[k] * log(probs[k])
        end
    end
    return s
end