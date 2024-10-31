#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#

"""
    entropy(probs::Vector{__FLOAT_TYPE})

It calculates the Shanon's entropy using the probabilities of our microstates. 
We call it as Recurrence Entropy and it can be used to find the best value of ε
when we want to use the probabilities as input for some machine learning.
"""
function entropy(probs::Vector{__FLOAT_TYPE})
    s = 0.0
    for p in probs
        if (p > 0)
            s += (-1) * p * log(p)
        end
    end
    return s
end