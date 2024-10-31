#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#
"""
    findthreshold(serie::AbstractArray{__FLOAT_TYPE,3}, n::Int; thres::Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}=(0.0, 1.0, 20), multiplier::Int=4)

Find a good value for `ε` using some samples. So, the parameter `serie` here is the same that for `microstates(...)`
for the dimensions 1 and 2, but here we have a third dimension that is the samples that we are using for
find the maximum entropy of recurrence.

This function first finds the points next to a maximum entropy without a good accuracy. For it,
we define a Tuple `thres::Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}` with 3 elements. The elements at
index 1 and 2 are the `ε` range that we are taking, and the third element is the number of elements
that we are taking to find "the points next to the maximum".

When we find those points, we use it as parameter to take an other range between them, with
a number of elements equals to the previous number multiply by the parameter `multiplier`.
Doing it we try to get a value of ε that has more accuracy =D

If the threshold's range that you pass to the function does not have the maximum entropy
into it, the function will to return a message saying it.
"""
function findthreshold(serie::AbstractArray{__FLOAT_TYPE,3}, n::Int; thres::Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}=(0.0, 1.0, 20), multiplier::Int=4)
    #
    #       How I am going to work with some loops, I define a constant
    #   power vector here.
    pvec = power_vector(n)
    #
    #       Get a range of values for our threshold...
    ε_big = range(thres[1], thres[2], thres[3])
    #   ... and I define too a vector here that will store the entropy
    #   that we calculate. We have here 3 elements. The element [1] is
    #   our point before the max entropy finded, [2] is the max entropy
    #   and [3], the point after the max entropy.
    stats = [(0.0, 0.0), (0.0, 0.0), 0.0]
    #       I will to alloc some temp memory here to store the entropy values
    #   for each sample...
    etr = zeros(__FLOAT_TYPE, size(serie, 3))
    #
    #   The loop...
    for ε in ε_big
        #       Reset our temp memory...
        etr = zeros(__FLOAT_TYPE, size(serie, 3))

        for i in eachindex(etr)
            #       Get the microstates probabilities...
            probs, _ = microstates(serie[:, :, i], ε, n; power_aux=pvec)
            etr[i] = entropy(probs)
        end

        s = mean(etr)
        #       Here I check if the mean entropy of our samples is smaller than
        #   the max entropy finded before; if it is small, we know that the
        #   previous value of entropy is the maximum local, so we can register it
        #   an break the loop.
        if (s < stats[2][2])
            stats[3] = (ε, s)
            break
        end
        #       If the max does not have found, we change the max value that we have in the
        #   positions [2] and move the previous value to index [1].
        stats[1] = stats[2]
        stats[2] = (ε, s)
    end

    #       If the loop ends and we do not have a max, the position [3] of stats will to be NaN...
    #   ... so, if we have it, print a error message and exit =3
    if (typeof(stats[3]) == Float64)
        println(string("The max entropy cannot to be find in the threhold range 'thres=", thres, "'."))
        return
    end

    #       Now, we do the same thing again to get more accuracy xD
    ε_small = range(stats[1][1], stats[3][1], thres[3] * multiplier)
    stats = [(0.0, 0.0), (0.0, 0.0), NaN]
    for ε in ε_small
        etr = zeros(__FLOAT_TYPE, size(serie, 3))

        for i in eachindex(etr)
            #       Get the microstates probabilities...
            probs, _ = microstates(serie[:, :, i], ε, n; power_aux=pvec)
            etr[i] = entropy(probs)
        end

        s = mean(etr)
        if (s < stats[2][2])
            stats[3] = (ε, s)
            break
        end

        stats[1] = stats[2]
        stats[2] = (ε, s)
    end

    return stats[2]
end