#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#
##      -- STD Recurrence
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


    findthreshold(serie::AbstractArray{__FLOAT_TYPE,3}, n::Int, thres::Tuple{Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int},Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}}; multiplier::Tuple{Int,Int}=(2, 2), pvec=power_vector(n))

This variation of the same functions do the same thing but using the corridor recurrence. If you want to use
corridor recurrence and this function to find the `ε` value, the ideia is equals but you need to pass the threshold
range as a parameter `thres::Tuple{Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int},Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}}`, i.e:

```
    julia> findthreshold(serie, n, ((0.0, 0.1, 10), (0.0, 1.0, 20)))
```

where the first tuple says the range of `ε_min` and the second of `ε_max`.
"""
function findthreshold(serie::AbstractArray{__FLOAT_TYPE,3}, n::Int; thres::Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}=(0.0, 1.0, 20), multiplier::Int=4, pvec=power_vector(n))
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

##      -- CRD Recurrence
function findthreshold(serie::AbstractArray{__FLOAT_TYPE,3}, n::Int, thres::Tuple{Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int},Tuple{__FLOAT_TYPE,__FLOAT_TYPE,Int}}; multiplier::Tuple{Int,Int}=(2, 2), pvec=power_vector(n))
    #
    #       Now we have 2 ranges, then...
    ε_min_big = range(thres[1][1], thres[1][2], thres[1][3])
    ε_max_big = range(thres[2][1], thres[2][2], thres[2][3])
    #       Our status now needs to have 2 dimensions, and only
    #   will to save the entropy values...
    stats = zeros(__FLOAT_TYPE, length(ε_min_big), length(ε_max_big))
    #       And one vector to store the entropy of each sample...
    for ε_max in eachindex(ε_max_big)
        Threads.@threads for ε_min in eachindex(ε_min_big)
            if (ε_min_big[ε_min] >= ε_max_big[ε_max])
                continue
            end

            etr = zeros(__FLOAT_TYPE, size(serie, 3))
            for i in eachindex(etr)
                #       Get the microstates probabilities...
                probs, _ = microstates(serie[:, :, i], (ε_min_big[ε_min], ε_max_big[ε_max]), n; power_aux=pvec, recurrence=crd_recurrence)
                etr[i] = entropy(probs)
            end

            #       The logic used to standard recurrence don't work here, so we need to
            #   compute all entropies and get the max of the matrix =<
            stats[ε_min, ε_max] = mean(etr)
        end
    end

    #       So.. get the max =D
    max = findmax(stats)[2]
    around = [(ε_min_big[max[1] > 1 ? max[1] - 1 : max[1]], ε_min_big[max[1] < length(ε_min_big) ? max[1] + 1 : max[1]]), (ε_max_big[max[2] > 1 ? max[2] - 1 : max[2]], ε_max_big[max[2] < length(ε_min_big) ? max[2] + 1 : max[2]])]

    #       Now, to accuracy...
    ε_min_small = range(around[1][1], around[1][2], thres[1][3] * multiplier[1])
    ε_max_small = range(around[2][1], around[2][2], thres[2][3] * multiplier[2])
    stats = zeros(__FLOAT_TYPE, length(ε_min_small), length(ε_max_small))

    for ε_max in eachindex(ε_max_small)
        Threads.@threads for ε_min in eachindex(ε_min_small)
            if (ε_min_small[ε_min] >= ε_max_small[ε_max])
                continue
            end

            etr = zeros(__FLOAT_TYPE, size(serie, 3))
            for i in eachindex(etr)
                #       Get the microstates probabilities...
                probs, _ = microstates(serie[:, :, i], (ε_min_small[ε_min], ε_max_small[ε_max]), n; power_aux=pvec, recurrence=crd_recurrence)
                etr[i] = entropy(probs)
            end

            #       The logic used to standard recurrence don't work here, so we need to
            #   compute all entropies and get the max of the matrix =<
            stats[ε_min, ε_max] = mean(etr)
        end
    end

    #       Get the max again...
    max = findmax(stats)
    return (max[1], (ε_min_small[max[2][1]], ε_max_small[max[2][2]]))
end