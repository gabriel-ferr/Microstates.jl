#
#       -- Microstates.jl: read more in https://github.com/gabriel-ferr/Microstates.jl
#


function microstates(serie::AbstractArray{__FLOAT_TYPE,2}, ε::Any, n::Int; sampling::Int=floor(Int, size(serie, 2) * 0.1), p_vector=power_vector(n), recurrence=standard_recurrence)

    #       Okay, I'm going to use a dict here to save some RAM, since most of
    #   the microstates for n > 3 don't exist. This way we can easily compute
    #   the probabilities for higher values of n, and we are less likely to
    #   run out of RAM. =3
    stats = Dict{Int64,__FLOAT_TYPE}()

    #       Well, some memory allocations =V
    p = 0
    add = 0
    counter = 0
    sz = size(serie, 2)

    row_counter = 0
    col_counter = 0

    #       Gets the random positions to set the microstates. =3
    index_row = sample(1:sz-(n-1), sampling)
    index_col = sample(1:sz-(n-1), sampling)

    #       Let's do it >.<
    #       Since index_row and index_col have the same number of elements, 
    #   we can use just one index to get them, so just one for here. =3
    for i in eachindex(index_row)
        add = 0

        row_counter = 0
        col_counter = 0

        #       Well, now we set one index to get the power vector position =D
        for n_index in eachindex(p_vector)
            add = add + p_vector[n_index] * recurrence(serie[:, index_row[i]+row_counter], serie[:, index_col[i]+col_counter], ε)

            col_counter += 1
            if (col_counter >= n)
                col_counter = 0
                row_counter += 1
            end
        end

        p = Int64(add) + 1
        stats[p] = get(stats, p, 0) + 1
        counter += 1
    end

    #       Ok, now I calculate the probalities...
    for k in keys(stats)
        stats[k] /= counter
    end

    return stats, counter
end