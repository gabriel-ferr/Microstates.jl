#
#       -- Microstates.jl: read more in https://github.com/gabriel-ferr/Microstates.jl
#

function microstates(serie::AbstractArray{__FLOAT_TYPE,2}, ε::Any, n::Int; sampling::Int=floor(Int, size(serie, 1) * 0.1), p_vector::Vector{Int64}=power_vector(n), recurrence::Function=standard_recurrence)
    if (n >= 8)
        println("How Microstates.jl uses Int64 you cannot use n >= 8 because Int64 doesn't support it.")
        return
    end

    #       I will to work with multithreading process here, because I think that
    #   the previous version of this code isn't doing the sampling correctly.
    #       So, first I define a task function...
    function __task_microstates(x_index)
        #
        #       a) Do a sampling about the columns...
        y_samples = sample(1:size(serie, 1)-n, sampling)
        #       b) Create a dict to store our result.
        result = Dict{Int64,__FLOAT_TYPE}()

        #       c) Alloc some memory =D
        p = 0
        add = 0
        counter = 0
        x_counter = 0
        y_counter = 0

        #       Ok, now...
        for x in x_index
            for y in y_samples
                add = 0
                x_counter = 0
                y_counter = 0

                for n in eachindex(p_vector)
                    add = add + p_vector[n] * recurrence(serie[x+x_counter, :], serie[y+y_counter, :], ε)

                    y_counter += 1
                    if (y_counter >= n)
                        x_counter += 1
                        y_counter = 0
                    end
                end

                p = Int64(add) + 1
                result[p] = get(result, p, 0) + 1
                counter += 1
            end
        end

        #       Return the dict =3
        return result, counter
    end

    #       Okay, I will to take the sample for x here...
    x_samples = sample(1:size(serie, 1)-n, sampling)

    #       Now, we partition it in the threads ...
    int_numb = trunc(Int, sampling / Threads.nthreads())
    par_numb = sampling - (int_numb * Threads.nthreads())

    #       Okay, now create the partition ranges...
    _numb_init = int_numb + (par_numb > 0 ? 1 : 0)
    if (par_numb > 0)
        par_numb -= 1
    end

    itr = [1:_numb_init]
    for i = 2:Threads.nthreads()
        _numb = int_numb + (par_numb > 0 ? 1 : 0)
        if (par_numb > 0)
            par_numb -= 1
        end
        push!(itr, (itr[i-1][end]+1):(itr[i-1][end]+_numb))
    end

    #       Free space...
    int_numb = Nothing
    par_numb = Nothing
    _numb_init = Nothing

    #       Okay, now we create the tasks...
    tasks = []
    for i = 1:Threads.nthreads()
        push!(tasks, Threads.@spawn __task_microstates(x_samples[itr[i]]))
    end

    #       Get the results...
    results = fetch.(tasks)
    #       Vector with our stats =3
    stats = zeros(Float64, 2^(n * n))
    cnt = 0
    #       To finish...
    for r in results
        cnt += r[2]
        stats[collect(keys(r[1]))] .= collect(values(r[1]))
    end
    #       Calc the probs...
    for i in eachindex(stats)
        stats[i] /= cnt
    end

    return stats, cnt
end

function microstates_lowsamples(serie::AbstractArray{__FLOAT_TYPE,2}, ε::Any, n::Int; sampling::Int=floor(Int, size(serie, 1) * 0.2), p_vector::Vector{Int64}=power_vector(n), recurrence::Function=standard_recurrence)

    if (n >= 8)
        println("How Microstates.jl uses Int64 you cannot use n >= 8 because Int64 doesn't support it.")
        return
    end

    #       Okay, I'm going to use a dict here to save some RAM, since most of
    #   the microstates for n > 3 don't exist. This way we can easily compute
    #   the probabilities for higher values of n, and we are less likely to
    #   run out of RAM. =3
    stats = Dict{Int64,__FLOAT_TYPE}()

    #       Well, some memory allocations =V
    p = 0
    add = 0
    counter = 0
    sz = size(serie, 1)

    row_counter = 0
    col_counter = 0

    #       Gets the random positions to set the microstates. =3
    index_row = sample(1:sz-(n-1), sampling)
    index_col = sample(1:sz-(n-1), sampling)

    #       Let's do it >.<
    #       Since index_row and index_col have the same number of elements, 
    #   we can use just one index to get them, so just one for here. =D
    for i in eachindex(index_row)
        add = 0

        row_counter = 0
        col_counter = 0

        #       Well, now we set one index to get the power vector position =D
        for n_index in eachindex(p_vector)
            add = add + p_vector[n_index] * recurrence(serie[index_row[i]+row_counter, :], serie[index_col[i]+col_counter, :], ε)

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

    #       Ok, now calculate the probalities...
    for k in keys(stats)
        stats[k] /= counter
    end

    return stats, counter
end