#
#       Microstates.jl
#       In this script, we have several functions to compute the microstates from various types of data.
#
# -------------------------------------------------------------------------------------------
#           The power vector helps convert microstates from a binary form to a decimal, which
#   we use as an index for a vector.
function power_vector(n::Int)
    return power_vector((n, n))
end
# -------------------------------------------------------------------------------------------
#           Here, I create a general function to compute a power vector of any size,
#   while limiting the microstate "area" to 64.
function power_vector(sz::Tuple{Vararg{Int}})
    #
    #       First, I will calculate the "area" of our microstate. The `sz` tuple receives 
    #   the size of each side of the microstate, so I calculate the area by multiplying the
    #   sizes of each side, because it is a rectangle =D
    microstate_area = 1
    for a in sz
        microstate_area *= a
    end
    #
    #       Now, since julia uses Int64, we cannot have an area greater than 64, so I will check it.
    if (microstate_area >= 64)
        throw("Julia uses Int64, so the microstate space cannot be greater than 64.")
    end
    #
    #       Finally, we create the power vector =D
    vect = zeros(Int, microstate_area)
    for i in eachindex(vect)
        vect[i] = 2^(i - 1)
    end
    #
    return vect
end
# -------------------------------------------------------------------------------------------
#           This is the default function used to compute the microstates of a time series.
function microstates(data::AbstractArray{Float64, 2}, threshold::Any, n::Int;
    samples_percent::Float64 = 0.2, vect::Vector{Int64} = power_vector(n), recurr::Function = std_recurrence)
    return microstates(data, data, threshold, (n, n); samples_percent = 0.2, vect = power_vector((n, n)), recurr  = std_recurrence)
end
# -------------------------------------------------------------------------------------------
#           And here, I create a general version to compute the microstates of "simple data".
#   More specifically, this function calculates the microstates of a CRP.
function microstates(data_x::AbstractArray{Float64, 2}, data_y::AbstractArray{Float64, 2}, threshold::Any, sz::Tuple{Int, Int}; 
    samples_percent::Float64 = 0.2, vect::Vector{Int64} = power_vector(sz), recurr::Function = std_recurrence)
    #
    side_samples = samples_percent ^ (1/length(sz))
    #
    x_samples = sample(1:size(data_x, 2)-(sz[1] - 1), Int(ceil(size(data_x, 2) * side_samples)))
    #
    stats = zeros(Float64, 2^(sz[1] * sz[2]), Threads.nthreads())
    #
    int_numb = trunc(Int, Int(ceil(size(data_x, 2) * side_samples)) / Threads.nthreads())
    par_numb = Int(ceil(size(data_x, 2) * side_samples)) - (int_numb * Threads.nthreads())

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

    int_numb = Nothing
    par_numb = Nothing
    _numb_init = Nothing

    function __async_compute(x_index, th_index)
        #
        y_samples = sample(1:size(data_y, 2)-(sz[2]-1), Int(ceil(size(data_y, 2) * side_samples)))
        #
        add = 0
        x_counter = 0
        y_counter = 0
        counter = 0
        #
        for x in x_index
            y_samples = sample(1:size(data_y, 2)-(sz[2]-1), Int(ceil(size(data_y, 2) * side_samples)))
            for y in y_samples
                add = 0
                x_counter = 0
                y_counter = 0

                for m in eachindex(vect)
                    add = add + vect[m] * recurr(data_x[:, x+x_counter], data_y[:, y+y_counter], threshold)
                    y_counter += 1
                    if (y_counter >= sz[2])
                        x_counter += 1
                        y_counter = 0
                    end
                end

                p = add + 1
                stats[p, th_index] += 1
                counter += 1
            end
        end

        return counter
    end

    tasks = []
    for i = 1:Threads.nthreads()
        push!(tasks, Threads.@spawn __async_compute(x_samples[itr[i]], i))
    end

    results = fetch.(tasks)
    cnt = 0

    res = zeros(Float64, 2^(sz[1] * sz[2]))
    for i = 1:Threads.nthreads()
        res .= res .+ stats[:, i]
        cnt += results[i]
    end

    stats ./= cnt

    return stats, cnt
end
# -------------------------------------------------------------------------------------------