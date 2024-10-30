#
#           Microstates.jl
#               Gabriel Ferreira
#               Orientation: Sérgio Roberto Lopes, Thiago de Lima Prado
#
"""
    power_vector(n::Int)

Returns a vector that helps to convert the microstates from a binary format to a decimal format.
If you want to run the microstates in a loop, consider defining a constant power vector and passing it to the microstates function, such as:
```
    julia> const pvec = power_vector(2)
    julia> microstates(...; power_aux = pvec)
```
"""
function power_vector(n::Int)
    vec = zeros(Int, (n * n))
    for i in eachindex(vec)
        vec[i] = 2^(i - 1)
    end
    return vec
end

"""
    microstates(serie::AbstractArray{__FLOAT_TYPE,2}, ε::Any, n::Int; samples=floor(Int, size(serie, 2)), power_aux::Vector{Int}=power_vector, recurrence::Function=std_recurrence)

Calculates the probabilities of microstates from a series and returns these probabilities and the number of samples.
Instead of computing the recurrence matrix and taking the microstates after it, we only compute the recurrence in each microstate block.

If you want to run this function in a loop, consider defining a constant power vector and passing it to the function, i.e:
```
    julia> const pvec = power_vector(2)
    julia> microstates(...; power_aux = pvec)
```
This will to free up some time with Julia's garbage collector.

The serie format is like the fomat that we usually use in the physics, i.e:
```
    julia> serie = rand(Float64, 2, 1000)
    2x1000 Matrix{Float64}:
    0.562823    0.796586    0.369896    ...     0.745746
    0.381401    0.544322    0.003962    ...     0.188805
```
Or, if you are using DifferentialEquations.jl and take the result from `data = solve(...)` as `data[:, :]` the format is similar.
Where each row is a dimension and each column is a timestamp.

The threshold (`ε`) value have a Any format because it is defined by recurrence function. If you use the standard recurrence `recurrence::Function=std_recurrence`
The threshold value (`ε`) has an Any format because it is defined by the recurrence function. If you use the default recurrence `recurrence::Function=std_recurrence' 
the ε value must be a Float64, but if you need to use other recurrence functions you can change the recurrence function and pass an appropriate ε value, i.e:
```
    julia> function corridor(x::Vector{Float64}, y::Vector{Float64}, ε::Tuple{Float64, Float64})
                return ((ε[2] - euclidean(x, y)) * (euclidean(x, y) - ε[1])) >= 0 ? Int8(1) : Int8(0)
           end
    julia> microstate(serie, (0.05, 0.2), n; recurrence=corridor)
```
The `n` parameter specifies the size of the microstate. As we are using Int64 the value of `n` must be between 2 and 7.

Finally, the `samples` defines the number of "rows"  that we take, and for each row we take the same number of "columns", both randomly.
"""
function microstates(serie::AbstractArray{__FLOAT_TYPE,2}, ε::Any, n::Int; samples=floor(Int, size(serie, 2)), power_aux::Vector{Int}=power_vector, recurrence::Function=std_recurrence)
    if (n < 2 || n > 7)
        println("As Microstates.jl uses Int64, you cannot use n > 7 because Int64 does not support it. And n < 2 makes no sense!")
        return
    end
    #
    #       As I need to use an async process to make it faster, I make here a "task function".
    function __task_microstates(x_index)
        #
        #       a) Does a sampling over the columns...
        y_samples = sample(1:size(serie, 2)-n, samples)
        #       b) Creates a dict to store our result...
        result = Dict{Int64,__FLOAT_TYPE}()
        #       c) Alloc some memory =D
        p = 0
        add = 0
        counter = 0
        x_counter = 0
        y_counter = 0
        #
        #       Ok...
        for x in x_index
            for y in y_samples
                add = 0
                x_counter = 0
                y_counter = 0

                for n in eachindex(power_aux)
                    add = add + power_aux[n] * recurrence(serie[:, x+x_counter], serie[:, y+y_counter], ε)

                    y_counter += 1
                    if (y_counter >= n)
                        x_counter += 1
                        y_counter = 0
                    end
                end

                p = add + 1
                result[p] = get(result, p, 0) + 1
                counter += 1
            end
        end
        #
        #       Returns the dict...
        return result, counter
    end

    #       Okay, first I take the samples for x...
    x_samples = sample(1:size(serie, 2)-n, samples)

    #       Now, we need to partition based on the number of threads xD
    int_numb = trunc(Int, samples / Threads.nthreads())
    par_numb = samples - (int_numb * Threads.nthreads())

    #       Creates the partition ranges and puts them into a vector...
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

    #       Free memory...
    int_numb = Nothing
    par_numb = Nothing
    _numb_init = Nothing

    #       Okay, now we create the tasks...
    tasks = []
    for i = 1:Threads.nthreads()
        push!(tasks, Threads.@spawn __task_microstates(x_samples[itr[i]]))
    end

    #       Wait and get the result of the tasks =3
    results = fetch.(tasks)
    #       Alloc memory for all microstate that we can have.
    stats = zeros(Float64, 2^(n * n))
    cnt = 0
    #       Groups the results of the tasks...
    for r in results
        cnt += r[2]
        stats[collect(keys(r[1]))] .= collect(values(r[1]))
    end
    #       To finish, calculates the probabilities =D
    for i in eachindex(stats)
        stats[i] /= cnt
    end
    #
    return stats, cnt
end