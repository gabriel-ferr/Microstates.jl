#
#       -- Microstates.jl: read more in https://github.com/gabriel-ferr/Microstates.jl
#

#           Okay... I'll try make the same logic that I used for CPU and adapt it to GPU...
#   But, first: I cannot use if in GPU, so I need compute the recorrences in CPU, organize
#   them and pass it to GPU convert to decimal. Okay... I'll try >.<
function gpu_microstates(serie::AbstractArray{__FLOAT_TYPE,2}, ε::Any, n::Int; sampling::Int=floor(Int, size(serie, 2) * 0.1))
    #       How I'll use Int32 in GPU, we cannot use n > 5 =3
    if (n > 5)
        println("How Microstates.jl GPU uses Int32 you cannot use n >= 6 because Int32 doesn't support it.")
        return
    end

    if (!CUDA.functional())
        println("CUDA isn't avaliable, please use CPU function: `microstates(..)`.")
        return
    end

    sz = size(serie, 2)
    index_row = sample(1:sz-(n-1), sampling)
    index_col = sample(1:sz-(n-1), sampling)

    states = []
    Threads.@threads for i in eachindex(index_row)
    end
end