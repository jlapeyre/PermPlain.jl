# surely, there is a more straightforward way
# These are meant to collect and convert arguments for use with PermutationsA

collect1(a::Vector{<:String}) =  map(BigInt,a)
collect1(a::Vector{<:Real}) = a

function collect1(::Type{T}, a::Vector{V}) where {T, V}
    return T == V ? a : collect(T, a)  # note that it may return a copy!
end

function tupcollect(::Type{T}, a::Tuple) where T
    aout = Array{Array{T,1}}(0)
    for x in a
        a1 = Array{T}(0)
        for x1 in x
            push!(a1,convert(T,x1))
        end
        push!(aout,a1)
    end
    return aout
end

for (f,t) in ((:int32, Int32), (:int64, Int64), (:BigInt, BigInt))
    @eval begin
        function tupcollect(::Type{$t}, a::Tuple)
            aout = Array{Array{$t,1}}(0)
            for x in a
                a1 = Array{$t}(0)
                for x1 in x
                    push!(a1, ($f)(x1))
                end
                push!(aout,a1)
            end
            return aout
        end
    end
end

function tupcollect(a::Tuple)
    tupcollect(Int,a)
end

function collect2(a::AbstractVector...)
    collect(map(collect1,a))
end

function collect2(::Type{T}, a::AbstractVector...) where T
    collect(map((x)->collect1(T,x) ,a))
end
