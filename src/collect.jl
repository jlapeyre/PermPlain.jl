# surely, there is a more straightforward way
# These are meant to collect and convert arguments for use with PermutationsA

function collect1{T<:String}(a::Array{T,1})
    map(BigInt,a)
end

function collect1{T<:Real}(a::Array{T,1})
    return a
end

function collect1{T<:Real, V<:Real}(::Type{T}, a::Array{V,1})
    return T == V ? a : collect(T,a)  # note that it may return a copy!
end

function tupcollect{T}(::Type{T}, a::Tuple)
    aout = Array(Array{T,1},0)
    for x in a
        a1 = Array(T,0)
        for x1 in x
            push!(a1,convert(T,x1))
        end
        push!(aout,a1)
    end
    return aout
end

for (f,t) in ((:int32, Int32), (:int64, Int64), (:int, Int) , (:BigInt, BigInt))
    @eval begin
        function tupcollect(::Type{$t}, a::Tuple)
            aout = Array(Array{$t,1},0)
            for x in a
                a1 = Array($t,0)
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

function collect2{T}(a::Array{T,1}...)
    collect(map(collect1,a))
end

function collect2{T,V}(::Type{T}, a::Array{V,1}...)
    collect(map((x)->collect1(T,x) ,a))
end

