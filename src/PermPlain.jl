module PermPlain

using DataStructures: counter
import LinearAlgebra
using SparseArrays: sparse
using Random
import Base: isperm

include("collect.jl")

export permlist, permcycles, permsparse, permmatrix # whether to export, and what ?

randperm(::Type{T}, n::Integer) where T = collect(T,randperm(n))

AbstractVectorVector{T} = AbstractVector{V} where {V <: (AbstractVector{T} where {T})}

## permcycles. Find cyclic decomposition  ##

# compute cyclic decomposition.
# builds a cycle list in the canonical order.
# See pari for possible efficiency improvement.
function permcycles(p::AbstractVector{T}) where T
    n = length(p)
    visited = falses(n)
    cycles = Array{Vector{T}}(undef, 0)
@inbounds for k in convert(T,1):convert(T,n)
        if ! visited[k]
            knext = k
            cycle = Array{T}(undef, 0)
            while ! visited[knext]
                push!(cycle,knext)
                visited[knext] = true
                knext = p[knext]
            end
            length(cycle) > 1 && push!(cycles,cycle)
        end
    end
    return cycles
end

permcycles(m::AbstractMatrix) = permcycles(permlist(m))
permcycles(sp::Dict{T,T}) where T = sparsetocycles(sp)

function sparsetocycles(sp::Dict{T,T}) where T
    cycs = Array{Vector{T}}(undef, 0)
    length(sp) == 0 && return cycs
    ks = collect(keys(sp))
    n = length(ks)
    seen = Dict{T,Bool}()
    for k in ks seen[k] = false end
@inbounds   k = ks[1]
    nseen = 0
@inbounds while nseen <= n
        foundunseen = false
        for i in ks  # could be more efficient starting search after last found
            if seen[i] == false
                k = i
                foundunseen = true
                break
            end
        end
        foundunseen == false && break
        kcyclestart = k
        cyc = Array{T}(undef, 0)
        while true
            push!(cyc,k)
            if seen[k] == true
                error("Algorithm error: double setting seen k=$k, nseen=$nseen, kcyclestart=$kcyclestart")
            end
            seen[k] = true
            k = sp[k]
            nseen = nseen + 1
            if k == kcyclestart
                break
            end
        end
        push!(cycs,cyc)
    end
    return cycs
end

## permlist. permutation in single-list form ##

permlist(m::AbstractMatrix) = mattoperm(m)
mattoperm(m::AbstractMatrix{T}) where T = mattoperm!(Vector{T}(undef, size(m, 1)), m)

# FIXME: p is sometimes an Array{:Bool}.
function mattoperm!(p, m::AbstractMatrix{T}) where T
    n = size(m)[1]
    maxk = zero(T)
@inbounds  for i in 1:n
        for j in 1:n
            if m[j,i] != 1
                continue
            end
            p[j] = i
        end
    end
    p
end

permlist(cycs::AbstractArray{Array{T, 1}, 1}, pmax::Real = 0) where T <: Real =  cycstoperm(cycs,pmax)

function cycstoperm(cycs::AbstractArray{Array{T,1},1}, pmax::Integer = 0) where T
    isempty(cycs) && return [one(T):convert(T,pmax)]
    cmaxes = [maximum(c) for c in cycs]
    cmax = maximum(cmaxes)  # must be a faster way
#    perm = [one(T): (pmax > cmax ? convert(T,pmax) : convert(T,cmax))]
    perm = collect(one(T): (pmax > cmax ? convert(T,pmax) : convert(T,cmax)))
@inbounds  for c in cycs
        for i in 2:length(c)
            perm[c[i-1]] = c[i]
        end
        perm[c[end]] = c[1]
    end
    perm[cmax+1:pmax] = collect(cmax+1:pmax)
    return perm
end

permlist(sp::Dict{T, T}) where T = sparsetolist(sp::Dict{T, T})

function sparsetolist(sp::Dict{T, T}) where T
    p = collect(one(T):convert(T,maximum(sp)[1]))
@inbounds for (i,v) in sp
        p[i] = v
    end
    return p
end

## permsparse ##

permsparse(m::AbstractMatrix) = mattosparse(m)

function mattosparse(m::AbstractMatrix{T}) where T
    p = Dict{T, T}()
    return mattoperm!(p, m), maximum(p)[1]
end

permsparse(p::AbstractVector) = listtosparse(p)

function listtosparse(p::AbstractVector{T}) where T
    data = Dict{T,T}()
    maxk = zero(T)
    length(p) == 0 && return (data,maxk)
@inbounds for i in p
        pv = p[i]
        if pv != i
            data[i] = pv
            pv > maxk ? maxk = pv : nothing
        end
    end
    return (data,maxk)
end

permsparse(cycs::AbstractArray{Array{T, 1}, 1}) where T = cycstosparse(cycs)

function cycstosparse(cycs::AbstractArray{Array{T, 1}, 1}) where T
    data = Dict{T,T}()
    maxk = zero(T)
@inbounds for c in cycs
        pv = c[1]
        data[c[end]] = pv
        c[end] > maxk ? maxk = c[end] : nothing
        for i in 1:length(c)-1
            pv = c[i]
            data[pv] = c[i+1]
            pv > maxk ? maxk = pv : nothing
        end
    end
    return (data,maxk)
end

## permmatrix ##

permmatrix(p::AbstractVector{<:Real}, sparse::Bool = false) = permtomat(p, sparse)
# Convert PLIST to PMAT
function permtomat(p::AbstractVector{T}, sparse::Bool = false) where T <: Real
    n::T = length(p)
    A = sparse ? sparse(LinearAlgebra.I*one(T), n, n) : Matrix{T}(LinearAlgebra.I, n, n)
    return A[p,:]
end

function permmatrix(sp::Dict{T, T}) where T
    n = convert(T,maximum(sp)[1])
    ot = one(T)
    z = zero(T)
    m = Matrix(LinearAlgebra.I, n, n)
@inbounds for (i, v) in sp
        m[i,i] = z
        m[i,v] = ot
    end
    return m
end

# Convert cyclic decomposition to canonical form
# used by gap, Mma, and Arndt thesis.
# Note that Arndt uses a different order internally to store the cycles as a single list.
function canoncycles(cycs::AbstractArray{Array{T,1},1}) where T
    ocycs = Array{Array{T,1}}(undef, 0)
    for cyc in cycs
        push!(ocycs,circshift(cyc, -argmin(cyc) + 1))
    end
    sort!(ocycs,by=(x)->x[1])
    return ocycs
end

## cyclelengths. Find cyclic decomposition, but only save cycle lengths ##

function cyclelengths(p::AbstractVector{T}) where T
    n = length(p)
    visited = falses(n)
    lengths = Array{Int}(undef, 0)
@inbounds for k in one(T):convert(T,n)
        if ! visited[k]
            knext = k
            len = 0
            while ! visited[knext]
                len += 1
                visited[knext] = true
                knext = p[knext]
            end
            len > 1 && push!(lengths,len)
        end
    end
    return lengths
end

# FIXME: Uh. fix this type specification
#cyclelengths(c::AbstractArray{Array{T,1},1}) where T = [length(x) for x in c]

# better
# cyclelengths(c::AbstractVector{V}) where {V <: AbstractVector} = [length(x) for x in c]

# best
cyclelengths(c::AbstractVector{<:AbstractVector}) = [length(x) for x in c]

# Gives cyclelengths in canonical order.
function cyclelengths(sp::Dict{T, T}) where T
    cyclens = Array{Int}(undef, 0)
    isempty(sp) && return cyclens
    ks = sort(collect(keys(sp)))
    n = length(ks)
    seen = Dict{T, Bool}()
    for k in ks seen[k] = false end
    k = ks[1]
    nseen = 0
    while nseen <= n
        didsee = false
        for i in ks
            if seen[i] == false
                k = i
                didsee = true
                break
            end
        end
        didsee == false && break
        k1 = k
        nincyc = 0
        while true
            nincyc += 1
            if seen[k] == true
                error("Algorithm error: double setting seen k=$k, nseen=$nseen, k1=$k1")
            end
            seen[k] = true
            k = sp[k]
            nseen = nseen + 1
            if k == k1
                break
            end
        end
        push!(cyclens,nincyc)  # inefficient
    end
    return cyclens
end

## numcycles ##

# These are not used by anything now.
numcycles(p) = length(cyclelengths(p))
# Wikipedia says some authors also require no fixed points.
# Here, we allow fixed points
iscyclic(p) = numcycles(p) == 1

## cycle type, sign ##

cycletype(p) = counter(cyclelengths(p))

permsgn_from_lengths(lens) = (-1)^(length(lens)+sum(lens))
permsgn(p) = permsgn_from_lengths(cyclelengths(p))

function permorder_from_lengths(clengths)
    result = 1
    for c in clengths
        result = lcm(result, c)
    end
    return result
end

permorder(p) = permorder_from_lengths(cyclelengths(p))

# distance between two PLISTs
# could use a macro for this and ==.
# is there a penalty for using swap macro on p and q instead of two branches ?
function permdistance(p::AbstractVector{T}, q::AbstractVector{T}) where T
    lp = length(p)
    lq = length(q)
    count = 0
    if lp < lq
        for i in 1:lp
            p[i] != q[i] ? count += 1 : nothing
        end
        for i in lp+1:lq
            q[i] != i ? count += 1 : nothing
        end
    else  # could factor code with refs, prbly not worth the trouble
        for i in 1:lq
            p[i] != q[i] ? count += 1 : nothing
        end
        for i in lq+1:lp
            p[i] != i ? count += 1 : nothing
        end
    end
    return count
end

# composition (multiplication) of two PLISTs
function permcompose(q::AbstractVector{T}, p::AbstractVector{<:Real}) where T
    lp = length(p)
    lq = length(q)
    lp == lq && return q[p]  # prbly not much time saved
    lr = lp < lq ? lq : lp
    r = Array{T}(lr)
    if lp <= lq
        r[1:lp] = q[p]
        r[lp+1:lq] = q[lp+1:lq]
    else
        for k in 1:lp
            i = p[k]
            r[k] = (i < lq ? q[i] : i) # parens for clarity
        end
    end
    return r
end

# inplace composition (multiplication) of two PLISTs
# first argument q is updated.
function permcompose!(q::AbstractVector, p::AbstractVector)
    lp = length(p)
    lq = length(q)
    if lp > lq
        error("Can't do inplace composition of permutation with longer permutation.")
    end
    for k in 1:lp
        q[k] = q[p[k]]
    end
    return q
end

function permcompose(q::Dict{T, T}, p::Dict{V, V}) where {T, V}
    dout = Dict{T,T}()
    z = zero(T)
    maxk = z
    seen = Dict{T,Bool}()
    for k in keys(q) seen[k] = false end  # wasteful!
    for (k,pofk) in p
        qofpofk = get(q,pofk,z)
        seen[pofk] = true
        k == qofpofk ? continue : nothing  # ignore 1-cycles
        dout[k] = (qofpofk == z ? pofk : qofpofk)
          qofpofk > maxk ? maxk = qofpofk : nothing
    end
    for j in keys(q)
        seen[j] ? nothing : dout[j] = q[j]
    end
    return dout, maxk
end

## permapply ##

# Vector ?
function permapply(q::Dict{T, T}, a::AbstractArray) where T
    aout = copy(a)
    len = length(aout)
    for (k,v) in q
        if k <= len && v <= len
            aout[k] = a[v]
        end
    end
    return aout
end

function permapply(q::AbstractVector, a::AbstractArray)
    aout = copy(a)
    lenq = length(q)
    lena = length(a)
    len =  lenq < lena ? lenq : lena
    for k in 1:len
        v = q[k]
        v <= len  && (aout[k] = a[v])
    end
    return aout
end

permapply(q::Union{Dict{T, T}, AbstractVector{T}}, a::String) where T = String(permapply(q,a.data))

# power of PLIST. output is PLIST
# The method permpower below is usually preferred because it does less allocation
function permpower2(p::AbstractVector{T}, n::Integer) where T
    n == 0 && return [one(T):convert(T,length(p))]
    n == 1 && return copy(p) # for consistency, don't return ref
    n < 0  && return permpower2(invperm(p),-n)
    q = permpower2(p, Int(floor(n/2)))
    q = q[q]
    return iseven(n) ? q : p[q]
end

function permpower!(p::AbstractVector{T},
                    pret::AbstractVector{T},
                    ptmp::AbstractVector{T},
                    n::Integer) where T
    onep = one(T)
    lenp = convert(T,length(p))
    n == 0 && (for i in onep:lenp pret[i] = i end; return )
    n < 0  && (permpower!(invperm(p), pret, ptmp, -n); return )
    n == 1 && (copy!(pret,p); return)
    permpower!(p, ptmp, pret,Int(floor(n/2)))
    if iseven(n)
        for i in onep:lenp pret[i] = ptmp[ptmp[i]] end
    else
        for i in onep:lenp pret[i] = p[ptmp[ptmp[i]]] end
    end
end

# This does less allocation (in general) than permpower2, and is faster in benchmarks
function permpower(p::AbstractVector{T}, n::Integer) where T
    n == 0 && return [one(T):convert(T,length(p))]
    n == 1 && return copy(p) # for consistency, don't return ref
    pret = similar(p)
    ptmp = similar(p)
    permpower!(p,pret,ptmp,n)
    return pret
end

# Compute power of permutation. Both input and output are PCYC
# Translated from pari perm.c
# Careful of degeneracy, and empty array may be returned.
function permpower(cyc::AbstractVector{<:AbstractVector{T}}, exp::Integer) where {T}
    r = 1
    for j in 1:length(cyc)
        r += gcd(length(cyc[j])-1,exp)
    end
    c = Vector{Vector{T}}(undef, 0)
    for j in 1:length(cyc)
        v = cyc[j]
        n = length(v)
        e = mod(exp,n)
        g = gcd(n,e)
        m = div(n,g)
        if m == 1 continue end
        for i in 1:g
            p = Array{Int}(undef, 0)
            l = i
            for k in 1:m
                push!(p,v[l])
                l += e
                if l > n  l -= n  end
            end
            push!(c,p)
        end
    end
    return c
end

permpower(q::Dict{T, T}, exp::Integer) where T = cycstosparse(permpower(sparsetocycles(q),exp))

# power of PCYC. output is PLIST
# see pari perm.c
function cyc_pow_perm(cyc::AbstractVector{<:AbstractVector{T}}, exp::Integer) where {T}
    n = 0
    cmaxes = [maximum(c) for c in cyc]
    n = maximum(cmaxes)
# This routine assumes cycles of length 1 are included. Still don't understand it.
#    for j = 1:length(cyc)
#        n += length(cyc[j])+1
#    end
    p = collect(T,1:n) # wasteful
    for j = 1:length(cyc)
        v = cyc[j]
        n = length(v)
        e = mod(exp,n)
        l = e
        for k in 1:n
            nv = v[l+1]
            ind = v[k]
            if (ind>length(p))
                error("ind is $ind length p is $(length(p))")
            end
            p[ind] = nv
            l += 1
            l == n ? l = 0 : nothing
        end
    end
    return p
end

# Test if two permutations commute
function permcommute(p::AbstractVector{T}, q::AbstractVector{T}) where T
    length(q) < length(p) ? (p,q) = (q,p) : nothing
    for i in length(p)
        q[p[i]] == p[q[i]] || return false
    end
    for i in length(p)+1:length(q)
        q[i] == i || return false
    end
    return true
end

# is m a permutation matrix
function isperm(m::AbstractMatrix{T}) where T
    sx,sy = size(m)
    sx == sy || return false
    z = zero(T)
    o = one(T)
    seen = falses(sx)
    for i in 1:sx
        for j in 1:sx
            val = m[j,i]
            if val == o
                seen[j] == true && return false
                seen[j] = true
            elseif val != z
                return false
            end
        end
    end
    true
end

# is cyclic decomposition a permutation
function isperm(cycs::AbstractVector{V}) where {V <: AbstractArray{T} where T}
    seen = counter(eltype(cycs[1])) # inefficient
    for c in cycs
        for i in c
            push!(seen,i) == 1 || return false
        end
    end
    return true
end

# is sparse representation a permutation
function isperm(sp::Dict{T, T}) where T
    sort(collect(keys(sp))) == sort(collect(values(sp)))  # inefficient
end

function isid(p::AbstractVector)
    for i in 1:length(p)
        i == p[i] || return false
    end
    return true
end

function permlistisequal(p::AbstractVector{T}, q::AbstractVector) where T
    lp = length(p)
    lq = length(q)
    s = one(T)
    if lp < lq
        for i in s:lp
            p[i] == q[i] || return false
        end
        for i in lp+s:lq
            q[i] == i || return false
        end
    else  # could factor code with refs, prbly not worth the trouble
        for i in s:lq
            p[i] == q[i] || return false
        end
        for i in lq+s:lp
            p[i] == i || return false
        end
    end
    return true
end

# is p < q with lexicographical ordering ?
for (f,ret) in ((:ltpermlist, :false), (:lepermlist, :true))
    @eval begin
function ($f)(p::AbstractVector{T}, q::AbstractVector{V}) where {T, V}
    o = one(T)
    lp = length(p)
    lq = length(q)
    if ( lp < lq )
        for k in o:lp
            if q[k] != p[k]
                return p[k] < q[k] ? true : false
            end
        end
        for k in lp+o:lq
            q[k] != k && return false
        end
    else
        for k in o:lq
            if q[k] != p[k]
                return p[k] < q[k] ? true : false
            end
        end
        for k in lq+o:lp
            p[k] != k && return false
        end
    end
    return $ret
end
end
end

# preimage of k under p
function preimage(p::AbstractVector, k::Int)
    k > length(p) && return k
    for i in 1:length(p)
        if p[i] == k
            return i
        end
    end
    return k  #  make preimage consistent with image
#    error("Can't find inverse image of $k.")
end

function preimage(p::Dict{T, T}, i::Int) where T
    for (k,v) in p
        v == i && return k
    end
    return i  #  make preimage consistent with image
#    error("Can't find inverse image of $k.")
end

# List of points mapped to same point by p and q
#function same(pin::PermList, qin::PermList)
function same(p::AbstractVector, q::AbstractVector)
    lp = length(p)
    lq = length(q)
    d = Array{eltype(p)}(undef, 0)
    if lp < lq
        for i in 1:lp
            p[i] == q[i] ? push!(d,p[i]) : nothing
        end
        for i in lp+1:lq
            q[i] == i ? push!(d,q[i]) : nothing
        end
    else  # could factor code with refs, prbly not worth the trouble
        for i in 1:lq
            p[i] == q[i] ? push!(d,p[i]) : nothing
        end
        for i in lq+1:lp
            p[i] != i ? push!(d,p[i]) : nothing
        end
    end
    return d
end

# The return type depends on value of input. How to get around this ?
# There is no integer Inf.
# agrees with gap (except definition, use of Inf is different)
function leastmoved(p::AbstractVector)
    lp = length(p)
    lm = lp+1
    for i in 1:lp
        k = p[i]
        k == i ? nothing :
           k < lm ? lm = k : nothing
    end
    return lm > lp ? Inf : lm
end

# agrees with gap
function greatestmoved(p::AbstractVector)
    lp = length(p)
    gm = 0
    for i in 1:lp
        k = p[i]
        k == i ? nothing :
           k > gm ? gm = k : nothing
    end
    return gm
end

function supportsize(p::AbstractVector)
    lp = length(p)
    count = 0
    for i in 1:lp
        k = p[i]
        k != i ? count += 1 : nothing
    end
    count
end

function support(p::AbstractVector)
    lp = length(p)
    mov = Array{eltype(p)}(undef, 0)
    for i in 1:lp
        k = p[i]
        k != i ? push!(mov,i) : nothing
    end
    return mov
end

function fixed(p::AbstractVector)
    lp = length(p)
    fixedel = Array{eltype(p)}(undef, 0)
    for i in 1:lp
        k = p[i]
        k == i ? push!(fixedel,i) : nothing
    end
    return fixedel
end

## Output ##

function cycleprint(io::IO, v::Array)
    len = length(v)
    print(io,"(")
    for i in 1:len
        print(io,v[i])
        i != len ? print(" ") : nothing
    end
    print(io,")")
end

function permarrprint(io::IO, v::Array)
    if isempty(v)
        print(io, "()")
    else
        pd = ndigits(maximum(v))
        len = length(v)
        print(io,"( ")
        for k in v
            print(io,rpad(k,pd,' ')," ")
        end
        print(io,")")
    end
end

cycleprint(v::Array) = cycleprint(STDOUT,v)
permarrprint(v::Array) = permarrprint(STDOUT,v)

# Copied from Base string.jl
function pprint_to_string(printfunc::Function, xs...)
    s = IOBuffer(Array{Uint8}(isa(xs[1],String) ? endof(xs[1]) : 0), true, true)
    truncate(s,0)
    for x in xs
        printfunc(s, x)
    end
    takebuf_string(s)
end

permarrstring(xs...) = pprint_to_string(permarrprint, xs...)
cyclestring(xs...) = pprint_to_string(cycleprint, x...)


include("matrixops.jl")
include("pivot.jl")

end # module Permplain
