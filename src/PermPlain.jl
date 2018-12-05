"""
    module PermPlain

Permuation types and operations.

Types representing permutations: `PList`, `PCycles`, `PMatrix`, `PDict`

Functions for converting among reprentations: `permlist`, `permcycles`, `permmatrix`, `permsparse`.

Functions: `isperm`, `permsgn`, `canoncycles`, `cycletype`, `cyclelengths`, `numcycles`, `permlength`,
      `permorder`, `permpower`, `permpower2`, `permcommute`, `permcompose`, `permcompose!`,
     `permkron`, `ltpermlist`
"""
module PermPlain

using DataStructures: counter
import LinearAlgebra
using SparseArrays: sparse
using Random
import Base: isperm, ==

include("collect.jl")

export permlist, permcycles, permsparse, permmatrix, canoncycles, permpower, permpower2,
    permsgn, cycletype
export cyclelengths, numcycles, permlength, permorder, getdata, cyc_pow_perm, permcompose, permcompose!,
    permcommute, permkron, ltpermlist
export PCycles, PList, PMatrix, PDict

randperm(::Type{T}, n::Integer) where T = collect(T, randperm(n))

eye(::Type{T}, n::Integer) where T =  Matrix{T}(LinearAlgebra.I, n, n)
speye(::Type{T}, n) where T = sparse(LinearAlgebra.I * one(T), n, n)

PCyclesVec{T} = AbstractVector{<:AbstractVector{T}} where {T <: Real}

"""
    PCycles{T}

Cycle representation of a permutation.
"""
mutable struct PCycles{T}
    data::PCyclesVec{T}
    len::Int
end

"""
    permlength(p)

The number of elements in the permutation.
"""
permlength(c::PCycles) = c.len
getdata(c::PCycles) = c.data
==(c1::PCycles, c2::PCycles) = c1.data == c2.data && c1.len == c2.len

"""
    PList{T} = AbstractVector{T} where {T <: Real}

List representation of a permtuation.
"""
PList{T} = AbstractVector{T} where {T <: Real}

permlength(p::PList) = length(p)

"""
    PMatrix{T} = AbstractMatrix{T} where {T <: Real}

Permutation matrix.
"""
PMatrix{T} = AbstractMatrix{T} where {T <: Real}
permlength(m::PMatrix) = size(m, 1)

PDictData{T} = Dict{T, T} where T

"""
    PDict{T}

Sparse representation of a permtuation as a `Dict`.
"""
mutable struct PDict{T}
    data::PDictData{T}
    len::Int
end

permlength(c::PDict) = c.len
getdata(c::PDict) = c.data
==(c1::PDict, c2::PDict) = c1.data == c2.data && c1.len == c2.len

# Compute cyclic decomposition.
# Builds a cycle list in the canonical order.
# See pari for possible efficiency improvement.
"""
    permcycles(p::AbstractVector{T}) where {T <: Real}

Return a cycle representation of the permuation `p`.

### Example
```jldoctest
julia> p = [10, 8, 5, 1, 9, 4, 6, 3, 2, 7];

julia> permcycles(p)
2-element Array{Array{Int64,1},1}:
 [1, 10, 7, 6, 4]
 [2, 8, 3, 5, 9]
```
"""
function permcycles(p::PList{T}) where T
    n = length(p)
    visited = falses(n)
    cycles = Vector{Vector{T}}(undef, 0)
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
    return PCycles(cycles, n)
end

"""
    permcycles(m::AbstractMatrix)

Return a cycle representation of the permutation `m`.

### Example
```jldoctest
julia> p = [6, 1, 5, 3, 4, 2];

julia> m = permmatrix(p)
6×6 Array{Int64,2}:
 0  0  0  0  0  1
 1  0  0  0  0  0
 0  0  0  0  1  0
 0  0  1  0  0  0
 0  0  0  1  0  0
 0  1  0  0  0  0

julia> permcycles(m)
2-element Array{Array{Int64,1},1}:
 [1, 6, 2]
 [3, 5, 4]
```
"""
permcycles(m::PMatrix) = permcycles(permlist(m))

"""
    permcycles(sp::Dict{T, T})

Return a cycle representation of the permutation `sp`.

### Example
```jldoctest
julia> p = [6, 1, 5, 3, 4, 2];

julia> sp = permsparse(p)
(Dict(4=>3,2=>1,3=>5,5=>4,6=>2,1=>6), 6)

julia> p = [10, 8, 5, 1, 9, 4, 6, 3, 2, 7];

julia> sp, maxn = permsparse(p)
(Dict(7=>6,9=>2,4=>1,10=>7,2=>8,3=>5,8=>3,5=>9,6=>4,1=>10…), 10)

julia> permcycles(sp)
2-element Array{Array{Int64,1},1}:
 [7, 6, 4, 1, 10]
 [9, 2, 8, 3, 5]
```
"""
permcycles(sp::PDict) = sparsetocycles(sp)

function sparsetocycles(dsp::PDict{T}) where T
    cycs = Array{Vector{T}}(undef, 0)
    permlength(dsp) == 0 && return PCycles(cycs, permlength(dsp))
    sp = getdata(dsp)
    ks = collect(keys(sp))
    n = length(ks)
    seen = Dict{T, Bool}()
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
    return PCycles(cycs, permlength(dsp))
end

"""
    permlist(m::PMatrix)

Return a list representation of the permutation `m`.

### Example
```jldoctest
julia> m = [0 0 0 0 0 1; 1 0 0 0 0 0; 0 0 0 0 1 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 1 0 0 0 0]
6×6 Array{Int64,2}:
 0  0  0  0  0  1
 1  0  0  0  0  0
 0  0  0  0  1  0
 0  0  1  0  0  0
 0  0  0  1  0  0
 0  1  0  0  0  0

julia> print(permlist(m))
[6, 1, 5, 3, 4, 2]
```
"""
permlist(m::PMatrix) = mattoperm(m)
mattoperm(m::PMatrix{T}) where T = mattoperm!(Vector{T}(undef, size(m, 1)), m)

function mattoperm!(p, m::PMatrix{T}) where T
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

"""
    permlist(cycs::PCycles, pmax::Real = 0)

Return the list representation of the permuatation `cycs`, where `pmax` is
the length of the permutation.

### Examples
```jldoctest
julia> pcyc = [[4, 3, 5], [2, 1, 6]]
2-element Array{Array{Int64,1},1}:
 [4, 3, 5]
 [2, 1, 6]

julia> print(permlist(pcyc))
[6, 1, 5, 3, 4, 2]
```
"""
permlist(cycs::PCycles) = _permlist(getdata(cycs), permlength(cycs))

_permlist(cycsvec::PCyclesVec, pmax::Real = 0) = cycstoperm(cycsvec, pmax)

function cycstoperm(cycsvec::PCyclesVec{T}, pmax::Integer = 0) where T
    isempty(cycsvec) && return collect(one(T):convert(T,pmax))
    cmaxes = [maximum(c) for c in cycsvec]
    cmax = maximum(cmaxes)  # must be a faster way
#    perm = [one(T): (pmax > cmax ? convert(T,pmax) : convert(T,cmax))]
    perm = collect(one(T): (pmax > cmax ? convert(T,pmax) : convert(T,cmax)))
@inbounds  for c in cycsvec
        for i in 2:length(c)
            perm[c[i-1]] = c[i]
        end
        perm[c[end]] = c[1]
    end
    perm[cmax+1:pmax] = collect(cmax+1:pmax)
    return perm
end

"""
    permlist(sp::PDict)

Return the list representation of the permuatation `cycs`.

### Examples
```jldoctest
julia> sp = Dict(4=>3,2=>1,3=>5,5=>4,6=>2,1=>6);

julia> print(permlist(sp))
[6, 1, 5, 3, 4, 2]
```
"""
permlist(sp::PDict) = sparsetolist(sp)

function sparsetolist(sp::PDict{T}) where T
    p = collect(one(T):convert(T, permlength(sp)))
    data = getdata(sp)
@inbounds for (i,v) in data
        p[i] = v
    end
    return p
end

## permsparse ##
"""
    permsparse(m::PMatrix)

Return a sparse (`PDict`) representation of the permutation `m`.
"""
permsparse(m::PMatrix) = mattosparse(m)

function mattosparse(m::PMatrix{T}) where T
    p = Dict{T, T}()
    return PDict{T}(mattoperm!(p, m), permlength(m))
end

permsparse(p::PList) = listtosparse(p)

function listtosparse(p::PList{T}) where T
    data = Dict{T,T}()
    maxk = zero(T)
    length(p) == 0 && return PDict{T}(data, maxk)
@inbounds for i in p
        pv = p[i]
        if pv != i
            data[i] = pv
            pv > maxk ? maxk = pv : nothing
        end
    end
    return PDict(data, permlength(p))
end

permsparse(cycs::PCycles) = cycstosparse(cycs)

function cycstosparse(cycs::PCycles{T}) where T
    data = Dict{T,T}()
    maxk = zero(T)
@inbounds for c in getdata(cycs)
        pv = c[1]
        data[c[end]] = pv
        c[end] > maxk ? maxk = c[end] : nothing
        for i in 1:length(c)-1
            pv = c[i]
            data[pv] = c[i+1]
            pv > maxk ? maxk = pv : nothing
        end
    end
    return PDict(data, permlength(cycs))
end

## permmatrix ##

"""
    permmatrix(p, sparse::Bool = false)

Convert `p` to a matrix representation.
"""
permmatrix(p::PList, sparse::Bool = false) = permtomat(p, sparse)
# Convert PLIST to PMAT
function permtomat(p::PList{T}, sparse::Bool = false) where T <: Real
    n::T = length(p)
    A = sparse ? speye(T, n) : eye(T, n)
    return A[p,:]
end

function permmatrix(sp::PDict{T}) where {T}
    n = convert(T, permlength(sp))
    ot = one(T)
    z = zero(T)
    m = eye(Int, n)
    data = getdata(sp)
@inbounds for (i, v) in data
        m[i,i] = z
        m[i,v] = ot
    end
    return m
end

# Convert cyclic decomposition to canonical form
# used by gap, Mma, and Arndt thesis.
# Note that Arndt uses a different order internally to store the cycles as a single list.
"""
    canoncycles(cycs::PCycles)

Convert cyclic decomposition `cycs` to a canonical form.
"""
function canoncycles(cycs::PCycles{T}) where T
    ocycsdata = Array{Array{T,1}}(undef, 0)
    for cyc in getdata(cycs)
        push!(ocycsdata, circshift(cyc, -argmin(cyc) + 1))
    end
    sort!(ocycsdata, by=(x)->x[1])
    return PCycles(ocycsdata, permlength(cycs))
end

## cyclelengths. Find cyclic decomposition, but only save cycle lengths ##

"""
    cyclelengths(p)

Return a `Vector` of the lengths of the cycles in the permutation `p`.

`cyclelengths(p)` is equivalent to, but more efficient than `lenth(permcycles(p))`.
"""
cyclelengths(c::PCycles) = [length(x) for x in getdata(c)]
cyclelengths(m::PMatrix) = cyclelengths(permlist(m))
function cyclelengths(p::PList{T}) where T
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

# Gives cyclelengths in canonical order.
function cyclelengths(spd::PDict{T}) where T
    sp = getdata(spd)
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

"""
    numcycles(p)

Return the number of cycles in the permutation `p`.

`numcycles(p)` is in general more efficient than
`length(cyclelengths(p))`.
"""
numcycles(p) = length(cyclelengths(p))
numcycles(c::PCycles) = length(getdata(c))
numcycles(m::PMatrix) = numcycles(permlist(m))
function numcycles(p::PList{T}) where T
    n = length(p)
    visited = falses(n)
    Ncycles = 0
@inbounds for k in one(T):convert(T, n)
        if ! visited[k]
            knext = k
            len = 0
            while ! visited[knext]
                len += 1
                visited[knext] = true
                knext = p[knext]
            end
            len > 1 ? Ncycles += 1 : nothing
        end
    end
    return Ncycles
end


# Wikipedia says some authors also require no fixed points.
# Here, we allow fixed points
"""
    iscyclic(p)

Return `true` if permutation `p` is cyclic.
"""
iscyclic(p) = numcycles(p) == 1

## cycle type, sign ##

"""
    cycletype(p)

The cycle type of permutation `p`.
"""
cycletype(p) = counter(cyclelengths(p))

permsgn_from_lengths(lens) = (-1)^(length(lens)+sum(lens))

"""
    permsgn(p)

The sign of permutation `p`.
"""
permsgn(p) = permsgn_from_lengths(cyclelengths(p))

function permorder_from_lengths(clengths)
    result = 1
    for c in clengths
        result = lcm(result, c)
    end
    return result
end

"""
    permorder(p)

Return the order of permutation `p`.
"""
permorder(p) = permorder_from_lengths(cyclelengths(p))

# distance between two PLISTs
# could use a macro for this and ==.
# is there a penalty for using swap macro on p and q instead of two branches ?

"""
    permdistance(p::PList, q::PList)

Compute the distance between permutations `p` and `q`.
"""
function permdistance(p::PList{T}, q::PList{T}) where T
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

"""
    permcompose(p::PList, q::PList)

Compute the composition of permutations `p` and `q`.
"""
function permcompose(p::PList{T}, q::PList) where {T}
    lq = length(q)
    lp = length(p)
    lq == lp && return p[q]  # Much faster to use the Base method
    lr = lq < lp ? lp : lq
    r = Array{T}(undef, lr)
    @inbounds  if lq <= lp
        r[1:lq] = p[q]
        r[lq+1:lp] = p[lq+1:lp]
    else
      @inbounds  for k in 1:lq
            i = q[k]
            r[k] = (i < lp ? p[i] : i)
        end
    end
    return r
end

function permcompose(q::PDict{T}, p::PDict{V}) where {T, V}
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
    return PDict(dout, max(permlength(q), permlength(p)))
end

permcompose(p::PMatrix, q::PMatrix) = p * q

## permapply ##

# Vector ?
function permapply(q::PDict{T}, a::AbstractArray) where T
    aout = copy(a)
    len = length(aout)
    for (k,v) in getdata(q)
        if k <= len && v <= len
            aout[k] = a[v]
        end
    end
    return aout
end

function permapply(q::PList, a::AbstractArray)
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

permapply(q::Union{PDict{T}, PList{T}}, a::String) where T = String(permapply(q, a.data))

# power of PLIST. output is PLIST
# The method permpower below is usually preferred because it does less allocation
function permpower2(p::PList{T}, n::Integer) where T
    n == 0 && return [one(T):convert(T,length(p))]
    n == 1 && return copy(p) # for consistency, don't return ref
    n < 0  && return permpower2(invperm(p),-n)
    q = permpower2(p, Int(floor(n/2)))
    q = q[q]
    return iseven(n) ? q : p[q]
end

function permpower!(p::PList{T},
                    pret::PList{T},
                    ptmp::PList{T},
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
function permpower(p::PList{T}, n::Integer) where T
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
function permpower(pcyc::PCycles{T}, exp::Integer) where {T}
    cyc = getdata(pcyc)
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
    return PCycles(c, permlength(pcyc))
end

permpower(q::PDict{T}, exp::Integer) where T = cycstosparse(permpower(sparsetocycles(q),exp))

# power of PCYC. output is PLIST
# see pari perm.c
function cyc_pow_perm(pcyc::PCycles{T}, exp::Integer) where {T}
    n = 0
    cyc = getdata(pcyc)
    # cmaxes = [maximum(c) for c in cyc]
    # n = maximum(cmaxes)
    n = permlength(pcyc)
# This routine assumes cycles of length 1 are included. Still don't understand it.
#    for j = 1:length(cyc)
#        n += length(cyc[j])+1
#    end
    p = collect(T, 1:n) # wasteful
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

"""
    permcommute(p, q)

Return `true` if permutations `p` and `q` commuted.
"""
function permcommute(p::PList, q::PList)
    length(q) < length(p) ? (p,q) = (q,p) : nothing
    for i in length(p)
        q[p[i]] == p[q[i]] || return false
    end
    for i in length(p)+1:length(q)
        q[i] == i || return false
    end
    return true
end

permcommute(p::PMatrix, q::PMatrix) = p * q == q * p

"""
    isperm(m::AbstractMatrix)

Return `true` if `m` is a permutation matrix.
"""
function isperm(m::PMatrix{T}) where T
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
function isperm(pcycs::PCycles)
    cycs = getdata(pcycs)
    isempty(cycs) && return true
    seen = counter(eltype(cycs[1])) # inefficient
    for c in cycs
        for i in c
            push!(seen,i) == 1 || return false
        end
    end
    return true
end

# is sparse representation a permutation
function isperm(dsp::PDict)
    sp = getdata(dsp)
    sort(collect(keys(sp))) == sort(collect(values(sp)))  # inefficient
end

# TODO: Could be isone ? No,..
function isid(p::PList)
    for i in 1:length(p)
        i == p[i] || return false
    end
    return true
end

function permlistisequal(p::PList{T}, q::PList) where T
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
function ($f)(p::PList{T}, q::PList{V}) where {T, V}
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
function preimage(p::PList, k::Int)
    k > length(p) && return k
    for i in 1:length(p)
        if p[i] == k
            return i
        end
    end
    return k  #  make preimage consistent with image
#    error("Can't find inverse image of $k.")
end

function preimage(p::PDict, i::Int)
    for (k,v) in getdata(p)
        v == i && return k
    end
    return i  #  make preimage consistent with image
#    error("Can't find inverse image of $k.")
end

# List of points mapped to same point by p and q
#function same(pin::PermList, qin::PermList)
function same(p::PList, q::PList)
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
function leastmoved(p::PList)
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
function greatestmoved(p::PList)
    lp = length(p)
    gm = 0
    for i in 1:lp
        k = p[i]
        k == i ? nothing :
           k > gm ? gm = k : nothing
    end
    return gm
end

function supportsize(p::PList)
    lp = length(p)
    count = 0
    for i in 1:lp
        k = p[i]
        k != i ? count += 1 : nothing
    end
    count
end

function support(p::PList)
    lp = length(p)
    mov = Array{eltype(p)}(undef, 0)
    for i in 1:lp
        k = p[i]
        k != i ? push!(mov,i) : nothing
    end
    return mov
end

function fixed(p::PList)
    lp = length(p)
    fixedel = Array{eltype(p)}(undef, 0)
    for i in 1:lp
        k = p[i]
        k == i ? push!(fixedel,i) : nothing
    end
    return fixedel
end

## Output ##

function cycleprint(io::IO, v::Vector)
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
