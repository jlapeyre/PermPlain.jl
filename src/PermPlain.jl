module PermPlain

using DataStructures.counter
import Base: isperm, randperm

include("collect.jl")

# Permutations implemented without defining new types.
# The types, PermList, PermCycs, use this code.

#  The following acronyms refer to the storage model, not the DataType.
#  Specifically, they are more-or-less plain julia types.
#  PLIST  means  permutation stored as in one line array form
#  PCYC   means  permutation stored as cyclic decomposition

export permcycles, cyclelengths, permsgn, permorder,
       permcompose, permcompose!, permpower, permtomat, mattoperm,
       permtotrans, cycletype, permlistisequal, isperm,
       canoncycles, cycstoperm, cycleprint, permarrprint,
       cyc_pow_perm, permcommute, permdistance, permordercyc,
       ltpermlist

randperm{T<:Real}(::Type{T}, n::Integer) = collect(T,randperm(n))

# Get the lengths of the cycles in cyclic decomposition
# from input permutation list (PLIST).
# Store only the lengths of the cycles, not the cycles
# themselves. We assume that the elements in p can be used as indices
# into p. This increases efficiency twofold for Int32,
# if Int64 is the default. Not sure what it does for efficiency
# for FloatingPoint types.
function cyclelengths{T<:Real}(p::AbstractVector{T})
    n = length(p)
    visited = falses(n)
    lengths = Array(Int,0)
    for k in convert(T,1):convert(T,n)
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

cyclelengths{T<:Real}(c::AbstractArray{Array{T,1},1}) = [length(x) for x in c]

# compute the cycletype property from PLIST
cycletype{T<:Real}(p::AbstractVector{T}) = counter(cyclelengths(p))
cycletype{T<:Real}(c::AbstractArray{Array{T,1},1}) = counter(cyclelengths(c))
cycletype{T<:Real}(p::Dict{T,T}) = counter(cyclelengths(p))

# Compute cyclic decomposition (PCYC) from input permutation list (PLIST).
# This builds a cycle list in the canonical order.
# See pari for possible efficiency improvement.
function permcycles{T<:Real}(p::AbstractVector{T})
    n = length(p)
    visited = falses(n)
    cycles = Array(Array{T,1},0)
    for k in convert(T,1):convert(T,n)
        if ! visited[k]
            knext = k
            cycle = Array(T,0)
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

# Convert cyclic decomposition (PCYC) to canonical form
# canonical order used by gap, Mma, and Arndt thesis
# Note that Arndt uses a different order internally to store the cycles as a single list.
#function canoncycles{T<:Real}(cycs::AbstractArray{Array{T,1},1})
function canoncycles{T}(cycs::AbstractArray{Array{T,1},1})
    ocycs = Array(Array{T,1},0)
    for cyc in cycs
        push!(ocycs,circshift(cyc,-indmin(cyc)+1))
    end
    sort!(ocycs,by=(x)->x[1])
    return ocycs
end

permsgn_from_lengths(lens) = (-1)^(length(lens)+sum(lens))

# return the signature (also called sign) of the permutation
# from permutation list (PLIST)
permsgn{T<:Real}(p::AbstractVector{T}) =  permsgn_from_lengths(cyclelengths(p))
# from PCYC
permsgn{T<:Real}(c::AbstractArray{Array{T,1},1}) = permsgn_from_lengths(cyclelengths(c))

function permorder_from_lengths(clengths)
    result = 1
    for c in clengths
        result = lcm(result, c)
    end
    return result
end    

# return the order of the permutation from PLIST
permorder{T<:Real}(p::AbstractVector{T}) = permorder_from_lengths(cyclelengths(p))
permorder{T<:Real}(c::AbstractArray{Array{T,1},1}) = permorder_from_lengths(cyclelengths(c))
permorder{T<:Real}(p::Dict{T,T}) = permorder_from_lengths(cyclelengths(p))

macro swap!(p,q)
    return quote
        begin
            tmp = $(esc(p))
            $(esc(p)) = $(esc(q))
            $(esc(q)) = tmp
        end
    end
end

# Test if two PLISTs commute
function permcommute{T<:Real}(p::AbstractVector{T}, q::AbstractVector{T})
    length(q) < length(p) ? @swap!(p,q) : nothing
    for i in length(p)
        q[p[i]] == p[q[i]] || return false
    end
    for i in length(p)+1:length(q)
        q[i] == i || return false
    end
    return true
end

# distance between two PLISTs
# could use a macro for this and ==.
# is there a penalty for using swap macro on p and q instead of two branches ?
function permdistance{T<:Real}(p::AbstractVector{T}, q::AbstractVector{T})
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
function permcompose{T<:Real, V<:Real}(q::AbstractVector{T}, p::AbstractVector{V})
    lp = length(p)
    lq = length(q)
    lp == lq && return q[p]  # prbly not much time saved
    lr = lp < lq ? lq : lp
    r = Array(T,lr)
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
function permcompose!{T<:Real, V<:Real}(q::AbstractVector{T}, p::AbstractVector{V})
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

function permcompose{T<:Real, V<:Real}(q::Dict{T,T}, p::Dict{V,V})
    dout = Dict{T,T}()
    maxk = zero(T)    
    for (k,v) in p
        qv = get(q,v,zero(T))
        k == qv ? continue : nothing  # ignore 1-cycles
        dout[k] = (qv == zero(T) ? v : qv)
        qv > maxk ? maxk = qv : nothing
    end
    return dout, maxk
end

function permapply{T<:Real, V}(q::Dict{T,T}, a::AbstractArray{V})
    aout = copy(a)
    len = length(aout)
    for (k,v) in q
        if k <= len && v <= len
            aout[k] = a[v]
        end
    end
    aout
end

function permapply{T<:Real, V}(q::AbstractVector{T}, a::AbstractArray{V})
    aout = copy(a)
    lenq = length(q)
    lena = length(a)
    len =  lenq < lena ? lenq : lena
    for k in 1:len
        v = q[k]
        v <= len  && (aout[k] = a[v])
    end
    aout
end

function permapply{T<:Real, V<:String}(q::Union(Dict{T,T},AbstractVector{T}), a::V)
    ASCIIString(permapply(q,a.data))
end

# power of PLIST. output is PLIST
# This is slow. Does too much allocation.
function permpower{T<:Real}(p::AbstractVector{T}, n::Integer)
    n == 0 && return [one(T):convert(T,length(p))]
    n < 0  && return permpower(invperm(p),-n)
    n == 1 && return copy(p) # for consistency, don't return ref
    q = permpower(p, int(floor(n/2)))
    q = q[q]
    return iseven(n) ? q : p[q]
end

# Compute power of permutation.
# Both input and output are PCYC
# see pari perm.c
function permpower{T<:Real}(cyc::AbstractArray{Array{T,1},1}, exp::Integer)
    r = 1
    for j in 1:length(cyc)
        r += gcd(length(cyc[j])-1,exp)
    end
    c = Array(Array{T,1},0)
    for j in 1:length(cyc)
        v = cyc[j]
        n = length(v)
        e = mod(exp,n)
        g = gcd(n,e)
        m = div(n,g)
        if m == 1 continue end
        for i in 1:g
            p = Array(Int,0)
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

permpower{T<:Real}(q::Dict{T,T}, exp::Integer) = cycstosparse(permpower(sparsetocycles(q),exp))

# power of PCYC. output is PLIST
# see pari perm.c
function cyc_pow_perm{T<:Real}(cyc::AbstractArray{Array{T,1},1}, exp::Integer)
    n = 0
    cmaxes = [maximum(c) for c in cyc]
    n = maximum(cmaxes)
# This is for a cycles of length 1    
#    for j = 1:length(cyc)
#        n += length(cyc[j])+1
#    end
#    p = Array(Int,n)
    p = [1:n] # wasteful
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

# Convert PLIST to PMAT
# convert permutation to matrix operator
function permtomat{T<:Real}(p::AbstractVector{T}, sparse::Bool = false)
    n::T = length(p)
    A = sparse ? speye(T,n) : eye(T,n)
    return A[p,:]
end

# convert matrix to perm in list form
function mattoperm{T<:Real}(m::AbstractArray{T,2})
    n = size(m)[1]
    p = Array(T,n)
    for i in 1:n
        for j in 1:n
            if m[j,i] != 1
                continue
            end
            p[j] = i
        end
    end
    p
end

# is m a permutation matrix
function isperm{T<:Real}(m::AbstractArray{T,2})
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
function isperm{T<:Real}(cycs::AbstractArray{Array{T,1},1})
    seen = counter(eltype(cycs[1])) # inefficient
    for c in cycs
        for i in c
            push!(seen,i) == 1 || return false
        end
    end
    return true
end

# is sparse "permutation" a permutation
function isperm{T}(sp::Dict{T,T})
    sort(collect(keys(sp))) == sort(collect(values(sp)))  # inefficient
end

function isid{T<:Real}(p::AbstractVector{T})
    for i in 1:length(p)
        i == p[i] || return false
    end
    return true
end

# The input cycles must be disjoint.
# somewhat inefficient
function cycstoperm{T<:Real}(cycs::AbstractArray{Array{T,1},1}, pmax::Integer = 0)  
    length(cycs) == 0 && return [one(T):convert(T,pmax)]
    cmaxes = [maximum(c) for c in cycs]
    cmax = maximum(cmaxes)  # must be a faster way
    perm = [one(T): (pmax > cmax ? convert(T,pmax) : convert(T,cmax))]
    for c in cycs
        for i in convert(T,2):convert(T,length(c))
            perm[c[i-1]] = c[i]
        end
        perm[c[end]] = c[1]
    end
    perm[cmax+1:pmax] = [cmax+1:pmax]
    return perm
end

function permlistisequal{T<:Real, V<:Real}(p::AbstractVector{T}, q::AbstractVector{V})
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
function ($f){T<:Real, V<:Real}(p::AbstractVector{T}, q::AbstractVector{V})
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
function preimage{T<:Real}(p::AbstractVector{T}, k::Int)
    k > length(p) && return k
    for i in 1:length(p)
        if p[i] == k
            return i
        end
    end
    return k  #  make preimage consistent with image
#    error("Can't find inverse image of $k.")
end

function preimage{T<:Real}(p::Dict{T,T}, i::Int)
    for (k,v) in p
        v == i && return k
    end
    return i  #  make preimage consistent with image
#    error("Can't find inverse image of $k.")
end

# List of points mapped to same point by p and q
#function same(pin::PermList, qin::PermList)
function same{T<:Real, V<:Real}(p::AbstractVector{T}, q::AbstractVector{V})
    lp = length(p)
    lq = length(q)
    d = Array(eltype(p),0)
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
# agrees with gap (except definition,use of inifinity is different)
function leastmoved{T<:Real}(p::AbstractVector{T})
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
function greatestmoved{T<:Real}(p::AbstractVector{T})
    lp = length(p)
    gm = 0
    for i in 1:lp
        k = p[i]
        k == i ? nothing :
           k > gm ? gm = k : nothing
    end
    return gm
end

function supportsize{T<:Real}(p::AbstractVector{T})
    lp = length(p)
    count = 0
    for i in 1:lp
        k = p[i]
        k != i ? count += 1 : nothing
    end
    count
end

function support{T<:Real}(p::AbstractVector{T})
    lp = length(p)
    mov = Array(eltype(p),0)
    for i in 1:lp
        k = p[i]
        k != i ? push!(mov,i) : nothing
    end
    return mov
end

function fixed{T<:Real}(p::AbstractVector{T})
    lp = length(p)
    fixedel = Array(eltype(p),0)
    for i in 1:lp
        k = p[i]
        k == i ? push!(fixedel,i) : nothing
    end
    return fixedel
end

function listtosparse{T<:Real}(p::AbstractVector{T})
    data = Dict{T,T}()
    maxk = zero(T)
    length(p) == 0 && return (data,maxk)
    for i in p
        pv = p[i]
        if pv != i
            data[i] = pv
#            data[pv] = i  # correct order!
            pv > maxk ? maxk = pv : nothing
        end
    end
    return (data,maxk)
end

function cycstosparse{T<:Real}(cycs::AbstractArray{Array{T,1},1})
    data = Dict{T,T}()
    maxk = zero(T)
    for c in cycs
        pv = c[1]
        data[c[end]] = pv
#        pv > maxk ? maxk = pv : nothing
        c[end] > maxk ? maxk = c[end] : nothing        
        for i in 1:length(c)-1
            pv = c[i]            
            data[pv] = c[i+1]
            pv > maxk ? maxk = pv : nothing
        end
    end
    return (data,maxk)
end

function sparsetolist{T}(sp::Dict{T,T})
    p = [one(T):convert(T,maximum(sp)[1])]
    for (i,v) in sp
        p[i] = v
    end
    return p
end

function sparsetocycles{T}(sp::Dict{T,T})
    cycs = Array(Array{T,1},0)
    length(sp) == 0 && return cycs
    ks = collect(keys(sp))
    n = length(ks)
    seen = Dict{T,Bool}()
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
        didsee == false ? (push!(cycs,cyc); break) : nothing
        k1 = k
        cyc = Array(T,0)
        nseen = nseen + 1
        while true
            push!(cyc,k)
            seen[k] = true
            k = sp[k]
            nseen = nseen + 1
            if k == k1
                break
            end
        end
#        reverse!(cyc)   # inefficient
        push!(cycs,cyc)
    end
    return cycs
end

function cyclelengths{T}(sp::Dict{T,T})
    cyclens = Array(Int,0)
    ks = collect(keys(sp))
    n = length(ks)
    seen = Dict{T,Bool}()
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
        didsee == false ? (push!(cyclens,nincyc); break) : nothing
        k1 = k
        nincyc = 0
        nseen = nseen + 1
        while true
            nincyc += 1
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
    
end # module Permplain
