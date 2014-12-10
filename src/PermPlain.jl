module PermPlain

using DataStructures.counter

# Permutations implemented without defining new types.
# The types, PermList, PermCycs, use this code.

#  The following acronyms refer to the storage model, not the DataType.
#  Specifically, they are more-or-less plain julia types.
#  PLIST  means  permutation stored as in one line array form
#  PCYC   means  permutation stored as cyclic decomposition

export permcycles, cyclelengths, permsgn, permorder,
       permcompose, permcompose!, permpower, permlisttomatrix,
       permtotrans, cycletype, permlistisequal, isperm,
       canoncycles, cycstoperm, cycleprint, permarrprint,
       cyc_pow_perm, permcommute, permdistance, permordercyc

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
function permsgn{T<:Real}(p::AbstractVector{T})
    clengths = cyclelengths(p)
    return permsgn_from_lengths(clengths)
end

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
function permlisttomatrix{T<:Real}(p::AbstractVector{T}, sparse::Bool = false)
    n::T = length(p)
    A = sparse ? speye(T,n) : eye(T,n)
    return A[p,:]
end

# Broken
# Just work directly on list form,
# convert permutation to equivalent array of transpositions.
# But we have to choose storage.
function permtotrans{T<:Real}(p::AbstractVector{T})
    cycs = permcycles(p)
    transes = Array(Array{T,1},0)    
    for cyc in cycs
        len = length(cyc)
        if len == 2
            push!(transes,cyc)
        elseif len > 2
            for i in 2:len
                push!(transes, [cyc[1],cyc[i]])
            end
        end
    end
    return transes
end

# Can't get a signature to match, so use generic for now
#function cycstoperm{T<:Real}(cycs::AbstractArray{AbstractVector{T},1})
# The input cycles must be disjoint.

# somewhat inefficient
function cycstoperm{T<:Real}(cycs::AbstractArray{Array{T,1},1}, pmax::Integer = 0)  
#function cycstoperm(cycs, pmax=0)  # somewhat inefficient
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

function isperm{T<:Real}(cycs::AbstractArray{Array{T,1},1})
    seen = counter(eltype(cycs[1]))
    for c in cycs
        for i in c
            push!(seen,i) == 1 || return false
        end
    end
    return true
end

function permlistisequal{T<:Real, V<:Real}(p::AbstractVector{T}, q::AbstractVector{V})
    lp = length(p)
    lq = length(q)
    if lp < lq
        for i in 1:lp
            p[i] == q[i] || return false
        end
        for i in lp+1:lq
            q[i] == i || return false
        end
    else  # could factor code with refs, prbly not worth the trouble
        for i in 1:lq
            p[i] == q[i] || return false
        end
        for i in lq+1:lp
            p[i] == i || return false
        end
    end
    return true
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
