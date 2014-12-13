module PermPlain

using DataStructures.counter
import Base: isperm, randperm

include("collect.jl")
include("util.jl")

export permlist, permcycles, permsparse, permmatrix # whether to export, and what ?

#  The following acronyms refer to the storage model, not the DataType.
#  Specifically, they are more-or-less plain julia types.
#  PLIST  means  permutation stored as in one line array form
#  PCYC   means  permutation stored as cyclic decomposition

randperm{T<:Real}(::Type{T}, n::Integer) = collect(T,randperm(n))

## permcycles. Find cyclic decomposition  ##

# Compute cyclic decomposition (PCYC) from input permutation list (PLIST).
# builds a cycle list in the canonical order.
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

permcycles{T<:Real}(m::AbstractArray{T,2}) = permcycles(permlist(m))

permcycles{T}(sp::Dict{T,T}) = sparsetocycles(sp)
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
        push!(cycs,cyc)
    end
    return cycs
end

# Convert cyclic decomposition to canonical form
# used by gap, Mma, and Arndt thesis.
# Note that Arndt uses a different order internally to store the cycles as a single list.

function canoncycles{T}(cycs::AbstractArray{Array{T,1},1})
    ocycs = Array(Array{T,1},0)
    for cyc in cycs
        push!(ocycs,circshift(cyc,-indmin(cyc)+1))
    end
    sort!(ocycs,by=(x)->x[1])
    return ocycs
end


## cyclelengths. Find cyclic decomposition, but only save cycle lengths ##

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
# The method permpower below is usually preferred because it does less allocation
function permpower2{T<:Real}(p::AbstractVector{T}, n::Integer)
    n == 0 && return [one(T):convert(T,length(p))]
    n == 1 && return copy(p) # for consistency, don't return ref
    n < 0  && return permpower2(invperm(p),-n)
    q = permpower2(p, int(floor(n/2)))
    q = q[q]
    return iseven(n) ? q : p[q]
end

function permpower!{T<:Real}(p::AbstractVector{T},
                             pret::AbstractVector{T},
                             ptmp::AbstractVector{T},                             
                             n::Integer)
    onep = one(T)
    lenp = convert(T,length(p))
    n == 0 && (for i in onep:lenp pret[i] = i end; return )
    n < 0  && (permpower!(invperm(p), pret, ptmp, -n); return )
    n == 1 && (copy!(pret,p); return)
    permpower!(p, ptmp, pret,int(floor(n/2)))
    if iseven(n)
        for i in onep:lenp pret[i] = ptmp[ptmp[i]] end
    else
        for i in onep:lenp pret[i] = p[ptmp[ptmp[i]]] end
    end
end

# This does less allocation (in general) than permpower2.
function permpower{T<:Real}(p::AbstractVector{T}, n::Integer)
    n == 0 && return [one(T):convert(T,length(p))]
    n == 1 && return copy(p) # for consistency, don't return ref    
    pret = similar(p)
    ptmp = similar(p)
    permpower!(p,pret,ptmp,n)
    return pret
end

# Compute power of permutation. Both input and output are PCYC
# see pari perm.c
# Careful of degeneracy, and empty array may be returned.
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
# This routine assumes cycles of length 1 are included. Still don't understand it.
#    for j = 1:length(cyc)
#        n += length(cyc[j])+1
#    end
    p = T[1:n] # wasteful
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

permmatrix{T<:Real}(p::AbstractVector{T}, sparse::Bool = false) = permtomat(p,sparse)
# Convert PLIST to PMAT
function permtomat{T<:Real}(p::AbstractVector{T}, sparse::Bool = false)
    n::T = length(p)
    A = sparse ? speye(T,n) : eye(T,n)
    return A[p,:]
end

function permmatrix{T}(sp::Dict{T,T})
    n = convert(T,maximum(sp)[1])
    ot = one(T)
    z = zero(T)
    m = eye(T,n,n);
    for (i,v) in sp
        m[i,i] = z
        m[i,v] = ot
    end
    return m
end
    
function permtomat{T<:Real}(p::AbstractVector{T}, sparse::Bool = false)
    n::T = length(p)
    A = sparse ? speye(T,n) : eye(T,n)
    return A[p,:]
end

function mattoperm!{T<:Real}(m::AbstractArray{T,2}, p)
    n = size(m)[1]
    maxk = zero(T)    
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

# convert matrix to perm in list form
permlist{T<:Real}(m::AbstractArray{T,2}) = mattoperm(m)
mattoperm{T<:Real}(m::AbstractArray{T,2}) = mattoperm!(m,Array(T,size(m)[1]))

permsparse{T<:Real}(m::AbstractArray{T,2}) = mattosparse(m)
function mattosparse{T<:Real}(m::AbstractArray{T,2})
    p = Dict{T,T}()
    return mattoperm!(m,p), maximum(p)[1]
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

# is sparse representation a permutation
function isperm{T}(sp::Dict{T,T})
    sort(collect(keys(sp))) == sort(collect(values(sp)))  # inefficient
end

function isid{T<:Real}(p::AbstractVector{T})
    for i in 1:length(p)
        i == p[i] || return false
    end
    return true
end

permlist{T<:Real}(cycs::AbstractArray{Array{T,1},1}, pmax::Real = 0) =  cycstoperm(cycs,pmax)
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

permsparse{T<:Real}(p::AbstractVector{T}) = listtosparse(p)
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

permsparse{T<:Real}(cycs::AbstractArray{Array{T,1},1}) = cycstosparse(cycs)
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

permlist{T}(sp::Dict{T,T}) = sparsetolist(sp::Dict{T,T})
function sparsetolist{T}(sp::Dict{T,T})
    p = [one(T):convert(T,maximum(sp)[1])]
    for (i,v) in sp
        p[i] = v
    end
    return p
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
