## permkronecker ##

# Kronecker product for matrices induces a Kronecker product on permutations.
"""
    permkron(p, q)


Compute the Kronecker product of permutations `p` and `q`
induced by their matrix representations.
"""
function permkron(p::PList{T}, q::PList{S}) where {T, S}
    np = length(p); nq = length(q)
    dc = Array{promote_type(T, S)}(undef, np*nq)
@inbounds for i in 1:np
    for k in 1:nq
            dc[nq*(i-1) + k] = nq*(p[i]-1)+q[k]
        end
    end
    dc
end

# This is partly broken
function permkron(pt::PDict{T}, qt::PDict) where {T<:Real}
    dout = Dict{T, T}()
    p = getdata(pt); q = getdata(qt)
    (isempty(p) || isempty(q)) && error("Can't yet do sparse kronecker product with identity permutations")
    np = convert(T,maximum(p)[1])
    nq = convert(T,maximum(q)[1])
    maxk = zero(T)
    #    for (i,pi) in p
    for i in 1:np
        pi = get(p,i,zero(T))
        pi = pi == 0 ? i : pi
        #        for (k,qk) in q
        for k in 1:nq
            qk = get(q,k,zero(T))
            qk = qk == 0 ? k : qk
            ind1 = nq *(pi-1) + qk
            ind2 = nq * (i-1) + k
            ind1 == ind2 ? continue : nothing
            z = get(dout,ind2,zero(T))
            z != zero(T) && continue
            dout[ind2] = ind1
            ind2 > maxk ? maxk = ind2 : nothing
            ind1 > maxk ? maxk = ind1 : nothing  # should be superfluous
        end
    end
    return PDict(dout, permlength(pt) * permlength(qt))
end

# This a bit more efficient than full matrix multiplication. Cuts out one loop.
function permkron(a::PList{T}, b::PMatrix{S}) where {T<:Real,S<:Real}
    (nrowa, ncola) = (length(a),length(a))
    (nrowb, ncolb) = size(b)
    R = zeros(promote_type(T,S), nrowa * nrowb, ncola * ncolb)
    d = invperm(a)
    for j = 1:ncola, l = 1:ncolb
        soff = ncola * (j-1)
        i = d[j]
        roff = ncola * (i-1)
        for k = 1:nrowb
            R[roff+k,soff+l] = b[k,l]
        end
    end
    return R
end

# For testing. If there is a difference, it is small
# function permkron2{T<:Real,S<:Real}(a::AbstractVector{T}, b::AbstractMatrix{S})
#     (nrowa, ncola) = size(a)
#     (nrowb, ncolb) = size(a)
#     R = zeros(promote_type(T,S), nrowa * nrowb, ncola * ncolb)
#     d = a.data
#     for j = 1:ncola, l = 1:ncolb
#         soff = ncola * (j-1)
#         i = d[j]
#         roff = ncola * (i-1)
#         for k = 1:nrowb
#             R[soff+l,roff+k] = b[l,k]
#         end
#     end
#     R
# end
