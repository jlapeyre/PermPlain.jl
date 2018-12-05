# Copied from linalg/lu.jl
function ipiv2perm(v::PList{T}, maxi::Integer) where T
    p = T[1:maxi]
    @inbounds for i in 1:length(v)
        p[i], p[v[i]] = p[v[i]], p[i]
    end
    return p
end

ipiv2perm(v::PList) = ipiv2perm(v, length(v))

# A naive algorithm, maybe there exists a faster one.
# Build the pivot while permuting 1:n. Make
# inverse perm at the same time, so we know where to
# find p[i]
function perm2ipiv(p::PList)
    n = length(p)
    m = [1:n]
    v = similar(m)
    x = copy(m)
    @inbounds for i in 1:n
        j = x[p[i]]
        mi = m[i]
        mj = m[j]
        m[i],m[j] = mj,mi
        v[i] = j
        x[mj] = i
        x[mi] = j
    end
    return v
end
