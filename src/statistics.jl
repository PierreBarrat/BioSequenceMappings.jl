function f1(A::AbstractAlignment, w::AbstractVector{Float64})
    @assert sum(w) ≈ 1 "Weights must sum to 1 - Instead $(sum(w))"
    L, M = size(A)
    q = length(Alphabet(A))
    X = A.data' # M x L - for faster iteration

    f = zeros(Float64, q, L)
    for i in 1:L, m in 1:M
        f[X[m,i], i] += w[m]
    end
    return f
end

function f1(X::AbstractAlignment)
    M = size(X, 2)
    return f1(X, ones(Float64, M)/M)
end
f1(X::AbstractAlignment, ::Nothing) = f1(X)
# f1(X::AbstractAlignment, weights::AbstractVector) = f1(X, weights) # pretty sure this does nothing?
# if weights comes directly from `compute_weights`
function f1(X::AbstractAlignment, weights::Tuple{Any,Any})
    return f1(X, weights[1])
end

"""
    site_specific_frequencies(X::AbstractAlignment, weights=nothing; as_vec=true)

Return the site specific frequencies of `X`.
If `as_vec` (default), the result is a vector of length `Lxq`.
Otherwise, it is a matrix of `q` rows and `L` columns
"""
function site_specific_frequencies(X, args...; as_vec=false)
    f = f1(X, args...)
    if as_vec
        q, L = size(f)
        return reshape(f, q*L)
    end
    return f
end

function consensus(X::AbstractAlignment, w=nothing)
    f = site_specific_frequencies(X, w; as_vec=false)
    cons = map(argmax, eachcol(f))
    return Alignment(cons, Alphabet(X); names = ["consensus"])
end

function f2(A::AbstractAlignment, w::AbstractVector)
    @assert sum(w) ≈ 1 "Weights must sum to 1 - Instead $(sum(w))"
    L, M = size(A)
    q = length(Alphabet(A))
    X = A.data'

    f = zeros(Float64, q, q, L, L)
    for i in 1:L, m in 1:M
        f[X[m,i], X[m,i], i, i] += w[m]
        for j in (i+1):L
            f[X[m, i], X[m, j], i, j] += w[m]
            f[X[m, j], X[m, i], j, i] += w[m]
        end
    end

    return f
end
f2(X) = f2(X, ones(Float64, size(X,2))/size(X,2))
f2(X, ::Nothing) = f2(X)
f2(X, weights::Tuple{Any,Any}) = f2(X, weights[1])

function pairwise_frequencies(X::AbstractAlignment, args...; as_mat=false)
    f = f2(X, args...)
    if as_mat
        q, L = size(f, 1), size(f, 3)
        permutedims!(f, [1, 3, 2, 4])
        return reshape(f, q*L, q*L)
    end
    return f
end
