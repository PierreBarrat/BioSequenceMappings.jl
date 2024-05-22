function f1(A::AbstractAlignment, w::AbstractVector{Float64})
    @assert sum(w) ≈ 1 "Weights must sum to 1 - Instead $(sum(w))"
    L, M = size(A)
    q = length(Alphabet(A))
    X = A.data' # M x L - for faster iteration

    f = zeros(Float64, L*q)
    for i in 1:L, m in 1:M
        f[(i-1)*q + X[m,i]] += w[m]
    end
    return f
end

function f1(X::AbstractAlignment)
    M = size(X, 2)
    return f1(X, ones(Float64, M)/M)
end
f1(X::AbstractAlignment, ::Nothing) = f1(X)
f1(X::AbstractAlignment, weights::AbstractVector) = f1(X, weights)
# if weights comes directly from `compute_weights`
function f1(X::AbstractAlignment, weights::Tuple{Any,Any})
    return f1(X, weights[1])
end

"""
    site_specific_frequencies(X::AbstractAlignment, weights=nothing)

Return the site specific frequencies of `X`.
If `as_vec` (default), the result is a vector of length `Lxq`.
Otherwise, it is a matrix of `q` rows and `L` columns
"""
function site_specific_frequencies(X, args...; as_vec=true)
    f = f1(X, args...)
    if !as_vec
        L = size(X, 1)
        q = length(Alphabet(X))
        return reshape(f, q, L)
    end
    return f
end

function consensus(X::AbstractAlignment, w=nothing)
    f = site_specific_frequencies(X, w; as_vec=false)
    cons = map(argmax, eachcol(f))
    return Alignment(cons, Alphabet(X); names = ["consensus"])
end
