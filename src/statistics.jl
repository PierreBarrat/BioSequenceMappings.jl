function f1(A::AbstractAlignment, w::AbstractVector{Float64})
    @assert sum(w) ≈ 1 "Weights must sum to 1 - Instead $(sum(w))"
    L, M = size(A)
    q = if isnothing(Alphabet(A))
        maximum(A.data)
    else
        length(Alphabet(A))
    end
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
    site_specific_frequencies(X::AbstractAlignment[, weights=X.weights]; as_vec=false)

Return the site specific frequencies of `X`.
If `as_vec`, the result is a vector of length `Lxq`.
Otherwise, it is a matrix of `q` rows and `L` columns (default).
"""
function site_specific_frequencies(X, w=X.weights; as_vec=false)
    f = f1(X, w)
    if as_vec
        q, L = size(f)
        return reshape(f, q*L)
    end
    return f
end

function consensus(X::AbstractAlignment, w=X.weights)
    f = site_specific_frequencies(X, w; as_vec=false)
    cons = map(argmax, eachcol(f))
    return Alignment(cons, Alphabet(X); names = ["consensus"])
end

function f2(A::AbstractAlignment, w::AbstractVector)
    @assert sum(w) ≈ 1 "Weights must sum to 1 - Instead $(sum(w))"
    L, M = size(A)
    q = isnothing(Alphabet(A)) ? maximum(A.data) : length(Alphabet(A))
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

"""
    pairwise_frequencies(X::AbstractAlignment, w=X.weights; as_mat=false)

Return a `q x q x L x L` tensor.
The `(a, b, i, j)` element is the fraction of sequences for which we see `a` at position
`i` and `b` at position `j`.

If `as_mat=true`, will return a `qL x qL` matrix, with `q x q` blocks representing
correlations between two specific columns.
"""
function pairwise_frequencies(X::AbstractAlignment, w=X.weights; as_mat=false)
    f = f2(X, w)
    if as_mat
        q, L = size(f, 1), size(f, 3)
        f = permutedims(f, [1, 3, 2, 4])
        return reshape(f, q*L, q*L)
    end
    return f
end

"""
    pairwise_correlations(X, w=X.weights; as_mat=false)

Compute connected correlations: the difference between the pairwise frequencies and the
product of the single site frequencies.
See `?pairwise_frequencies` for the shape of the output.
"""
function pairwise_correlations(X::AbstractAlignment, w=X.weights; as_mat=false)
    f1 = site_specific_frequencies(X, w)
    f2 = pairwise_frequencies(X, w)

    q, L = size(f1)
    C = zeros(Float64, q, q, L, L)

    for i in 1:L, j in (i+1):L, a in 1:q, b in 1:q
        C[a, b, i, j] = f2[a, b, i, j] - f1[a, i]*f1[b, j]
        C[b, a, j, i] = C[a, b, i, j]
    end
    if as_mat
        C = permutedims(C, [1, 3, 2, 4])
        return reshape(C, q*L, q*L)
    end
    return C
end
