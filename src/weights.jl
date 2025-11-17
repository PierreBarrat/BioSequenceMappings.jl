"""
    compute_weights(X::AbstractAlignment, θ = 0.2; normalize = true)

Compute phylogenetic correction weights for sequences of `X`.
    The weight sequence `S` is `1/N`, where `N` is the number of sequences in `X` at
    hamming distance less than `H` from `S` (including `S` itself).
    The threshold `H` is `floor(θ⋅L)` where `L` is the sequence length.


The return value is a tuple `(weights, Meff)`, where `Meff` is the sum of weights
(pre-normalization). If `normalize`, weights are normalized to sum to one. .
"""
function compute_weights(X::AbstractAlignment, θ = 0.2; normalize = true)
    return if θ == 0.
        M = sequence_number(X)
        w = ones(M)
        normalize ? (w/sum(w), M) : (w, M)
    else
        compute_weights(X.data, θ; normalize)
    end
end

"""
    compute_weights!(X, θ; kwargs...)

Compute and set weights for `X`. See `compute_weights`.
"""
function compute_weights!(X, θ = 0.2; kwargs...)
    w, Meff = compute_weights(X, θ; kwargs...)
    X.weights = w
    X.Meff = Meff
    return w, Meff
end


function compute_weights(
    Y::AbstractMatrix{<:Integer}, theta::Float64; normalize=true,
)
    L, M = size(Y)
    threshold = Int(floor(L * theta))
    weights = ones(Int, M);
    d = 0
    @inbounds @showprogress for m in 1:M, l in (m+1):M
        d = 0 # distance
        i = 1 # index in sequence
        while d < threshold && i <= L
            (Y[i,m] != Y[i,l]) && (d += 1)
            i += 1
        end
        if d < threshold
            weights[m]+=1
            weights[l]+=1
        end
    end
    weights = 1 ./ weights
    Meff = sum(weights)
    return (normalize ? weights / Meff : weights, Meff)
end
