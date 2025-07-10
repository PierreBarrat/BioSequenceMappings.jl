"""
    hamming(x, y; normalize=true, positions=nothing)

Hamming distance between `Vectors` `x` and `y`.
Only sites in vector `positions` will be considered.
"""
function hamming(
    X::AbstractVector{<:Integer}, Y::AbstractVector{<:Integer};
    normalize=true, positions=nothing, exclude_state=nothing,
)
    @assert length(X) == length(Y) """Expect vectors of same length.
        Instead $(length(X)) != $(length(Y))"""

    positions = isnothing(positions) ? (1:length(X)) : positions
    H = 0
    Z = 0
    if isnothing(exclude_state)
        for i in positions
            Z += 1
            (X[i] != Y[i]) && (H += 1)
        end
    else
        for i in positions
            if X[i] != exclude_state && Y[i] != exclude_state
                Z += 1
                (X[i] != Y[i]) && (H += 1)
            end
        end
    end
    return normalize ? H/Z : H
end

function hamming(X::Alignment, i::Integer, Y::Alignment, j::Integer)
    error("Not implemented yet")
end


"""
    pairwise_hamming(X, Y; step=1, step_left, step_right, as_vec=true, kwargs...)
    pairwise_hamming(X; as_vec, step, kwargs...)

Return all hamming distances between sequences of `X` and `Y`.
In the second form, consider pairs of sequences in `X`.

Only consider sequences every `step`.
`step_left` and `step_right` can be used to skip sequence either in `X` or in `Y`.
This is useful for large alignment, as the number of computations grows with the product
    of the size of the alignments

By default, the return value is a vector organized like
`[H(1,2), H(1,3), ..., H(M-1, M)]` with `H` standing for hamming distance and `M` for the
number of sequences.
If a matrix is prefered, use `as_vec=false`

Extra keyword arguments are passed to `hamming`.
"""
function pairwise_hamming(
    X::AbstractAlignment, Y::AbstractAlignment;
    step=1, step_left=step, step_right=step, as_vec=true, kwargs...
)
    if Alphabet(X) != Alphabet(Y)
        @warn """Alignments do not have the same alphabet. Are you sure?
            Left aln: $(Alphabet(X))
            Right aln: $(Alphabet(Y))
        """
    end
    X_sequences = eachsequence(X; skip = step_left)
    Y_sequences = eachsequence(Y; skip = step_right)
    D = [hamming(x, y; kwargs...) for x in X_sequences, y in Y_sequences]
    return if as_vec
        M, N = size(D)
        [D[i,j] for i in 1:M for j in 1:N]
    else
        D
    end
end
function pairwise_hamming(X::AbstractAlignment; as_vec=true, step=1, kwargs...)
    n = 0
    M = size(X, 2)
    for i in 1:step:M, j in (i+1):step:M
        n += 1
    end

    if as_vec
        H = zeros(Float64, n)
        n = 1
        for i in 1:step:M, j in (i+1):step:M
            H[n] = hamming(X[i], X[j]; kwargs...)
            n += 1
        end
        return H
    else
        H = zeros(Float64, M, M)
        for i in 1:step:M, j in (i+1):step:M
            H[j,i] = hamming(X[i], X[j]; kwargs...)
            H[i,j] = H[j,i]
        end
        return H
    end
end
