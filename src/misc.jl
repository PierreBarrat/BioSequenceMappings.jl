"""
    hamming(x, y; normalize=true, positions=nothing)

Hamming distance between `Vectors` `x` and `y`.
Only sites in vector `positions` will be considered.
"""
function hamming(
    X::AbstractVector{<:Integer}, Y::AbstractVector{<:Integer};
    normalize=true, positions = nothing,
)
    @assert length(X) == length(Y) """Expect vectors of same length.
        Instead $(length(X)) != $(length(Y))"""
    H = if isnothing(positions)
        sum(zip(X, Y)) do (x,y)
            x != y
        end
    else
        sum(positions) do i
            X[i] != Y[i]
        end
    end
    Z = (isnothing(positions) ? length(X) : length(positions))
    return normalize ? H / Z : H
end


function hamming(X::Alignment, i::Integer, Y::Alignment, j::Integer)
    @warn "Not implemented yet"
    return nothing
end


"""
    pairwise_hamming(X, Y; step=1, step_left, step_right, as_vec=true, kwargs...)
    pairwise_hamming(X; kwargs...)

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
        M = size(D, 1)
        [D[i,j] for i in 1:M for j in (i+1):M]
    else
        D
    end
end
pairwise_hamming(X::AbstractAlignment; kwargs...) = pairwise_hamming(X, X; kwargs...)
