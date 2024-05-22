function hamming(X::AbstractVector{<:Integer}, Y::AbstractVector{<:Integer}; normalize=true)
    @assert length(X) == length(Y) """Expect vectors of same length.
        Instead $(length(X)) != $(length(Y))"""
    H = sum(zip(X, Y)) do (x,y)
        x != y
    end
    return normalize ? H / length(X) : H
end


function hamming(X::Alignment, i::Integer, Y::Alignment, j::Integer)
    @warn "Not implemented yet"
    return nothing
end


"""
    pairwise_hamming(X, Y; step=1, step_left, step_right)

Return matrix of all hamming distances between sequences of `X` and `Y`.
Only consider sequences every `step`.
`step_left` and `step_right` can be used to skip sequence either in `X` or in `Y`:
"""
function pairwise_hamming(
    X::Alignment, Y::Alignment; step=1, step_left=step, step_right=step, normalize=true,
)
    if Alphabet(X) != Alphabet(Y)
        @warn """Alignments do not have the same alphabet. Are you sure?
            Left aln: $(Alphabet(X))
            Right aln: $(Alphabet(X))
        """
    end
    X_sequences = eachsequence(X; skip = step_left)
    Y_sequences = eachsequence(Y; skip = step_right)
    return [hamming(x, y; normalize) for x in X_sequences, y in Y_sequences]
end

"""
    pairwise_hamming(X; step, as_vec=true)

Vector of pairwise hamming distances of sequences in `X`, ordered as
`[H(1,2), H(1,3), ..., H(M-1, M)]` with `H` standing for hamming distance.
If `as_vec=false`, will return a `Matrix` instead.
"""
function pairwise_hamming(X::Alignment; step=1, normalize=true, as_vec = true)
    D = pairwise_hamming(X, X; step, normalize)
    return if as_vec
        M = size(D, 1)
        [D[i,j] for i in 1:M for j in (i+1):M]
    else
        D
    end
end
