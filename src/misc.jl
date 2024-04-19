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

function pairwise_hamming(X::Alignment; step=1, normalize=true, as_vec = true)
    D = pairwise_hamming(X, X; step, normalize)
    return if as_vec
        M = size(D, 1)
        [D[i,j] for i in 1:M for j in (i+1):M]
    else
        D
    end
end
