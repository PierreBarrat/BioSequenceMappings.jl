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
