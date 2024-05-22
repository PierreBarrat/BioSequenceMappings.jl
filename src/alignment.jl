abstract type AbstractAlignment end

###################################################################################
#################################### Alignment ####################################
###################################################################################

"""
    mutable struct Alignment{T} where T <: Integer

```
    dat::Matrix{T}
    alphabet::Union{Nothing, Alphabet{T}}
    weights::Vector{Float64} = ones(size(dat,1))/size(dat,1) # phylogenetic weights of sequences
    names::Vector{String} = fill("", size(dat, 1))
```

Biological sequences as vectors of type `T<:Integer`.
`dat` stores sequences in *columns*: `size(dat)` returns a tuple `(L, M)` with `L` the
length and `M` the number of sequences.

`alphabet` represents the mapping between integers in `dat` and biological symbols (nucleotides, amino acids...).
If `nothing`, the alignment cannot be mapped to biological sequences.

`weights` represent phylogenetic weights, and are initialized to `1/M`. They must sum to 1.
`names` are the label of sequences, and are expected to be in the same order as the columns
of `dat`. They do not have to be unique, and can be ignored

# **Important**: When built from a matrix, will *transpose* the input; if `size(dat) = (M, L)`,
# `X=Alignment(dat, alphabet)` will return an object with `size(X.dat) = (L, M)`.
# In other words, assumes that the input matrix has sequences as rows.

## Methods

- `getindex(X::Alignment, i)` returns a matrix/vector `X.dat[:, i]`.
- `for s in X::Alignment` iterates over sequences.
- `eachsequence(X::Alignment)` returns an iterator over sequences (`Vector{Int}`).
- `eachsequence_weighted(X::Alignment)` returns an iterator over sequences and weights as tuples
- `subaln(X::Alignment, idx)` constructs the subaln defined by index `idx`.
"""
@kwdef mutable struct Alignment{T<:Integer} <: AbstractAlignment
    data::Matrix{T}
    alphabet::Union{Nothing, Alphabet{T}}
    weights::Vector{Float64} = ones(size(data, 2))/size(data, 2) # phylogenetic weights of sequences
    names::Vector{String} = fill("", size(data, 2))

    function Alignment{T}(data, alphabet, weights, names) where T
        @assert length(names) == length(weights) == size(data, 2) """\
            Inconsistent sizes between `data`, `weight` and `names` \
            - got $(size(data,2)), $(length(weights)), $(length(names))
            """

        # Check data and alphabet are consistent
        @assert isnothing(alphabet) || all(i -> in(i, alphabet), data) """\
            Some elements of `data` are not in `alphabet`
            """

        # Check weights
        @assert all(>(0), weights) "Weights must be positive"
        @assert isapprox(sum(weights), 1; rtol = 1e-8) """
            Weights must sum to 1 - got $(sum(weights))
            """

        return new{T}(Matrix(data), alphabet, weights, names)
    end
end

################################################
################# Constructors #################
################################################

#=
- from data matrix alone - kwargs to determine alphabet
- from data + alphabet
- from data + any alphabet constructor input (should be easy)
=#

function autofind_alphabet(data::AbstractMatrix{T}; verbose=true) where T <: Integer
    verbose && @info "Finding alphabet automatically from data ..."
    q = maximum(data)
    A = default_alphabet(q, T)
    verbose && @info "Found $A"
    return A
end


#=
Different options for alphabet
    1. an `Alphabet` of the right type --> main constructor
    2. `nothing` --> special case
    3. an `Alphabet` of the wrong type --> convert the data matrix if possible and fall back to 1
    4. an input `X` to the `Alphabet` constructor: string or symbol --> call `Alphabet(X, T)` and fall back to 1
=#
"""
    Alignment(data::AbstractMatrix, alphabet; kwargs...)

`data` is a matrix of integers.
`alphabet` can be either
    - an `Alphabet`
    - `nothing`
    - something to build an alphabet from: the constructor `Alphabet` will be called.
    (*e.g.* a symbol like `:aa`, a string, ...)

If the types of `alphabet` and `data` mismatch, `data` is converted.
"""
function Alignment(data::AbstractMatrix{T}, alphabet::Alphabet{T}; kwargs...) where T
    return Alignment{T}(;data, alphabet, kwargs...)
end
function Alignment(data::AbstractMatrix{T}, ::Nothing; kwargs...) where T
    return Alignment{T}(; data, alphabet=nothing, kwargs...)
end
function Alignment(D::AbstractMatrix{T}, alphabet::Alphabet{U}; kwargs...) where {T,U}
    data = convert(Matrix{U}, D)
    return Alignment(data, alphabet; kwargs...) # go to the first constructor
end
function Alignment(data::AbstractMatrix{T}, alphabet; kwargs...) where T
    return Alignment(data, Alphabet(alphabet, T); kwargs...) # go to first constructor
end
function Alignment(data::AbstractVector{<:AbstractVector{T}}, alphabet; kwargs...) where T
    # each element of `data` is one sequence
    return Alignment(reduce(hcat, data), alphabet; kwargs...)
end

function Alignment(data::AbstractVector{T}, alphabet::Alphabet{T}; kwargs...) where T
    return Alignment(reshape(data, length(data), 1), alphabet; kwargs...)
end

"""
    Alignment(data::AbstractMatrix{T}; alphabet = :auto, kwargs...)

Keyword argument `alphabet` can be `:auto`, `:none` or `nothing`, or an input to `Alphabet`.
Other keyword arguments are passed to the default constructor of `Alignment`.
"""
function Alignment(
    data::AbstractMatrix{T}; alphabet = :auto, verbose = true, kwargs...
) where T <: Integer
    A = if alphabet == :auto
        autofind_alphabet(data; verbose)
    elseif isnothing(alphabet) || in(alphabet, (:none, :no))
        nothing
    else
        alphabet(alphabet, T)
    end
    return Alignment(data, A; kwargs...)
end

################################################
##################### Misc #####################
################################################

function Base.show(io::IO, X::Alignment)
    L, M = size(X)
    print(io, "Alignment of M=$M sequences of length L=$L - Shown as `MxL` matrix")
    show(io, X.data')
end
function Base.show(io::IO, x::MIME"text/plain", X::Alignment)
    L, M = size(X)
    println(io, "Alignment of M=$M sequences of length L=$L - Shown as `MxL` matrix")
    show(io, x, X.data')
end

function Base.copy(X::Alignment{T}) where T
    return Alignment{T}(;
        data = copy(X.data),
        alphabet = copy(X.alphabet),
        weights = copy(X.weights),
        names = copy(X.names),
    )
end

Base.convert(::Type{T}, X::Alignment{T}) where T = X
Base.convert(::Type{Alignment{T}}, X::Alignment{T}) where T = X
function Base.convert(::Type{T}, X::Alignment) where T <: Integer
    return Alignment{T}(;
        data = convert(Matrix{T}, X.data),
        alphabet = convert(T, X.alphabet),
        weights = copy(X.weights),
        names = copy(X.names),
    )
end
function Base.convert(::Type{Alignment{T}}, X::Alignment) where T <: Integer
    return convert(T, X)
end


###################################################################################
################################# OneHotAlignment #################################
###################################################################################

#=
to implement...
=#
@kwdef mutable struct OneHotAlignment{T<:Integer} <: AbstractAlignment
    data::OneHotArray{T, 2, 3, Matrix{T}}
    alphabet::Union{Nothing, Alphabet{T}}
    weights::Vector{Float64} = ones(size(data, 2))/size(data, 2) # phylogenetic weights of sequences
    names::Vector{String} = fill("", size(data, 2))

    function OneHotAlignment{T}(data, alphabet, weights, names) where T
        @assert length(names) == length(weights) == size(data, 2) """\
            Inconsistent sizes between `data`, `weight` and `names` \
            - got $(size(data,2)), $(length(weights)), $(length(names))
            """

        # Check data and alphabet are consistent
        @assert isnothing(alphabet) || all(i -> in(i, alphabet), data) """\
            Some elements of `data` are not in `alphabet`
            """

        # Check weights
        @assert all(>(0), weights) "Weights must be positive"
        @assert isapprox(sum(weights), 1; rtol = 1e-8) """
            Weights must sum to 1 - got $(sum(weights))
            """

        return new{T}(Matrix(data), alphabet, weights, names)
    end
end

###################################################################################
################################# AbstractAlignment ###############################
###################################################################################

################################################
############# Iterating / Indexing #############
################################################

Base.size(aln::AbstractAlignment) = size(aln.data)
Base.size(aln::AbstractAlignment, dim) = size(aln.data, dim)
Base.length(aln::AbstractAlignment) = size(aln.data, ndims(aln.data))

Base.iterate(X::AbstractAlignment) = iterate(eachslice(X.data, dims=ndims(X.data)))
function Base.iterate(X::AbstractAlignment, state)
    return iterate(eachslice(X.data, dims=ndims(X.data)), state)
end

function eachsequence(X::AbstractAlignment, indices)
    return Iterators.map(i -> selectdim(X.data, ndims(X.data), i), indices)
end
function eachsequence(X::AbstractAlignment; skip::Integer = 1)
    @assert skip > 0 "`skip` kwarg must be positive - instead $skip"
    return if skip == 1
        eachslice(X.data, dims=ndims(X.data))
    else
        eachsequence(X, 1:skip:size(X.data)[end])
    end
end


# Different for OneHot and normal alignment
Base.eltype(X::Alignment{T}) where T = AbstractVector{T}

function Base.iterate(rX::Iterators.Reverse{<:AbstractAlignment})
    iterate(Iterators.Reverse(eachslice(rX.itr.data, dims = ndims(rX.itr.data))))
end
function Base.iterate(rX::Iterators.Reverse{<:AbstractAlignment}, state)
    iterate(Iterators.Reverse(eachslice(rX.itr.data, dims = ndims(rX.itr.data))), state)
end

Base.getindex(X::AbstractAlignment, i) = selectdim(X.data, ndims(X.data), i) # returns a view!!!
Base.firstindex(X::AbstractAlignment) = 1
Base.lastindex(X::AbstractAlignment) = length(X)
Base.view(X::AbstractAlignment, i) = getindex(X, i)
Base.keys(X::AbstractAlignment) = LinearIndices(1:length(X))

"""
    subsample(X::AbstractAlignment, indices)

Return an `Alignment` containing only the sequences of `X` at `indices`.
"""
subsample(X::AbstractAlignment, i::Int) = subsample(X, i:i)
function subsample(X::AbstractAlignment, indices)
    data_copy = copy(X[indices])
    Y = copy(X)
    Y.data = data_copy
    Y.weights = X.weights[indices] / sum(X.weights[indices])
    Y.names = X.names[indices]
    return Y
end

"""
    subsample_random(X::AbstractAlignment, m::Int)

Return an `Alignment` with `m` sequences taking randomly from `X`.
Sampling is done without replacement, meaning the `m` sequences are all at different
positions in `X`.
"""
function subsample_random(X::AbstractAlignment, m::Int)
    M = length(X)
    @assert m < M "Cannot take $m different sequences from alignment of size $M"
    return subsample(X, randperm(M)[1:m])
end

################################################
##################### Misc #####################
################################################

Alphabet(A::AbstractAlignment) = A.alphabet

Base.unique(X::AbstractAlignment) = subsample(X, unique(i -> X[i], eachindex(X)))

sequence_length(X::Alignment) = size(X, 1)
sequence_length(X::OneHotAlignment) = size(X, 2)
sequence_number(X::AbstractAlignment) = last(size(X))

function Random.rand(rng::AbstractRNG, X::Random.SamplerTrivial{<:AbstractAlignment})
    M = sequence_number(X[])
    return X[][rand(rng, 1:M)]
end
