abstract type AbstractAlignment end

###################################################################################
#################################### Alignment ####################################
###################################################################################

"""
    mutable struct Alignment{A,T} where {A, T<:Integer}

```
    data::Matrix{T}
    alphabet::Union{Nothing, Alphabet{A,T}}
    weights::Vector{Float64} = ones(size(dat,1))/size(dat,1) # phylogenetic weights of sequences
    names::Vector{String} = fill("", size(dat, 1))
```

Biological sequences as vectors of type `T<:Integer`.
`data` stores sequences in *columns*: `size(dat)` returns a tuple `(L, M)` with `L` the
length and `M` the number of sequences.
When displayed, shows `data` as an `MxL` matrix to match with traditional alignments.

`alphabet{A,T}` represents the mapping between integers in `data` and biological symbols of type `A` (nucleotides, amino acids...).
If `nothing`, the alignment cannot be mapped to biological sequences.

`weights` represent phylogenetic weights, and are initialized to `1/M`. They must sum to 1.
`names` are the label of sequences, and are expected to be in the same order as the columns
of `data`. They do not have to be unique, and can be ignored

**Important**: When built from a matrix, assumes that the sequences are stored in *columns*.

## Methods

- `getindex(X::Alignment, i)` returns a matrix/vector `X.data[:, i]`.
- `for s in X::Alignment` iterates over sequences.
- `eachsequence(X::Alignment)` returns an iterator over sequences (`Vector{Int}`).
- `eachsequence_weighted(X::Alignment)` returns an iterator over sequences and weights as tuples
- `subaln(X::Alignment, idx)` constructs the subaln defined by index `idx`.
"""
@kwdef mutable struct Alignment{A, T<:Integer} <: AbstractAlignment
    data::Matrix{T}
    alphabet::Union{Nothing, Alphabet{A,T}}
    weights::Vector{Float64} = ones(size(data, 2))/size(data, 2) # phylogenetic weights of sequences
    Meff::Float64 = size(data, 2)
    names::Vector{String} = fill("", size(data, 2))
    function Alignment{A,T}(data, alphabet, weights, Meff, ::Nothing) where {A, T}
        return Alignment{A,T}(data, alphabet, weights, Meff, fill("", size(data, 2)))
    end
    function Alignment{A,T}(data, alphabet, weights, Meff, names) where {A,T}
        @assert length(names) == length(weights) == size(data, 2) """\
            Inconsistent sizes between `data`, `weight` and `names` \
            - got $(size(data,2)), $(length(weights)), $(length(names))
            """

        # Check data and alphabet are consistent
        @assert isnothing(alphabet) || all(i -> in(i, alphabet), data) """\
            Some elements of `data` are not in `alphabet`
            Alphabet: $alphabet
            Problematic data: $(data[findall(x -> !in(x, alphabet), data)])
            """

        # Check weights
        @assert all(>(0), weights) "Weights must be positive"
        @assert isapprox(sum(weights), 1; rtol = 1e-8) """
            Weights must sum to 1 - got $(sum(weights))
            """
        if Meff - size(data, 2) > 1e-6
            error("Effective M larger than number of sequences $(Meff) > $(size(data,2))")
        elseif Meff - size(data, 2) > 0
            Meff = round(Int, Meff)
        end

        alphabet_copy = isnothing(alphabet) ? nothing : copy(alphabet)
        return new{A,T}(Matrix(data), alphabet_copy, copy(weights), Meff, string.(names))
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

`data` is a matrix of integers, with sequences stored in columns.
`alphabet` can be either
- an `Alphabet`
- `nothing`: no conversion from integers to biological symbols.
- something to build an alphabet from (*e.g.* a symbol like `:aa`, a string, ...).
    The constructor `Alphabet` will be called like so: `Alphabet(alphabet)`.

If the types of `alphabet` and `data` mismatch, `data` is converted.

`data` can also have the following shape:
- vector of integer vectors, *e.g.* [[1,2], [3,4]]: each element is considered as a sequence
- vector of integers: single sequence alignment
"""
function Alignment(data::AbstractMatrix{T}, alphabet::Alphabet{A,T}; kwargs...) where {A,T}
    return Alignment{A,T}(;data, alphabet, kwargs...)
end
function Alignment(data::AbstractMatrix{T}, ::Nothing; kwargs...) where T
    return Alignment{Nothing,T}(; data, alphabet=nothing, kwargs...)
end
function Alignment(D::AbstractMatrix{T}, alphabet::Alphabet{A,U}; kwargs...) where {A,T,U}
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

function Alignment(data::AbstractVector{T}, alphabet; kwargs...) where T<:Integer
    return Alignment(reshape(data, length(data), 1), alphabet; kwargs...)
end

"""
    Alignment(data::AbstractMatrix{T}; alphabet = :auto, kwargs...)

Keyword argument `alphabet` can be `:auto`, `:none`/`nothing`, or an input to the
constructor `Alphabet`.
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
        Alphabet(alphabet, T)
    end
    return Alignment(data, A; kwargs...)
end
function Alignment(data::AbstractVector{<:AbstractVector}; kwargs...)
    return Alignment(reduce(hcat, data); kwargs...)
end
function Alignment(data::AbstractVector{<:Integer}; kwargs...)
    return Alignment(reshape(data, length(data), 1); kwargs...)
end

################################################
##################### Misc #####################
################################################

function Base.show(io::IO, X::Alignment)
    L, M = size(X)
    print(io, "Alignment of M=$M sequences of length L=$L")
end
function Base.show(io::IO, x::MIME"text/plain", X::Alignment)
    L, M = size(X)
    println(io, "Alignment of M=$M sequences of length L=$L - Shown as `MxL` matrix")
    show(io, x, X.data')
end

function Base.copy(X::Alignment{A,T}) where {A,T}
    return Alignment{A,T}(;
        data = copy(X.data),
        alphabet = copy(X.alphabet),
        weights = copy(X.weights),
        names = copy(X.names),
    )
end

Base.convert(::Type{T}, X::Alignment{A,T}) where {A,T} = X
Base.convert(::Type{Alignment{A,T}}, X::Alignment{A,T}) where {A,T} = X
function Base.convert(::Type{I}, X::Alignment{A,J}) where {I<:Integer,A,J}
    return Alignment{A,I}(;
        data = convert(Matrix{I}, X.data),
        alphabet = convert(I, X.alphabet),
        weights = copy(X.weights),
        names = copy(X.names),
    )
end
function Base.convert(::Type{Alignment{A,T}}, X::Alignment) where {A,T<:Integer}
    return convert(T, X)
end


###################################################################################
################################# OneHotAlignment #################################
###################################################################################

#=
to implement...
=#
@kwdef mutable struct OneHotAlignment{A,T} <: AbstractAlignment
    data::OneHotArray{UInt32, 2, 3, Matrix{UInt32}}
    alphabet::Union{Nothing, Alphabet{A,T}}
    weights::Vector{Float64} = ones(size(data, 2))/size(data, 2) # phylogenetic weights of sequences
    names::Vector{String} = fill("", size(data, 2))

    function OneHotAlignment{A,T}(data, alphabet, weights, names) where {A,T}
        @assert length(names) == length(weights) == size(data, 3) """\
            Inconsistent sizes between `data`, `weight` and `names` \
            - got $(size(data,2)), $(length(weights)), $(length(names))
            """

        # Check data and alphabet are consistent
        @assert isnothing(alphabet) || all(i -> in(i, alphabet), 1:size(data,1)) """\
            Some elements of `data` are not in `alphabet`
            """

        # Check weights
        @assert all(>(0), weights) "Weights must be positive"
        @assert isapprox(sum(weights), 1; rtol = 1e-8) """
            Weights must sum to 1 - got $(sum(weights))
            """

        return new{A,T}(data, alphabet, weights, names)
    end
end

function onehot(X::Alignment{A,T}) where {A,T}
    return OneHotAlignment{A,T}(;
        data = onehotbatch(X.data, 1:length(Alphabet(X))),
        alphabet = Alphabet(X),
        weights = X.weights,
        names = X.names,
    )
end

###################################################################################
################################# AbstractAlignment ###############################
###################################################################################

################################################
############# Iterating / Indexing #############
################################################

"""
    size(A::AbstractAlignment)

Return a tuple with (in order) the length and the number of sequences.
"""
Base.size(aln::AbstractAlignment) = size(aln.data)
Base.size(aln::AbstractAlignment, dim) = size(aln.data, dim)
"""
    length(A::AbstractAlignment)

Return the number of sequences in `A`.
"""
Base.length(aln::AbstractAlignment) = size(aln.data, ndims(aln.data))

Base.iterate(X::AbstractAlignment) = iterate(eachslice(X.data, dims=ndims(X.data)))
function Base.iterate(X::AbstractAlignment, state)
    return iterate(eachslice(X.data, dims=ndims(X.data)), state)
end

"""
    eachsequence(X::AbstractAlignment[, indices]; skip)

Return an iterator over the sequences in `X`.
If `indices` is specified, consider only sequences at the corresponding indices.
Use the integer argument `skip` to return only one sequence every `skip` (~ `1:skip:end`).
"""
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

"""
    named_sequences(X::AbstractAlignment; skip)

Return an iterator of the form `(name, sequence)` over `X`.
"""
function named_sequences(X::AbstractAlignment, indices)
    return zip(X.names[indices], eachsequence(X, indices))
end
function named_sequences(X::AbstractAlignment; skip::Integer = 1)
    @assert skip > 0 "`skip` kwarg must be positive - instead $skip"
    return if skip == 1
        zip(X.names, eachslice(X.data, dims=ndims(X.data)))
    else
        indices = 1:skip:size(X.data)[end]
        zip(X.names[indices], eachsequence(X, indices))
    end
end

# Different for OneHot and normal alignment
Base.eltype(X::Alignment{A,T}) where {A,T} = AbstractVector{T}

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
    Y = Alignment(data_copy, copy(X.alphabet))
    Y.weights = X.weights[indices] / sum(X.weights[indices])
    Y.Meff = sum(X.weights[indices]) * X.Meff
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

"""
    find_sequence(label::AbstractString, aln::AbstractAlignment)

Find sequence with name `label` in `aln`, and return `(index, sequence)`.
Scales as the number of sequences.
Return the first sequence that matches the label.

!!! Return a *view* of the sequence.
"""
function find_sequence(label::AbstractString, aln::AbstractAlignment)
    i = findfirst(==(label), aln.names)
    return (i, isnothing(i) ? nothing : aln[i])
end
"""
    match_sequences(pattern, aln::AbstractAlignment)

Find sequences whose name matches `label` in `aln`, and return `(indices, sequences)`.
Sequences are returned as columns.

!!! Return a *view* of the sequences.
"""
function match_sequences(pattern, aln::AbstractAlignment)
    idx = findall(x -> occursin(pattern, x), aln.names)
    return idx, eachsequence(aln, idx)
end

"""
    filter(f, aln::AbstractAlignment)

Filter sequences of `aln` using boolean function `f`. Return another `Alignment`.
Use `filter(f, eachsequence(aln))` to obtain an array of sequences.
"""
function filter(f, aln::A) where A <: AbstractAlignment
    idx = findall(f, eachsequence(aln))

    return A(
        data = copy(aln.data[:, idx]),
        alphabet = aln.alphabet,
        names = aln.names[idx],
        weights = aln.weights[idx] / sum(aln.weights[idx]),
        Meff = sum(aln.weights[idx])*aln.Meff,
    )
end

#=============================================#
################ Concatenating ################
#=============================================#

function Base.cat(A::Alignment, B::Alignment, C::Vararg{<:Alignment})
    if !allequal(x -> x.alphabet, (A, B, C...))
        error("All alignments must share the same alphabet")
    end
    if !allequal(x -> size(x, 1), (A, B, C...))
        error("All alignments must have the same length")
    end

    data = hcat(A.data, B.data, map(x -> x.data, C)...)
    names = vcat(A.names, B.names, map(x -> x.names, C)...)
    weights = vcat(A.weights, B.weights, map(x -> x.weights, C)...)
    weights = weights / sum(weights)
    return Alignment(data, A.alphabet; names, weights)
end

#===========================#
########## Sorting ##########
#===========================#
"""
    sort!(aln::AbstractAlignment; kwargs...)

Sort `aln` in place.
The permutation is obtained by calling `sortperm(aln.names; kwargs...)`.
"""
function Base.sort!(aln::AbstractAlignment; with_name=true, with_seq=false, kwargs...)
    @argcheck with_name || with_seq """
    Expect either `with_name` or `with_seq`, but not both. Instead $with_name - $with_name.
    """
    return if with_name
        sort_by_name!(aln; kwargs...)
    else
        error("Not yet implemented - use `with_name=true` and `with_seq=false`.")
        # sort_by_seq!(aln; kwargs...)
    end
end

function sort_by_name!(aln; kwargs...)
    p = sortperm(aln.names; kwargs...)
    aln.data = aln.data[:, p]
    aln.names = aln.names[p]
    aln.weights = aln.weights[p]
    return aln
end

"""
    sort(aln::AbstractAlignment; kwargs...)

Copy `aln` and call `sort!`.
"""
Base.sort(aln::AbstractAlignment; kwargs...) = sort!(copy(aln); kwargs...)

#==================#
####### Misc #######
#==================#


Alphabet(A::AbstractAlignment) = A.alphabet

Base.unique(X::AbstractAlignment) = subsample(X, unique(i -> X[i], eachindex(X)))

sequence_length(X::Alignment) = size(X, 1)
sequence_length(X::OneHotAlignment) = size(X, 2)
sequence_number(X::AbstractAlignment) = last(size(X))

function Random.rand(rng::AbstractRNG, X::Random.SamplerTrivial{<:AbstractAlignment})
    M = sequence_number(X[])
    return X[][rand(rng, 1:M)]
end
