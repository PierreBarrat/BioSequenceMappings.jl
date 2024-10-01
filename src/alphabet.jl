#=
A mapping is a `Dict{A, <:Integer}` with `A` the type of symbols used (usually `Char`).
For example, `Dict('A' => 1, 'C' => 2, ...)
Below is a utility function for mappings
=#
"""
    compute_mapping(s::AbstractString)

Return a `Dict{Int, Char}`: `Dict(i => c for (i,c) in enumerate(s))`.
"""
function compute_mapping(s::AbstractVector{A}, ::Type{I}=Int) where {A,I<:Integer}
    @assert allunique(s) "Cannot compute mapping from vector with duplicates - got $s"
    return Dict{A, I}(c => i for (i, c) in enumerate(s))
end
function reverse_mapping(mapping::Dict{T,U}) where {T, U}
    return Dict{U,T}(y => x for (x, y) in mapping)
end

#######################################################################################
####################################### Alphabet ######################################
#######################################################################################


"""
    struct Alphabet{A,I}
        characters::Vector{A}
        char_to_index::Dict{A, I}
        index_to_char::Dict{I, A}
        default_char = nothing
        default_index
    end

Structure allowing the mapping from biological symbols of type `A` to integers of type `I`.
    The typical use case would be `Alphabet{Char, Int}`.
`Alphabet` can be constructed
- from a `Vector` of symbols and an optional type `I`, *e.g.* `Alphabet(['A','C','G','T'], UInt8)::Alphabet{Char, UInt8}`
- from a `String` and an optional type, *e.g.* `Alphabet("ACGT")`
- from a mapping `Dict{A, I}` where `I<:Integer`: `Alphabet(Dict('A'=>1, 'C'=>2))`
- from a `Symbol`, using default alphabets, *e.g.* `Alphabet(:nt)`
- from an integer, using default alphabets (see `?default_alphabets`).
"""
@kwdef struct Alphabet{A, I<:Integer}
    characters::Vector{A} # Alphabet characters
    char_to_index::Dict{A, I} = compute_mapping(characters)
    index_to_char::Dict{I, A} = reverse_mapping(char_to_index)
    default_char::Union{Nothing, A} = nothing
    default_index::Union{Nothing, I} = _default_index_from_char(default_char, characters)

    # Constructor with checks
    function Alphabet{A,I}(
        characters, char_to_index, index_to_char, default_char, default_index,
    ) where {A,I}
        @assert isconcretetype(I) && isconcretetype(A) """\
            Use concrete type to build `Alphabet` - got `($A,$I)`
            """
        if !isbitstype(A)
            @warn "Building alphabet using non bit-type symbols $(A). \
            May lead to some issues, for instance when copying."
        end
        @assert !(A <: Integer) """\
            Symbol type should not be `Integer`. Got $characters of type $A
            """
        @assert length(characters) == length(char_to_index) == length(index_to_char) """\
            Inconsistent lengths between `characters`, `char_to_index` and `index_to_char`:\
            $(length(characters)) vs $(length(char_to_index)) vs $(length(index_to_char))
            """

        # are characters and char_to_index consistent
        @assert all(c -> haskey(char_to_index, c), characters) """\
            Incomplete `char_to_index`: some symbols in $characters not in $char_to_index
            """
        @assert all(x -> char_to_index[x[2]] == x[1], enumerate(characters)) """\
            Inconsistency between `characters` and `char_to_index`: \
            $characters -- $(char_to_index)
            """

        # are char_to_index and index_to_char consistent
        @assert all(c -> haskey(char_to_index, c), values(index_to_char)) """\
            Some symbols in `index_to_char` not in `char_to_index`
            """
        @assert all(i -> haskey(index_to_char, i), values(char_to_index)) """\
            Some indices in `char_to_index` not in `index_to_char`
            """
        @assert all(c -> c == index_to_char[char_to_index[c]], keys(char_to_index)) """\
            Inconsistent `char_to_index` and `index_to_char`.
            """

        # are defaults ok
        @assert !xor(isnothing(default_char), isnothing(default_index)) """\
            Got `default_char=$(default_char)` and `default_index=$(default_index)`.
            These should either both be `nothing`, or both something.
            This error can happen if the proposed `default_char` is not in the alphabet.
            """
        @assert isnothing(default_char) || in(default_char, characters) """\
            Default char $(default_char) not found in alphabet characters $(string)
            Maybe use `default_char = c` with `c` in your alphabet?
            """
        @assert isnothing(default_index) || haskey(index_to_char, default_index) """\
            Got `default_index=$(default_index)` \
            but alphabet `characters` of length $(length(characters))
            """
        @assert (isnothing(default_index) && isnothing(default_char)) ||
            index_to_char[default_index] == default_char """\
            Defaults are inconsistent with mapping.\
            `default_index=$(default_index)` maps to $(char_to_index[default_index]),\
             but `default_char=$(default_char)`
            """

        # Checking for duplicates
        @assert allunique(characters) """Cannot form alphabet from `characters` with duplicates\
            - got $characters
            """
        return new{A,I}(characters, char_to_index, index_to_char, default_char, default_index)
    end
end

_default_index_from_char(c, characters::AbstractVector) = findfirst(==(c), characters)
_default_index_from_char(::Nothing, ::AbstractVector) = nothing

################################################
############## Equality, copying ###############
################################################

function Base.:(==)(α::Alphabet{A,I}, β::Alphabet{B,J}) where {A,B,I,J}
    return A == B &&
        I == J &&
        α.characters == β.characters &&
        α.default_char == β.default_char
        α.default_index == β.default_index
        α.index_to_char == β.index_to_char &&
        α.char_to_index == β.char_to_index
end

function Base.hash(alphabet::Alphabet{A,I}, h::UInt) where {A,I}
    h = hash(Alphabet{A,I}, h)
    h = hash(alphabet.characters, h)
    h = hash(alphabet.index_to_char, h)
    h = hash(alphabet.char_to_index, h)
    h = hash(alphabet.default_char, h)
    h = hash(alphabet.default_index, h)
    return h
end

function Base.copy(alphabet::Alphabet{A,I}) where {A,I}
    return Alphabet{A,I}(;
        characters = copy(alphabet.characters), # assumes A is bit type, otherwise deepcopy
        char_to_index = copy(alphabet.char_to_index),
        index_to_char = copy(alphabet.index_to_char),
        default_char = alphabet.default_char,
        default_index = alphabet.default_index,
    )
end

Base.convert(::Type{I}, alphabet::Alphabet{A,I}) where {A,I<:Integer} = alphabet
Base.convert(::Type{Alphabet{A,I}}, alphabet::Alphabet{A,I}) where {A,I} = alphabet
function Base.convert(::Type{J}, α::Alphabet{A,I}) where {A,I,J<:Integer}
    return Alphabet{A,J}(;
        characters = copy(α.characters),
        char_to_index = convert(Dict{A, J}, α.char_to_index),
        index_to_char = convert(Dict{J, A}, α.index_to_char),
        default_char = α.default_char,
        default_index = isnothing(α.default_index) ? nothing : convert(J, α.default_index)
    )
end
function Base.convert(::Type{Alphabet{A,J}}, alphabet::Alphabet{A,I}) where {A,I,J<:Integer}
    return convert(J, alphabet)
end

################################################
################# Constructors #################
################################################


function Alphabet(
    characters::AbstractVector{<:A}, ::Type{I}=Int; kwargs...
) where {A, I<:Integer}
    return Alphabet{A,I}(;
        characters,
        char_to_index = compute_mapping(characters, I),
        kwargs...
    )
end
function Alphabet(S::AbstractString, ::Type{I}=Int; kwargs...) where I
    return Alphabet(collect(S), I; kwargs...)
end
Alphabet(A::Alphabet) = A
Alphabet(alphabet::Alphabet{A,I}, ::Type{J}) where {A,I,J} = convert(J, alphabet)
function Alphabet(characters::AbstractVector{A}, mapping::Dict{A,I}; kwargs...) where {A,I}
    return Alphabet{A,I}(;
        characters,
        char_to_index = mapping,
        kwargs...
    )
end
function Alphabet(S::AbstractString, mapping::Dict{<:AbstractChar, I}; kwargs...) where I
    return Alphabet(collect(S), mapping; kwargs...)
end
function Alphabet(char_to_index::AbstractDict{A, I}; kwargs...) where {A,I}
    # cannot do collect(keys(...)) because of undefined ordering of Dict
    characters = Vector{A}(undef, length(char_to_index))
    for (c, i) in char_to_index
        characters[i] = c
    end
    return Alphabet{A,I}(; characters, char_to_index, kwargs...)
end

################################################
################### Defaults ###################
################################################

const _DEFAULT_AA_ALPHABET_STRING = "-ACDEFGHIKLMNPQRSTVWY"
const _DEFAULT_NT_ALPHABET_STRING = "-ACGT"
const _DEFAULT_NT_ALPHABET_STRING_NOGAP = "ACGT"
const _DEFAULT_BINARY_ALPHABET_STRING = "01"

const _DEFAULT_ALPHABET_STRING = _DEFAULT_AA_ALPHABET_STRING

const aa_alphabet = Alphabet(_DEFAULT_AA_ALPHABET_STRING)
const aa_alphabet_names = (:aa, :AA, :aminoacids, :amino_acids)

const nt_alphabet = Alphabet(_DEFAULT_NT_ALPHABET_STRING)
const nt_alphabet_names = (:nt, :nucleotide, :dna)

const binary_alphabet = Alphabet(_DEFAULT_BINARY_ALPHABET_STRING)
const binary_alphabet_names = (:binary, :spin)

function Alphabet(name::Symbol, ::Type{T}=Int) where T <: Integer
    return if name in aa_alphabet_names
        # Alphabet(_DEFAULT_AA_ALPHABET_STRING, T; default_char = '-')
        convert(Alphabet{Char, T}, aa_alphabet)
    elseif name in nt_alphabet_names
        # Alphabet(_DEFAULT_NT_ALPHABET_STRING, T; default_char = '-')
        convert(Alphabet{Char, T}, nt_alphabet)
    elseif name in binary_alphabet_names
        # Alphabet(_DEFAULT_BINARY_ALPHABET_STRING, T)
        convert(Alphabet{Char, T}, binary_alphabet)
    else
        names = convert(
            Vector{Any}, vcat(aa_alphabet_names, nt_alphabet_names, binary_alphabet_names)
        )
        error("Unrecognized alphabet name $name - Possible names $names")
    end
end

"""
    default_alphabet(q::Int, T::Type)

- if `q==21`, amino acids
- if `q==5`, nucleotides
- if `q==4`, nucleotides without gaps
- if `q==2`, binary (0, 1)
- else, if `q<21`, return the restriction of amino acids to the first q sites
- if `q>21`, fails
"""
function default_alphabet(q::Integer, ::Type{T}=Int) where T <: Integer
    @assert q > 0 "`q` must be strictly positive - got $q"
    return if q == 2
        Alphabet(:binary, T)
    elseif q == 4
        Alphabet(_DEFAULT_NT_ALPHABET_STRING_NOGAP, T)
    elseif q == 5
        Alphabet(:dna, T)
    elseif 5 < q <= 21
        Alphabet(:aa, T)
    else
        error("No defined default alphabet for q = $q (<2 or > 21) - provide your own or use `nothing`")
    end
end

################################################
############ Transforming sequences ############
################################################

## Symbol to Int
function (alphabet::Alphabet{A,T})(c::A) where {A,T}
    i = get(alphabet.char_to_index, c, alphabet.default_index)
    if isnothing(i)
        error("Symbol $c not in alphabet, and no defaults set.")
    end
    return i
end
function (alphabet::Alphabet{A,T})(S::AbstractVector{<:A}) where {A,T}
    return map(x -> alphabet(x), S)
end
function (alphabet::Alphabet{Char,T})(S::AbstractString) where T
    return map(x -> alphabet(x), collect(S))
end

## Int to Symbol
function (alphabet::Alphabet)(x::Integer)
    c = get(alphabet.index_to_char, x, alphabet.default_char)
    if isnothing(c)
        error("$x is not in alphabet range, and no defaults set.")
    end
    return c
end
# If it makes sense, convert a `Vector{Int}` to a string
function (alphabet::Alphabet{<:AbstractChar,T} where T<:Integer)(X::AbstractVector{<:Integer})
    return prod(map(alphabet, X))
end
# Otherwise return an array of symbols
(alphabet::Alphabet)(X::AbstractVector{<:Integer}) = map(alphabet, X)

(alphabet::Alphabet)(::Missing) = missing

# Equivalent to alphabet(X) for Char based alphabets
# Otherwise, can be overloaded for the desired effect
# Used when writing fasta
function to_string(
    X::AbstractVector{<:Integer},
    alphabet::Alphabet{<:AbstractChar,<:Integer},
)
    return alphabet(X)
end

"""
    translate(x, original_alphabet::Alphabet, new_alphabet::Alphabet)

Return the translation in `new_alphabet` of an integer or a vector of integers `x` that is
expressed in `original_alphabet`.
"""
function translate(x::Integer, original_alphabet::Alphabet, new_alphabet::Alphabet)
    return x |> original_alphabet |> new_alphabet
end
function translate(X::AbstractVector{<:Integer}, A::Alphabet, B::Alphabet)
    return map(x -> translate(x, A, B), X)
end

################################################
##################### Misc #####################
################################################

Base.length(alphabet::Alphabet) = length(alphabet.characters)

Base.in(i::Integer, alphabet::Alphabet) = haskey(alphabet.index_to_char, i)
Base.in(c::A, alphabet::Alphabet{A,I}) where {A,I} = haskey(alphabet.char_to_index, c)

function name(alphabet::Alphabet{Char,I}) where I
    return if prod(alphabet.characters) == _DEFAULT_AA_ALPHABET_STRING
        :aa
    elseif prod(alphabet.characters) == _DEFAULT_NT_ALPHABET_STRING
        :dna
    elseif prod(alphabet.characters) == _DEFAULT_BINARY_ALPHABET_STRING
        :binary
    else
        :custom
    end
end
function name(alphabet::Alphabet{A,I}) where {A,I}
    return :custom
end

Base.getindex(alphabet::Alphabet, i::Integer) = alphabet.index_to_char[i]
Base.getindex(alphabet::Alphabet{A,I}, c::A) where {A,I} = alphabet.char_to_index[c]

function Base.show(io::IO, alphabet::Alphabet{A,I}) where {A,I}
    print(io, "Alphabet{$A,$I}: $(alphabet.characters)")
end
function Base.show(io::IO, x::MIME"text/plain", alphabet::Alphabet{A,I}) where {A,I}
    println(io, "$(name(alphabet)) Alphabet{$A,$I} with mapping $(alphabet.characters)")
end

"""
    symbols(alphabet)

Return the vector of symbols/characters used by `alphabet`.
"""
symbols(alphabet::Alphabet) = alphabet.characters

