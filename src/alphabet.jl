#=
A mapping is a `Dict{Char, Integer}`
For example, `Dict('A' => 1, 'C' => 2, ...)
Below is a utility function for mappings
=#
"""
    compute_mapping(s::AbstractString)

Return a `Dict{Int, Char}`: `Dict(i => c for (i,c) in enumerate(s))`.
"""
function compute_mapping(s::AbstractString, ::Type{T}=Int) where T <: Integer
    @assert allunique(s) "Cannot compute mapping from string with duplicates - got $s"
    return Dict{Char, T}(c => i for (i,c) in enumerate(s))
end
function reverse_mapping(mapping::Dict{T,U}) where {T, U}
    return Dict{U,T}(y => x for (x, y) in mapping)
end

#######################################################################################
####################################### Alphabet ######################################
#######################################################################################


"""
    struct Alphabet{T}
        string::String
        char_to_index::Dict{Char, T}
        index_to_char::Dict{T, Char}
    end

`Alphabet` can be constructed
- from a `String` and an optional type: `Alphabet(string::String[, T<:Integer])`
- from a mapping `Dict{Char, T}` where `T<:Integer`: `Alphabet(mapping)`
- from a `Symbol`, using default alphabets
- from an integer, using default alphabets (see `default_alphabets`).
"""
@kwdef struct Alphabet{T<:Integer}
    string::String # Alphabet string
    char_to_index::Dict{Char, T} = compute_mapping(string)
    index_to_char::Dict{T, Char} = reverse_mapping(char_to_index)
    default_char::Union{Nothing, Char} = nothing
    default_index::Union{Nothing, T} = _default_index_from_char(default_char, string)

    # Constructor with checks
    function Alphabet{T}(
        string, char_to_index, index_to_char, default_char, default_index
    ) where T
        @assert isconcretetype(T) "Use concrete type to build `Alphabet` - got $T"
        @assert length(string) == length(char_to_index) == length(index_to_char) """\
            Inconsistent lengths between string, char_to_index and index_to_char:\
            $(length(string)) vs $(length(char_to_index)) vs $(length(index_to_char))
            """

        # are string and char_to_index consistent
        @assert all(c -> haskey(char_to_index, c), string) """\
            Incomplete char_to_index: some symbols in $string not in $char_to_index
            """
        @assert all(x -> char_to_index[x[2]] == x[1], enumerate(string)) """\
            Inconsistency between string and mapping: $string -- $mapping
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
            """
        @assert isnothing(default_char) || occursin(default_char, string) """\
            Default char $(default_char) not found in alphabet string $(string)
            Maybe use `default_char = c` with `c` in your alphabet?
            """
        @assert isnothing(default_index) || haskey(index_to_char, default_index) """\
            Got `default_index=$(default_index)` but alphabet string of length $(length(string))
            """
        @assert (isnothing(default_index) && isnothing(default_char)) ||
            index_to_char[default_index] == default_char """\
            Defaults are inconsistent with mapping.\
            `default_index=$(default_index)` maps to $(char_to_index[default_index]),\
             but `default_char=$(default_char)`
            """

        # Checking for duplicates
        @assert allunique(string) """Cannot form alphabet from string with duplicates\
            - got $s
            """
        return new{T}(string, char_to_index, index_to_char, default_char, default_index)
    end
end

_default_index_from_char(c::AbstractChar, S::AbstractString) = findfirst(==(c), S)
_default_index_from_char(::Nothing, ::AbstractString) = nothing

################################################
############## Equality, copying ###############
################################################

function Base.:(==)(A::Alphabet{T}, B::Alphabet{U}) where {T,U}
    return T == U &&
    A.string == B.string &&
    A.default_char == B.default_char
    A.default_index == B.default_index
    A.index_to_char == B.index_to_char &&
    A.char_to_index == B.char_to_index
end

function Base.hash(alphabet::Alphabet{T}, h::UInt) where T
    h = hash(Alphabet{T}, h)
    h = hash(alphabet.string, h)
    h = hash(alphabet.index_to_char, h)
    h = hash(alphabet.char_to_index, h)
    h = hash(alphabet.default_char, h)
    h = hash(alphabet.default_index, h)
    return h
end

function Base.copy(A::Alphabet{T}) where T
    return Alphabet{T}(;
        string = A.string, # strings are immutable so I don't need to copy
        char_to_index = copy(A.char_to_index),
        index_to_char = copy(A.index_to_char),
        default_char = A.default_char,
        default_index = A.default_index,
    )
end
Base.convert(::Type{T}, A::Alphabet{T}) where T = A
Base.convert(::Type{Alphabet{T}}, A::Alphabet{T}) where T = A
function Base.convert(::Type{T}, A::Alphabet) where T <: Integer
    return Alphabet{T}(;
        string = A.string,
        char_to_index = convert(Dict{Char, T}, A.char_to_index),
        index_to_char = convert(Dict{T, Char}, A.index_to_char),
        default_char = A.default_char,
        default_index = isnothing(A.default_index) ? nothing : convert(T, A.default_index),
    )
end
function Base.convert(::Type{Alphabet{T}}, A::Alphabet) where T <: Integer
    return convert(T, A)
end

################################################
################# Constructors #################
################################################


function Alphabet(S::AbstractString, ::Type{T}=Int; kwargs...) where T <: Integer
    return Alphabet{T}(; string = S, char_to_index = compute_mapping(S, T), kwargs...)
end
Alphabet(A::Alphabet) = A
Alphabet(A::Alphabet{T}, ::Type{U}) where {T,U} = convert(U, A)
function Alphabet(S::AbstractString, mapping::Dict{<:AbstractChar, T}; kwargs...) where T<:Integer
    return Alphabet{T}(; string = S, char_to_index = mapping, kwargs...)
end
function Alphabet(mapping::AbstractDict{Char, T}; kwargs...) where T
    str = Vector{Char}(undef, length(mapping))
    for (c, i) in mapping
        str[i] = c
    end
    return Alphabet{T}(; string = prod(str), char_to_index = mapping, kwargs...)
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
        Alphabet(_DEFAULT_AA_ALPHABET_STRING, T; default_char = '-')
    elseif name in nt_alphabet_names
        Alphabet(_DEFAULT_NT_ALPHABET_STRING, T; default_char = '-')
    elseif name in binary_alphabet_names
        Alphabet(_DEFAULT_BINARY_ALPHABET_STRING, T)
    else
        names = vcat(aa_alphabet_names..., nt_alphabet_names..., binary_alphabet_names...)
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

    return if q == 21
        Alphabet(:aa, T)
    elseif q == 5
        Alphabet(:dna, T)
    elseif q == 4
        Alphabet(_DEFAULT_NT_ALPHABET_STRING_NOGAP, T)
    elseif q == 2
        Alphabet(:binary, T)
    elseif q < 21
        Alphabet(_DEFAULT_AA_ALPHABET_STRING[1:q], T)
    else
        error("No defined default alphabet for q = $q > 21 - provide your own or use `nothing`")
    end
end

################################################
############ Transforming sequences ############
################################################


function (alphabet::Alphabet{T})(c::Char) where T
    i = get(alphabet.char_to_index, c, alphabet.default_index)
    if isnothing(i)
        error("Char $c not in alphabet, and no defaults set.")
    end
    return i
end
function (alphabet::Alphabet{T})(S::AbstractString) where T
    return map(x -> alphabet(x), collect(S))
end

function (alphabet::Alphabet)(x::Integer)
    c = get(alphabet.index_to_char, x, alphabet.default_char)
    if isnothing(c)
        error("Index $x does not represent any characters, and not defaults set.")
    end
    return c
end
function (alphabet::Alphabet)(X::AbstractVector{<:Integer})
    return string(map(a -> alphabet(a), X)...)
end

(alphabet::Alphabet)(::Missing) = missing

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

Base.length(alphabet::Alphabet) = length(alphabet.string)

Base.in(i::Integer, alphabet::Alphabet) = haskey(alphabet.index_to_char, i)
Base.in(c::Char, alphabet::Alphabet) = haskey(alphabet.char_to_index, c)

function name(A::Alphabet)
    return if A.string == _DEFAULT_AA_ALPHABET_STRING
        :aa
    elseif A.string == _DEFAULT_NT_ALPHABET_STRING
        :dna
    elseif A.string == _DEFAULT_BINARY_ALPHABET_STRING
        :binary
    else
        :custom
    end
end

Base.getindex(alphabet::Alphabet, i::Integer) = alphabet.index_to_char[i]
Base.getindex(alphabet::Alphabet, c::Char) = alphabet.char_to_index[c]

function Base.show(io::IO, A::Alphabet{T}) where T
    print(io, "Alphabet{$T}: \"$(A.string)\"")
end
function Base.show(io::IO, x::MIME"text/plain", A::Alphabet{T}) where T
    println(io, "$(name(A)) Alphabet{$T} with mapping \"$(A.string)\"")
end

