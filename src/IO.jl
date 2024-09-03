"""
    read_fasta(fastafile::AbstractString; alphabet = :auto, kwargs...)
    read_fasta(
        fastafile::AbstractString, alphabet;
        weights = false, theta = 0.2, verbose = false,
    )
"""
function read_fasta(fastafile::AbstractString; alphabet = :auto, kwargs...)
    return read_fasta(fastafile, alphabet; kwargs...)
end

function read_fasta(fastafile::AbstractString, alphabet::Symbol; kwargs...)
    if alphabet == :auto
        alphabet = auto_alphabet_from_fasta(fastafile)
    else
        error("`alphabet` must be `:auto` or an `Alphabet` object")
    end
    return read_fasta(fastafile, alphabet; kwargs...)
end

function read_fasta(
    fastafile::AbstractString, alphabet::Alphabet;
    weights = false, theta = 0.2, verbose = false,
)
    verbose && @info "Reading sequences from $fastafile using Alphabet $alphabet"
    data = map(FASTAReader(open(fastafile))) do rec
        description(rec), sequence(rec)
    end
    @assert allequal(Iterators.map(x -> length(x[2]), data)) """
        All sequences must have the same length
    """

    aln = Alignment(
        mapreduce(x -> alphabet(x[2]), hcat, data), alphabet;
        names = map(first, data)
    )
    L, M = size(aln)
    verbose && @info "Found $M sequence of length $M"

    if weights
        compute_weights!(aln, theta)
    end

    return aln
end

function auto_alphabet_from_fasta(fastafile::AbstractString; n = 5)
    # first n sequences (at most)
    sequences = map(sequence, first(FASTAReader(open(fastafile)), n))
    return auto_alphabet_from_sequences(sequences)
end
function auto_alphabet_from_sequences(sequences::AbstractVector{<:AbstractString})
    characters = sort(unique(prod(sequences)))
    return if all(in(_DEFAULT_NT_ALPHABET_STRING), characters)
        Alphabet(:dna)
    elseif all(in(_DEFAULT_AA_ALPHABET_STRING), characters)
        Alphabet(:aa)
    elseif all(in(_DEFAULT_BINARY_ALPHABET_STRING), characters)
        Alphabet(:binary)
    else
        alphabet = Alphabet(prod(characters))
        @warn "Could not find a default alphabet for characters $characters\
        \n Using $alphabet"
        alphabet
    end
end

function Base.write(file::AbstractString, X::Alignment)
    return open(file, "w") do io
        write(io, X)
    end
end
function Base.write(io::IO, X::Alignment)
    return FASTAWriter(io) do fw
        for (i, seq) in enumerate(X)
            header = isempty(X.names[i]) ? "$i" : X.names[i]
            rec = FASTARecord(header, to_string(seq, X.alphabet))
            write(fw, rec)
        end
    end
end



