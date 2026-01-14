"""
    read_fasta(fastafile::AbstractString; alphabet = :auto, kwargs...)
    read_fasta(
        fastafile::AbstractString, alphabet;
        weights=false, theta=0.2, verbose=false, safe=false
    )

If not provided, the alphabet is determined from the first sequences of the fasta file,
using the `BioSequenceMappings.auto_alphabet_from_fasta` procedure.
If unknown characters are found later on, an error will be raised.
This can be avoided by using the `safe` kwarg: sequences will unknown characters are then
excluded from the alignment and a warning is raised.
"""
function read_fasta(fastafile::AbstractString; alphabet = :auto, kwargs...)
    return read_fasta(fastafile, alphabet; kwargs...)
end

function read_fasta(fastafile::AbstractString, alphabet::Symbol; kwargs...)
    alphabet = if alphabet == :auto
        alphabet = auto_alphabet_from_fasta(fastafile)
    else
        Alphabet(alphabet)
    end
    return read_fasta(fastafile, alphabet; kwargs...)
end

function read_fasta(
    fastafile::AbstractString, alphabet::Alphabet;
    weights=false, theta=0.2, verbose=false, safe=false,
)
    verbose && @info "Reading sequences from $fastafile using Alphabet $alphabet"
    data = map(FASTAReader(open(fastafile))) do rec
        description(rec), sequence(rec)
    end
    @check allequal(Iterators.map(x -> length(x[2]), data)) Error("""
        All sequences must have the same length
    """)
    @check length(data) > 0 Error("Alignment $fastafile empty")

    aln = if safe
        M = length(data)
        sequences = []
        names = []
        errors = []
        n_invalid = 0
        for (i, (header, seq)) in enumerate(data)
            try
                seq = alphabet(seq)
                push!(sequences, seq)
                push!(names, header)
            catch err
                push!(errors, err)
                n_invalid += 1
                verbose && @warn "Could not read sequence $seq - Raised error $err"
            end
        end
        if n_invalid > 0
            @warn "Could not read $n_invalid sequences out of $M"
            # @warn """Could not read $n_invalid sequences out of $M
            # Error messages were $(errors[1:min(n_invalid, 5)]) (truncated)
            # """
        end

        Alignment(hcat(sequences...), alphabet; names)
    else
        aln = Alignment(
            mapreduce(x -> alphabet(x[2]), hcat, data), alphabet;
            names = map(first, data)
        )
    end
    L, M = size(aln)
    verbose && @info "Found $M sequence of length $L"

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

"""
    write(io, A::Alignment)

Write `A` in fasta format.
"""
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



