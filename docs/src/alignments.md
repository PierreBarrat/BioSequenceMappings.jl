```@meta 
DocTestSetup = quote
	using BioSequenceMappings
end
```

# The `Alignment` type

An [`Alignment`](@ref) essentially contains a set of aligned biological sequences mapped to integers. 
It is a subtype of the more general `AbstractAlignment`. 
The precise type is of the form `Alignment{A,T}` where `A` is the type of biological symbols (see [The `Alphabet` type](@ref)) and `T` the type of integers. 
Its fields are simple and self-explanatory:
- `data::Matrix{T}`: the mapped sequences. Each column of the matrix is one sequence. 
- `alphabet::Alphabet{A,T}`: the alphabet defining the mapping from symbols to integers.
- `weights::Vector{Float64}`: weights associated with each sequence. DCA-like methods assign phylogenetic weights to sequences when infering models. 
- `names::Vector{String}`: names of the sequences, corresponding to the description string if the alignment was built from a FASTA file. 

The dimensions of an alignment can be obtained using `size`: length $\cross$ number of sequences

## Reading & writing

For now, the only file format that this package interacts with is FASTA. 
Reading an alignment is simple: 
```@repl align1
using BioSequenceMappings # hide
A = read_fasta("../../example/PF00014.fasta")
size(A) # length and number of sequences
```

When reading from FASTA, the choice of the alphabet is made by reading the first five sequences, and comparing the observed characters with the list of default alphabets (see [The `Alphabet` type](@ref)). 
If they fit one of the defaults, it will be used. 
Otherwise, an alphabet will be created *ad hoc*: 
```@repl strange_characters
using BioSequenceMappings # hide
A = read_fasta("../../example/strange_characters.fasta"); # warning produced because no default alphabet was found
A.alphabet |> symbols |> prod
```

Writing to a FASTA file is just as easy: 
```@repl strange_characters
write("new_fasta_file.fasta", A) # or...
open("new_fasta_file.fasta", "w") do io
	write(io, A)
end
```

## Accessing & iterating

Sequences can be accessed by indexing. 
Indexing using a range will return a view in the underlying `data` matrix. 
```@repl align1
A[1] # the first sequence of the alignment
size(A[1:5]) # 5 sequences of length L
```

!!! tip
	When indexing or when iterating, the return value is a *reference* to the sequence and not a copy. 
	For this reason, if you modify the sequence, the alignment itself will be modified
	```@repl 
	using BioSequenceMappings # hide
	A = read_fasta("../../example/PF00014.fasta") # hide
	s = A[1]'
	s[1] = 21
	A[1]' # the first element has been modified
	```

The fact that calls like `A[1:5]` return a matrix-like object can be inconvenient. 
To iterate over sequences as vectors of integers, one can use the [`eachsequence`](@ref) function. 
It takes the alignment as a first argument and optionally indices, and return an iterator. 
```@repl align1
for s in eachsequence(A)
	# s: vector of integers
end
map(length, eachsequence(A, 1:5:16)) # length of sequences 1, 6, 11. 16
collect(eachsequence(A; skip=25)) # collect as list of vectors, taking one every 25 sequences
```

If the name of the sequence is needed, the iterator [`named_sequences`](@ref) can be used instead, taking the same arguments and iterating over tuples of the form `(name, sequence)`. 

## Finding sequences

The package provides [`find_sequence`](@ref) to find sequences by name in an alignment, and [`match_sequences`](@ref) to match all sequences with a particular name pattern. 

```jldoctest align2; setup = :(example_dir = joinpath(dirname(pathof(BioSequenceMappings)), "../example"))
julia> A = read_fasta(joinpath(example_dir, "toy_fasta_dna.fasta"));

julia> n = A.names[1] # name of the first sequence in the alignment
"sequence_1"

julia> find_sequence(n, A) # index and sequence
(1, [5, 3, 4, 2, 4, 1, 5, 1, 5, 5])

julia> indices, sequences = match_sequences(r"sequence_[0-9]", A); # using a regex 

julia> indices
5-element Vector{Int64}:
 1
 2
 3
 4
 5
```

!!! note
	`find_sequence` searches for an exact match. 
	It is based on `findfirst`, and will thus return the first match it finds. 
	If nothing is found, then the result is `nothing`. 
	On the other hand, `match_sequences` is based on `occursin` and `findall`. 
	The returned sequences are references to original objects, as when indexing. 


## Creating subalignments

The call `subsample(A, indices)` will create an `AbstractAlignment` of the same type as `A` by taking only the sequences at `indices`. 
It uses the same alphabet as `A` and copies over the names of the sequences. 
Note that [`subsample`](@ref) *copies* the underlying data, creating a completely independent object.
```@repl align3
using BioSequenceMappings # hide
A = read_fasta("../../example/toy_fasta_dna.fasta")
A[1][1] # first character of the first sequence
B = subsample(A, 1:2:5)
B[1][1] = 5
A[1][1] # A remains unchanged
```

With [`subsample_random`](@ref), it is also possible to create a random subalignment by picking sequences from the original one. 
For now, this is only possible without replacement, *i.e.* the same sequence cannot be picked twice. 
To just pick one sequence at random without creating a new alignment object, just call `rand`. 
```@repl align3
subsample_random(A, 3) # new alignment using three random sequences from A
subsample_random(A, 12) # sampling without replacement: this will error since size(A, 1) < 12
rand(A) # one random sequence from A (returns a view)
```

## OneHotAlignment

TBA