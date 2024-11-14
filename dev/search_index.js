var documenterSearchIndex = {"docs":
[{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"using BioSequenceMappings","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"DocTestSetup = quote\n\tusing BioSequenceMappings\n\tA = Alphabet(['A', 'C', 'G', 'T', '-'])\nend","category":"page"},{"location":"alphabets/#The-Alphabet-type","page":"Alphabets","title":"The Alphabet type","text":"","category":"section"},{"location":"alphabets/#Basics","page":"Alphabets","title":"Basics","text":"","category":"section"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"An Alphabet contains information necessary to map biological symbols to integers and inversely.  The full type is Alphabet{A,I}, where A is the type of biological symbols (typically Char) and I is a subtype of Integer. ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"The simplest way to create an alphabet is from a list of symbols: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"A = Alphabet(['A', 'C', 'G', 'T', '-'])","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"The created alphabet A associates Char correspondings to nucleotides to Int according to the index at which they appear in the input vector: 'A' => 1, 'C' => 2, etc... Note that we could have created the same alphabet from a string (since it's based on Chars) or from a dictionary: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"julia> B = Alphabet(\"ACGT-\");\n\njulia> C = Alphabet(Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, '-'=>5));\n\njulia> A == B\ntrue\n\njulia> A == C\ntrue","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"The alphabet is used to map symbols to integers and inversely.  This is done by calling the object directly, as a function: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"julia> A('A') # mapping Char to Int\n1\n\njulia> A(1) # mapping Int to Char\n'A': ASCII/Unicode U+0041 (category Lu: Letter, uppercase)\n\njulia> A(\"A-GT\") # mapping a string to a vector \n4-element Vector{Int64}:\n 1\n 5\n 3\n 4\n\njulia> A([1,2,3]) # mapping a vector to a string\n\"ACG\"","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"If needed, the mapping is accessible using the symbolsfunction: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"julia> symbols(A)\n5-element Vector{Char}:\n 'A': ASCII/Unicode U+0041 (category Lu: Letter, uppercase)\n 'C': ASCII/Unicode U+0043 (category Lu: Letter, uppercase)\n 'G': ASCII/Unicode U+0047 (category Lu: Letter, uppercase)\n 'T': ASCII/Unicode U+0054 (category Lu: Letter, uppercase)\n '-': ASCII/Unicode U+002D (category Pd: Punctuation, dash)\n\njulia> A |> symbols |> prod # as a string\n\"ACGT-\"","category":"page"},{"location":"alphabets/#Default-alphabets","page":"Alphabets","title":"Default alphabets","text":"","category":"section"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"The package comes with three default alphabets: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"an amino-acid alphabet Alphabet(:aa) using the mapping \"-ACDEFGHIKLMNPQRSTVWY\";\na nucleotide alphabet Alphabet(:dna) using the mapping \"-ACGT\";\na \"binary\" alphabet Alphabet(:binary), which I found useful for simulations, with the mapping: \"01\". ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"They can be accessed by calling Alphabet(name) where name is a symbol corresponding to any of the default alphabets.  The symbolic names can be easily be found:","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"julia> BioSequenceMappings.aa_alphabet_names # also works with nt and binary alphabets\n(:aa, :AA, :aminoacids, :amino_acids)\n\njulia> Alphabet(:aa) == Alphabet(\"-ACDEFGHIKLMNPQRSTVWY\")\ntrue\n\njulia> Alphabet(:amino_acids)([1,2,3])\n\"-AC\"","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"Each default alphabet is also associated to a specific cardinality of biological symbols through the function default_alphabet.  This means that an integer vector with elements ranging from 1 to q will be associated to the following alphabets: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"julia> default_alphabet(2) == Alphabet(:binary) # q == 2\ntrue\n\njulia> default_alphabet(5) == Alphabet(:nt) # q == 5\ntrue\n\njulia> default_alphabet(21) == Alphabet(:aa) # 5 < q <= 21\ntrue\n\njulia> default_alphabet(15) == Alphabet(:aa) # 5 < q <= 21","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"This association is useful to create Alignment objects from a matrix of integers without having to specify the alphabet manually. ","category":"page"},{"location":"alphabets/#Default-characters","page":"Alphabets","title":"Default characters","text":"","category":"section"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"When reading biological sequences, it can be convenient to associate all unexpected characters to a default symbol, for instance the gap.  This can be achieved by providing the default_char keyword argument when constructing the alphabet: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"A_default = Alphabet(\"ACGT-\"; default_char = '-')\nA_default(\"ABCDEF\") # 'unknown' chars are mapped to '-', in turn mapped to 5\nA(\"ABCDEF\") # if no defaults are provided, fails","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"This also works the other way around: integers that are not in the range of the alphabet are mapped to the default symbol: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"A_default(1:10) # indices larger than 5 are mapped to the gap\nA(1:10) # if no defaults are provided, fails","category":"page"},{"location":"alphabets/#Using-specific-integer-types","page":"Alphabets","title":"Using specific integer types","text":"","category":"section"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"When created as above, the alphabet will default to using Int as the integer type.  If dealing with large amounts of data, it can be beneficial to use a more sober type.  This is done by providing an extra argument of the desired type when constructing the alphabet: ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"B = Alphabet(\"ACGT-\", UInt8)\nB == A\nB(\"A-\")","category":"page"},{"location":"alphabets/#Translating-between-alphabets","page":"Alphabets","title":"Translating between alphabets","text":"","category":"section"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"It often happens to me that I have an integer vector X representing a sequence, but with a mapping different from the one I am used to.  The translate function lets me convert it to another integer vector with the right mapping. ","category":"page"},{"location":"alphabets/","page":"Alphabets","title":"Alphabets","text":"julia> strange_alphabet = Alphabet(\"TCGA-\"); # the default is \"-ACGT\"\n\njulia> nt_alphabet = Alphabet(:dna); # the default for nucleotides\n\njulia> X = Int[2, 2, 5, 4, 5]; # representing the sequence \"CC-A-\" according to the above\n\njulia> strange_alphabet(X)\n\"CC-A-\"\n\njulia> nt_alphabet(X) # this is obviously wrong - nt_alphabet uses \"-ACGT\"\n\"AATGT\"\n\njulia> Y = translate(X, strange_alphabet, nt_alphabet)\n5-element Vector{Int64}:\n 3\n 3\n 1\n 2\n 1\n\njulia> nt_alphabet(Y) # now this works\n\"CC-A-\"","category":"page"},{"location":"utilities/#Utilities","page":"Utilities","title":"Utilities","text":"","category":"section"},{"location":"utilities/#Hamming-distance","page":"Utilities","title":"Hamming distance","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"hamming(x, y) returns the normalized Hamming distance between two integer vectors. ","category":"page"},{"location":"utilities/#Alignment-statistics","page":"Utilities","title":"Alignment statistics","text":"","category":"section"},{"location":"utilities/#Pairwise-Hamming-distance","page":"Utilities","title":"Pairwise Hamming distance","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"pairwise_hamming(A) returns a vector containing all Hamming distances between pairs of sequences in A.  For large alignments, it is often practical to consider only a subset of Hamming distances: the step=n keyword can be used to only consider every nth sequence.  The function is also adapted to two alignments: [pairwise_hamming(A,B)] will consider all pairs of sequences with one member in A and the other one in B. ","category":"page"},{"location":"utilities/#Statistics","page":"Utilities","title":"Statistics","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"the profile of the alignment: site_specific_frequencies(A) returns a q x L matrix where q is the size of the alphabet and L the length of sequences, with element (a, i) being the fraction of sequences in which character a was found at position i. \nthe pairwise frequencies: pairwise_frequencies(A)\nthe pairwise correlations: pairwise_correlations(A)","category":"page"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"All these functions accept a vector of weights as a second argument, and will by default use A.weights if it is not provided. ","category":"page"},{"location":"utilities/#Weights","page":"Utilities","title":"Weights","text":"","category":"section"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"In DCA-like algorithms, it is customary to weight sequences by their degree of phylogenetic relations.  Typically, the weight of a sequence is inversely proportional to the number of other sequences with a Hamming distance smaller than some threshold.  For an alignment X, computing the weights is done using compute_weights: ","category":"page"},{"location":"utilities/","page":"Utilities","title":"Utilities","text":"using BioSequenceMappings # hide\nA = read_fasta(\"../../example/toy_fasta_dna.fasta\")\ncompute_weights(A) # compute and return the weight vector, as well as the effective number of sequences\ncompute_weights!(A); # same, but sets the weights field in A\nA.weights","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"DocTestSetup = quote\n\tusing BioSequenceMappings\nend","category":"page"},{"location":"alignments/#The-Alignment-type","page":"Alignments","title":"The Alignment type","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"An Alignment essentially contains a set of aligned biological sequences mapped to integers.  It is a subtype of the more general AbstractAlignment.  The precise type is of the form Alignment{A,T} where A is the type of biological symbols (see The Alphabet type) and T the type of integers.  Its fields are simple and self-explanatory:","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"data::Matrix{T}: the mapped sequences. Each column of the matrix is one sequence. \nalphabet::Alphabet{A,T}: the alphabet defining the mapping from symbols to integers.\nweights::Vector{Float64}: weights associated with each sequence. DCA-like methods assign phylogenetic weights to sequences when infering models. \nnames::Vector{String}: names of the sequences, corresponding to the description string if the alignment was built from a FASTA file. ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"The dimensions of an alignment can be obtained using size: length cross number of sequences","category":"page"},{"location":"alignments/#Reading-and-writing","page":"Alignments","title":"Reading & writing","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"For now, the only file format that this package interacts with is FASTA.  Reading an alignment is simple: ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"using BioSequenceMappings # hide\nA = read_fasta(\"../../example/PF00014.fasta\")\nsize(A) # length and number of sequences","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"When reading from FASTA, the choice of the alphabet is made by reading the first five sequences, and comparing the observed characters with the list of default alphabets (see The Alphabet type).  If they fit one of the defaults, it will be used.  Otherwise, an alphabet will be created ad hoc: ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"using BioSequenceMappings # hide\nA = read_fasta(\"../../example/strange_characters.fasta\"); # warning produced because no default alphabet was found\nA.alphabet |> symbols |> prod","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"Writing to a FASTA file is just as easy: ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"write(\"new_fasta_file.fasta\", A) # or...\nopen(\"new_fasta_file.fasta\", \"w\") do io\n\twrite(io, A)\nend","category":"page"},{"location":"alignments/#Accessing-and-iterating","page":"Alignments","title":"Accessing & iterating","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"Sequences can be accessed by indexing.  Indexing using a range will return a view in the underlying data matrix. ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"A[1] # the first sequence of the alignment\nsize(A[1:5]) # 5 sequences of length L","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"tip: Tip\nWhen indexing or when iterating, the return value is a reference to the sequence and not a copy.  For this reason, if you modify the sequence, the alignment itself will be modifiedusing BioSequenceMappings # hide\nA = read_fasta(\"../../example/PF00014.fasta\") # hide\ns = A[1]'\ns[1] = 21\nA[1]' # the first element has been modified","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"The fact that calls like A[1:5] return a matrix-like object can be inconvenient.  To iterate over sequences as vectors of integers, one can use the eachsequence function.  It takes the alignment as a first argument and optionally indices, and return an iterator. ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"for s in eachsequence(A)\n\t# s: vector of integers\nend\nmap(length, eachsequence(A, 1:5:16)) # length of sequences 1, 6, 11. 16\ncollect(eachsequence(A; skip=25)) # collect as list of vectors, taking one every 25 sequences","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"If the name of the sequence is needed, the iterator named_sequences can be used instead, taking the same arguments and iterating over tuples of the form (name, sequence). ","category":"page"},{"location":"alignments/#Finding-sequences","page":"Alignments","title":"Finding sequences","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"The package provides find_sequence to find sequences by name in an alignment, and match_sequences to match all sequences with a particular name pattern. ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"julia> A = read_fasta(joinpath(example_dir, \"toy_fasta_dna.fasta\"));\n\njulia> n = A.names[1] # name of the first sequence in the alignment\n\"sequence_1\"\n\njulia> find_sequence(n, A) # index and sequence\n(1, [5, 3, 4, 2, 4, 1, 5, 1, 5, 5])\n\njulia> indices, sequences = match_sequences(r\"sequence_[0-9]\", A); # using a regex \n\njulia> indices\n5-element Vector{Int64}:\n 1\n 2\n 3\n 4\n 5","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"note: Note\nfind_sequence searches for an exact match.  It is based on findfirst, and will thus return the first match it finds.  If nothing is found, then the result is nothing.  On the other hand, match_sequences is based on occursin and findall.  The returned sequences are references to original objects, as when indexing. ","category":"page"},{"location":"alignments/#Creating-subalignments","page":"Alignments","title":"Creating subalignments","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"The call subsample(A, indices) will create an AbstractAlignment of the same type as A by taking only the sequences at indices.  It uses the same alphabet as A and copies over the names of the sequences.  Note that subsample copies the underlying data, creating a completely independent object.","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"using BioSequenceMappings # hide\nA = read_fasta(\"../../example/toy_fasta_dna.fasta\")\nA[1][1] # first character of the first sequence\nB = subsample(A, 1:2:5)\nB[1][1] = 5\nA[1][1] # A remains unchanged","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"With subsample_random, it is also possible to create a random subalignment by picking sequences from the original one.  For now, this is only possible without replacement, i.e. the same sequence cannot be picked twice.  To just pick one sequence at random without creating a new alignment object, just call rand. ","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"subsample_random(A, 3) # new alignment using three random sequences from A\nsubsample_random(A, 12) # sampling without replacement: this will error since size(A, 1) < 12\nrand(A) # one random sequence from A (returns a view)","category":"page"},{"location":"alignments/#Misc.","page":"Alignments","title":"Misc.","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"Sorting / concatenating alignments.","category":"page"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"julia> n = A.names # Sorting an alignment\n5-element Vector{String}:\n \"sequence_1\"\n \"sequence_2\"\n \"sequence_3\"\n \"sequence_4\"\n \"sequence_5\"\n\njulia> B = sort(A; with_name=true, rev=true); # could also use `by=...`\n\njulia> B.names\n5-element Vector{String}:\n \"sequence_5\"\n \"sequence_4\"\n \"sequence_3\"\n \"sequence_2\"\n \"sequence_1\"\n\njulia> C = cat(A, B); # concatenate two alignments\n\njulia> size(C) # 10 sequences of length 10\n(10, 10)","category":"page"},{"location":"alignments/#OneHotAlignment","page":"Alignments","title":"OneHotAlignment","text":"","category":"section"},{"location":"alignments/","page":"Alignments","title":"Alignments","text":"TBA","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"CurrentModule = BioSequenceMappings","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Documentation for BioSequenceMappings.","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [BioSequenceMappings]\nPrivate = false","category":"page"},{"location":"reference/#BioSequenceMappings.Alignment","page":"Reference","title":"BioSequenceMappings.Alignment","text":"mutable struct Alignment{A,T} where {A, T<:Integer}\n\n    data::Matrix{T}\n    alphabet::Union{Nothing, Alphabet{A,T}}\n    weights::Vector{Float64} = ones(size(dat,1))/size(dat,1) # phylogenetic weights of sequences\n    names::Vector{String} = fill(\"\", size(dat, 1))\n\nBiological sequences as vectors of type T<:Integer. data stores sequences in columns: size(dat) returns a tuple (L, M) with L the length and M the number of sequences. When displayed, shows data as an MxL matrix to match with traditional alignments.\n\nalphabet{A,T} represents the mapping between integers in data and biological symbols of type A (nucleotides, amino acids...). If nothing, the alignment cannot be mapped to biological sequences.\n\nweights represent phylogenetic weights, and are initialized to 1/M. They must sum to 1. names are the label of sequences, and are expected to be in the same order as the columns of data. They do not have to be unique, and can be ignored\n\nImportant: When built from a matrix, assumes that the sequences are stored in columns.\n\nMethods\n\ngetindex(X::Alignment, i) returns a matrix/vector X.data[:, i].\nfor s in X::Alignment iterates over sequences.\neachsequence(X::Alignment) returns an iterator over sequences (Vector{Int}).\neachsequence_weighted(X::Alignment) returns an iterator over sequences and weights as tuples\nsubaln(X::Alignment, idx) constructs the subaln defined by index idx.\n\n\n\n\n\n","category":"type"},{"location":"reference/#BioSequenceMappings.Alignment-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:Integer","page":"Reference","title":"BioSequenceMappings.Alignment","text":"Alignment(data::AbstractMatrix{T}; alphabet = :auto, kwargs...)\n\nKeyword argument alphabet can be :auto, :none/nothing, or an input to the constructor Alphabet. Other keyword arguments are passed to the default constructor of Alignment.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.Alignment-Union{Tuple{T}, Tuple{A}, Tuple{AbstractMatrix{T}, Alphabet{A, T}}} where {A, T}","page":"Reference","title":"BioSequenceMappings.Alignment","text":"Alignment(data::AbstractMatrix, alphabet; kwargs...)\n\ndata is a matrix of integers, with sequences stored in columns. alphabet can be either\n\nan Alphabet\nnothing: no conversion from integers to biological symbols.\nsomething to build an alphabet from (e.g. a symbol like :aa, a string, ...).   The constructor Alphabet will be called like so: Alphabet(alphabet).\n\nIf the types of alphabet and data mismatch, data is converted.\n\ndata can also have the following shape:\n\nvector of integer vectors, e.g. [[1,2], [3,4]]: each element is considered as a sequence\nvector of integers: single sequence alignment\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.Alphabet","page":"Reference","title":"BioSequenceMappings.Alphabet","text":"struct Alphabet{A,I}\n    characters::Vector{A}\n    char_to_index::Dict{A, I}\n    index_to_char::Dict{I, A}\n    default_char = nothing\n    default_index\nend\n\nStructure allowing the mapping from biological symbols of type A to integers of type I.     The typical use case would be Alphabet{Char, Int}. Alphabet can be constructed\n\nfrom a Vector of symbols and an optional type I, e.g. Alphabet(['A','C','G','T'], UInt8)::Alphabet{Char, UInt8}\nfrom a String and an optional type, e.g. Alphabet(\"ACGT\")\nfrom a mapping Dict{A, I} where I<:Integer: Alphabet(Dict('A'=>1, 'C'=>2))\nfrom a Symbol, using default alphabets, e.g. Alphabet(:nt)\nfrom an integer, using default alphabets (see ?default_alphabets).\n\n\n\n\n\n","category":"type"},{"location":"reference/#BioSequenceMappings.compute_weights","page":"Reference","title":"BioSequenceMappings.compute_weights","text":"compute_weights(X::AbstractAlignment, θ = 0.2; normalize = true)\n\nCompute phylogenetic correction weights for sequences of X.     The weight sequence S is 1/N, where N is the number of sequences in X at     hamming distance less than H from S (including S itself).     The threshold H is floor(θ⋅L) where L is the sequence length.\n\nThe return value is a tuple (weights, Meff), where Meff is the sum of weights (pre-normalization). If normalize, weights are normalized to sum to one. .\n\n\n\n\n\n","category":"function"},{"location":"reference/#BioSequenceMappings.compute_weights!","page":"Reference","title":"BioSequenceMappings.compute_weights!","text":"compute_weights!(X, θ; kwargs...)\n\nCompute and set weights for X. See compute_weights.\n\n\n\n\n\n","category":"function"},{"location":"reference/#BioSequenceMappings.default_alphabet-Union{Tuple{Integer}, Tuple{T}, Tuple{Integer, Type{T}}} where T<:Integer","page":"Reference","title":"BioSequenceMappings.default_alphabet","text":"default_alphabet(q::Int, T::Type)\n\nif q==2, binary (0, 1)\nif 3 <= q <= 4, nucleotides without gaps\nif q==5, nucleotides\nif 5 < q <= 21, amino acids\nif q>21, fails\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.eachsequence-Tuple{AbstractAlignment, Any}","page":"Reference","title":"BioSequenceMappings.eachsequence","text":"eachsequence(X::AbstractAlignment[, indices]; skip)\n\nReturn an iterator over the sequences in X. If indices is specified, consider only sequences at the corresponding indices. Use the integer argument skip to return only one sequence every skip (~ 1:skip:end).\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.find_sequence-Tuple{AbstractString, AbstractAlignment}","page":"Reference","title":"BioSequenceMappings.find_sequence","text":"find_sequence(label::AbstractString, aln::AbstractAlignment)\n\nFind sequence with name label in aln, and return (index, sequence). Scales as the number of sequences. Return the first sequence that matches the label.\n\n!!! Return a view of the sequence.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.hamming-Tuple{AbstractVector{<:Integer}, AbstractVector{<:Integer}}","page":"Reference","title":"BioSequenceMappings.hamming","text":"hamming(x, y; normalize=true, positions=nothing)\n\nHamming distance between Vectors x and y. Only sites in vector positions will be considered.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.match_sequences-Tuple{Any, AbstractAlignment}","page":"Reference","title":"BioSequenceMappings.match_sequences","text":"match_sequences(pattern, aln::AbstractAlignment)\n\nFind sequences whose name matches label in aln, and return (indices, sequences). Sequences are returned as columns.\n\n!!! Return a view of the sequences.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.named_sequences-Tuple{AbstractAlignment, Any}","page":"Reference","title":"BioSequenceMappings.named_sequences","text":"named_sequences(X::AbstractAlignment; skip)\n\nReturn an iterator of the form (name, sequence) over X.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.pairwise_correlations","page":"Reference","title":"BioSequenceMappings.pairwise_correlations","text":"pairwise_correlations(X, w=X.weights; as_mat=false)\n\nCompute connected correlations: the difference between the pairwise frequencies and the product of the single site frequencies. See ?pairwise_frequencies for the shape of the output.\n\n\n\n\n\n","category":"function"},{"location":"reference/#BioSequenceMappings.pairwise_frequencies","page":"Reference","title":"BioSequenceMappings.pairwise_frequencies","text":"pairwise_frequencies(X::AbstractAlignment, w=X.weights; as_mat=false)\n\nReturn a q x q x L x L tensor. The (a, b, i, j) element is the fraction of sequences for which we see a at position i and b at position j.\n\nIf as_mat=true, will return a qL x qL matrix, with q x q blocks representing correlations between two specific columns.\n\n\n\n\n\n","category":"function"},{"location":"reference/#BioSequenceMappings.pairwise_hamming-Tuple{AbstractAlignment, AbstractAlignment}","page":"Reference","title":"BioSequenceMappings.pairwise_hamming","text":"pairwise_hamming(X, Y; step=1, step_left, step_right, as_vec=true, kwargs...)\npairwise_hamming(X; as_vec, step, kwargs...)\n\nReturn all hamming distances between sequences of X and Y. In the second form, consider pairs of sequences in X.\n\nOnly consider sequences every step. step_left and step_right can be used to skip sequence either in X or in Y. This is useful for large alignment, as the number of computations grows with the product     of the size of the alignments\n\nBy default, the return value is a vector organized like [H(1,2), H(1,3), ..., H(M-1, M)] with H standing for hamming distance and M for the number of sequences. If a matrix is prefered, use as_vec=false\n\nExtra keyword arguments are passed to hamming.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.read_fasta-Tuple{AbstractString}","page":"Reference","title":"BioSequenceMappings.read_fasta","text":"read_fasta(fastafile::AbstractString; alphabet = :auto, kwargs...)\nread_fasta(\n    fastafile::AbstractString, alphabet;\n    weights = false, theta = 0.2, verbose = false,\n)\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.site_specific_frequencies","page":"Reference","title":"BioSequenceMappings.site_specific_frequencies","text":"site_specific_frequencies(X::AbstractAlignment[, weights=X.weights]; as_vec=false)\n\nReturn the site specific frequencies of X. If as_vec, the result is a vector of length Lxq. Otherwise, it is a matrix of q rows and L columns (default).\n\n\n\n\n\n","category":"function"},{"location":"reference/#BioSequenceMappings.subsample-Tuple{AbstractAlignment, Int64}","page":"Reference","title":"BioSequenceMappings.subsample","text":"subsample(X::AbstractAlignment, indices)\n\nReturn an Alignment containing only the sequences of X at indices.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.subsample_random-Tuple{AbstractAlignment, Int64}","page":"Reference","title":"BioSequenceMappings.subsample_random","text":"subsample_random(X::AbstractAlignment, m::Int)\n\nReturn an Alignment with m sequences taking randomly from X. Sampling is done without replacement, meaning the m sequences are all at different positions in X.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.symbols-Tuple{Alphabet}","page":"Reference","title":"BioSequenceMappings.symbols","text":"symbols(alphabet)\n\nReturn the vector of symbols/characters used by alphabet.\n\n\n\n\n\n","category":"method"},{"location":"reference/#BioSequenceMappings.translate-Tuple{Integer, Alphabet, Alphabet}","page":"Reference","title":"BioSequenceMappings.translate","text":"translate(x, original_alphabet::Alphabet, new_alphabet::Alphabet)\n\nReturn the translation in new_alphabet of an integer or a vector of integers x that is expressed in original_alphabet.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Quickstart","title":"Quickstart","text":"CurrentModule = BioSequenceMappings","category":"page"},{"location":"#BioSequenceMappings","page":"Quickstart","title":"BioSequenceMappings","text":"","category":"section"},{"location":"","page":"Quickstart","title":"Quickstart","text":"The aim of the package is to facilitate the task of converting biological sequences (nucleotides, amino acids) to integers or to onehot representation.  It also provides simple function to manipulate alignments or to compute simple statistics. ","category":"page"},{"location":"","page":"Quickstart","title":"Quickstart","text":"Note: I do not frequently use the onehot format, so for now it is not documented / not well tested. I'll develop it more when I get the time. ","category":"page"},{"location":"#Installation","page":"Quickstart","title":"Installation","text":"","category":"section"},{"location":"","page":"Quickstart","title":"Quickstart","text":"From the Julia REPL: ","category":"page"},{"location":"","page":"Quickstart","title":"Quickstart","text":"using Pkg\nPkg.add(\"BioSequenceMappings\")\nusing BioSequenceMappings","category":"page"},{"location":"#Usage","page":"Quickstart","title":"Usage","text":"","category":"section"},{"location":"#Filter-sequences-by-Hamming-distance","page":"Quickstart","title":"Filter sequences by Hamming distance","text":"","category":"section"},{"location":"","page":"Quickstart","title":"Quickstart","text":"Load an alignment, find all sequences with a Hamming distance smaller than 66% to the first one, create a new alignment object from them and save it to a file. ","category":"page"},{"location":"","page":"Quickstart","title":"Quickstart","text":"using BioSequenceMappings\nA = read_fasta(\"../../example/PF00014.fasta\");\nsize(A) # 100 sequences of length 53\ns0 = A[1]; \ns0' # transpose for legibility\nindices = findall(s -> hamming(s, s0; normalize=true) < 0.66, A)\nB = subsample(A, indices)\nwrite(\"PF00014_subsample.fasta\", B)","category":"page"},{"location":"#Remove-columns-with-too-many-gaps","page":"Quickstart","title":"Remove columns with too many gaps","text":"","category":"section"},{"location":"","page":"Quickstart","title":"Quickstart","text":"Load an alignment, find columns with more than 5% gaps and remove them, then create a new alignment object. ","category":"page"},{"location":"","page":"Quickstart","title":"Quickstart","text":"using BioSequenceMappings # hide\nA = read_fasta(\"../../example/PF00014.fasta\"); # uses the default alphabet `Alphabet(:aa)`\ngap_digit = A.alphabet('-') \nf1 = site_specific_frequencies(A)\nnon_gapped_colums = findall(j -> f1[gap_digit, j] <= 0.05, 1:size(f1, 2))' # transpose for legibility\nB = Alignment(A.data[non_gapped_colums', :], A.alphabet; A.names) # sequences are stored as columns","category":"page"}]
}
