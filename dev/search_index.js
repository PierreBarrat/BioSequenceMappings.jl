var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = BioSequenceMappings","category":"page"},{"location":"#BioSequenceMappings","page":"Home","title":"BioSequenceMappings","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BioSequenceMappings.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [BioSequenceMappings]","category":"page"},{"location":"#BioSequenceMappings.Alignment","page":"Home","title":"BioSequenceMappings.Alignment","text":"mutable struct Alignment{T} where T <: Integer\n\n    dat::Matrix{T}\n    alphabet::Union{Nothing, Alphabet{T}}\n    weights::Vector{Float64} = ones(size(dat,1))/size(dat,1) # phylogenetic weights of sequences\n    names::Vector{String} = fill(\"\", size(dat, 1))\n\nBiological sequences as vectors of type T<:Integer. dat stores sequences in columns: size(dat) returns a tuple (L, M) with L the length and M the number of sequences.\n\nalphabet represents the mapping between integers in dat and biological symbols (nucleotides, amino acids...). If nothing, the alignment cannot be mapped to biological sequences.\n\nweights represent phylogenetic weights, and are initialized to 1/M. They must sum to 1. names are the label of sequences, and are expected to be in the same order as the columns of dat. They do not have to be unique, and can be ignored\n\nImportant: When built from a matrix, will transpose the input; if size(dat) = (M, L),\n\nX=Alignment(dat, alphabet) will return an object with size(X.dat) = (L, M).\n\nIn other words, assumes that the input matrix has sequences as rows.\n\nMethods\n\ngetindex(X::Alignment, i) returns a matrix/vector X.dat[:, i].\nfor s in X::Alignment iterates over sequences.\neachsequence(X::Alignment) returns an iterator over sequences (Vector{Int}).\neachsequence_weighted(X::Alignment) returns an iterator over sequences and weights as tuples\nsubaln(X::Alignment, idx) constructs the subaln defined by index idx.\n\n\n\n\n\n","category":"type"},{"location":"#BioSequenceMappings.Alignment-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T<:Integer","page":"Home","title":"BioSequenceMappings.Alignment","text":"Alignment(data::AbstractMatrix{T}; alphabet = :auto, kwargs...)\n\nKeyword argument alphabet can be :auto, :none or nothing, or an input to Alphabet. Other keyword arguments are passed to the default constructor of Alignment.\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.Alignment-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Alphabet{T}}} where T","page":"Home","title":"BioSequenceMappings.Alignment","text":"Alignment(data::AbstractMatrix, alphabet; kwargs...)\n\ndata is a matrix of integers. alphabet can be either     - an Alphabet     - nothing     - something to build an alphabet from: the constructor Alphabet will be called.     (e.g. a symbol like :aa, a string, ...)\n\nIf the types of alphabet and data mismatch, data is converted.\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.Alphabet","page":"Home","title":"BioSequenceMappings.Alphabet","text":"struct Alphabet{T}\n    string::String\n    char_to_index::Dict{Char, T}\n    index_to_char::Dict{T, Char}\nend\n\nAlphabet can be constructed\n\nfrom a String and an optional type: Alphabet(string::String[, T<:Integer])\nfrom a mapping Dict{Char, T} where T<:Integer: Alphabet(mapping)\nfrom a Symbol, using default alphabets\nfrom an integer, using default alphabets (see default_alphabets).\n\n\n\n\n\n","category":"type"},{"location":"#BioSequenceMappings.compute_mapping-Union{Tuple{AbstractString}, Tuple{T}, Tuple{AbstractString, Type{T}}} where T<:Integer","page":"Home","title":"BioSequenceMappings.compute_mapping","text":"compute_mapping(s::AbstractString)\n\nReturn a Dict{Int, Char}: Dict(i => c for (i,c) in enumerate(s)).\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.compute_weights","page":"Home","title":"BioSequenceMappings.compute_weights","text":"compute_weights(X::AbstractAlignment, θ = 0.2; normalize = true)\n\nCompute phylogenetic correction weights for sequences of X.     The weight sequence S is 1/N, where N is the number of sequences in X at     hamming distance less than H from S (including S itself).     The threshold H is floor(θ⋅L) where L is the sequence length.\n\nThe return value is a tuple (weights, Meff), where Meff is the sum of weights (pre-normalization). If normalize, weights are normalized to sum to one. .\n\n\n\n\n\n","category":"function"},{"location":"#BioSequenceMappings.compute_weights!-Tuple{Any, Any}","page":"Home","title":"BioSequenceMappings.compute_weights!","text":"compute_weights!(X, θ; kwargs...)\n\nCompute and set weights for X. See compute_weights.\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.default_alphabet-Union{Tuple{Integer}, Tuple{T}, Tuple{Integer, Type{T}}} where T<:Integer","page":"Home","title":"BioSequenceMappings.default_alphabet","text":"default_alphabet(q::Int, T::Type)\n\nif q==21, amino acids\nif q==5, nucleotides\nif q==4, nucleotides without gaps\nif q==2, binary (0, 1)\nelse, if q<21, return the restriction of amino acids to the first q sites\nif q>21, fails\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.pairwise_hamming-Tuple{Alignment, Alignment}","page":"Home","title":"BioSequenceMappings.pairwise_hamming","text":"pairwise_hamming(X, Y; step=1, step_left, step_right)\n\nReturn matrix of all hamming distances between sequences of X and Y. Only consider sequences every step. step_left and step_right can be used to skip sequence either in X or in Y:\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.pairwise_hamming-Tuple{Alignment}","page":"Home","title":"BioSequenceMappings.pairwise_hamming","text":"pairwise_hamming(X; step, as_vec=true)\n\nVector of pairwise hamming distances of sequences in X, ordered as [H(1,2), H(1,3), ..., H(M-1, M)] with H standing for hamming distance. If as_vec=false, will return a Matrix instead.\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.read_fasta-Tuple{AbstractString}","page":"Home","title":"BioSequenceMappings.read_fasta","text":"read_fasta(fastafile::AbstractString; alphabet = :auto, kwargs...)\nread_fasta(\n    fastafile::AbstractString, alphabet;\n    weights = false, theta = 0.2, verbose = false,\n)\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.site_specific_frequencies-Tuple{Any, Vararg{Any}}","page":"Home","title":"BioSequenceMappings.site_specific_frequencies","text":"site_specific_frequencies(X::AbstractAlignment, weights=nothing)\n\nReturn the site specific frequencies of X. If as_vec (default), the result is a vector of length Lxq. Otherwise, it is a matrix of q rows and L columns\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.subsample-Tuple{BioSequenceMappings.AbstractAlignment, Int64}","page":"Home","title":"BioSequenceMappings.subsample","text":"subsample(X::AbstractAlignment, indices)\n\nReturn an Alignment containing only the sequences of X at indices.\n\n\n\n\n\n","category":"method"},{"location":"#BioSequenceMappings.subsample_random-Tuple{BioSequenceMappings.AbstractAlignment, Int64}","page":"Home","title":"BioSequenceMappings.subsample_random","text":"subsample_random(X::AbstractAlignment, m::Int)\n\nReturn an Alignment with m sequences taking randomly from X. Sampling is done without replacement, meaning the m sequences are all at different positions in X.\n\n\n\n\n\n","category":"method"}]
}
