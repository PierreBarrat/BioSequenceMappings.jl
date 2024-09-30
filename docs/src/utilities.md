# Utilities

## Hamming distance

[`hamming(x, y)`](@ref) returns the normalized Hamming distance between two integer vectors. 

## Alignment statistics

### Pairwise Hamming distance

[`pairwise_hamming(A)`](@ref) returns a vector containing all Hamming distances between pairs of sequences in `A`. 
For large alignments, it is often practical to consider only a subset of Hamming distances: the `step=n` keyword can be used to only consider every nth sequence. 
The function is also adapted to two alignments: [`pairwise_hamming(A,B)`] will consider all pairs of sequences with one member in `A` and the other one in `B`. 

### Statistics

- the profile of the alignment: [`site_specific_frequencies(A)`](@ref) returns a `q x L` matrix where `q` is the size of the alphabet and `L` the length of sequences, with element `(a, i)` being the fraction of sequences in which character `a` was found at position `i`. 
- the pairwise frequencies: [`pairwise_frequencies(A)`](@ref)
- the pairwise correlations: [`pairwise_correlations(A)`](@ref)

All these functions accept a vector of weights as a second argument, and will by default use `A.weights` if it is not provided. 

### Weights

In DCA-like algorithms, it is customary to weight sequences by their degree of phylogenetic relations. 
Typically, the weight of a sequence is inversely proportional to the number of other sequences with a Hamming distance smaller than some threshold. 
For an alignment `X`, computing the weights is done using [`compute_weights`](@ref): 
```@repl
using BioSequenceMappings # hide
A = read_fasta("../../example/toy_fasta_dna.fasta")
compute_weights(A) # compute and return the weight vector, as well as the effective number of sequences
compute_weights!(A); # same, but sets the weights field in A
A.weights
```



