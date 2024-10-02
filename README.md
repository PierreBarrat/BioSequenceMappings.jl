# BioSequenceMappings

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PierreBarrat.github.io/BioSequenceMappings.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PierreBarrat.github.io/BioSequenceMappings.jl/dev/)

The aim of the package is to facilitate the task of converting biological sequences (nucleotides, amino acids) to integers or to onehot representation. 
It also provides simple functions to manipulate alignments or to compute simple statistics. 

*Note*: I do not frequently use the onehot format, so for now it is not documented / not well tested. I'll develop it more when I get the time. 

*Note*: The package is not aimed at being highly computationally efficient when dealing with sequences. The purpose is simply to help interfacing with other algorithms (I have DCA-like models in mind). If you are dealing with very large alignments and very large sequences, it is possible that this does not fit your needs. 

## Installation

From a Julia session: 
```julia
using Pkg
Pkg.add("BioSequenceMappings")
using BioSequenceMappings
```


## Usage

Check the [documentation](https://pierrebarrat.github.io/BioSequenceMappings.jl/dev/) for more information. Below are some examples. 

### Filter sequences by Hamming distance

Load an alignment, find all sequences with a Hamming distance smaller than 66% to the first one, create a new alignment object from them and save it to a file. 

```julia
using BioSequenceMappings
A = read_fasta("example/PF00014.fasta");
size(A) # 100 sequences of length 53
s0 = A[1]; 
indices = findall(s -> hamming(s, s0; normalize=true) < 0.66, A)
B = subsample(A, indices)
write("example/PF00014_subsample.fasta", B)
```

### Remove columns with too many gaps

Load an alignment, find columns with more than 5% gaps and remove them, then create a new alignment object. 

```julia
using BioSequenceMappings
A = read_fasta("example/PF00014.fasta"); # uses the default alphabet Alphabet(:aa)
gap_digit = A.alphabet('-') 
f1 = site_specific_frequencies(A)
non_gapped_colums = findall(j -> f1[gap_digit, j] <= 0.05, 1:size(f1, 2))
B = Alignment(A.data[non_gapped_colums, :], A.alphabet; A.names) # sequences are stored as columns
```
