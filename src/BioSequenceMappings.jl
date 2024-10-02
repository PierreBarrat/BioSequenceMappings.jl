module BioSequenceMappings

using FASTX
using OneHotArrays
using Random

import Base: length, size
import Base: in, ==, hash, convert, copy
import Base: getindex, firstindex, lastindex, eachindex, view, keys
import Base: iterate, eltype
import Base: unique
import Base: write

include("alphabet.jl")
export Alphabet
export default_alphabet, symbols, translate

include("alignment.jl")
export Alignment, AbstractAlignment
export eachsequence, sequence_length, sequence_number, subsample, subsample_random
export find_sequence, match_sequences, named_sequences

include("weights.jl")
export compute_weights, compute_weights!

include("IO.jl")
export read_fasta

include("misc.jl")
export hamming, pairwise_hamming

include("statistics.jl")
export site_specific_frequencies, consensus, pairwise_frequencies, pairwise_correlations

end
