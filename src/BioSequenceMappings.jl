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
export default_alphabet

include("alignment.jl")
export Alignment
export sequence_length, sequence_number, subsample

include("weights.jl")
export compute_weights

include("IO.jl")
export read_fasta

end
