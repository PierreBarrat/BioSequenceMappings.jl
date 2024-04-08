# ###################################################################################
# #################################### MappedAln ####################################
# ###################################################################################

# """
#     mutable struct MappedAln

# ```
#     dat::Matrix{T}
#     q::T = 21
#     mapping::Dict{T, Char} = default_mapping(q)
#     weights::Vector{Float64} = ones(size(dat,1))/size(dat,1) # phylogenetic weights of sequences
#     names::Vector{String} = fill("", size(dat, 1))
# ```

# Stores sequences or a sample of a DCA model. `dat` stores sequences/samples in *columns*:
# `eachcol(X.dat)` will iterate over sequences.

# `mapping` is a dictionary mapping integers in `dat` to `Char` (amino acids, nucleotides, etc...).
# A string can be passed instead of a `Dict`, in which case `DCATools.compute_mapping` will
# be called.

# The constraints on the fields are:
# - `dat` must contain strictly positive integers, smaller than `q`.;
# - the length of `mapping` must be equal to `q` or to `0` (empty mapping is okay);
# - `weights` must sum to one and have positive elements.

# **Important**: When built from a matrix, will *transpose* the input; if `size(dat) = (M, L)`,
# `X=MappedAln(dat)` will return an object with `size(X.dat) = (L, M)`. In other words, assumes
# that the input matrix has sequences as rows.

# ## Methods

# - `getindex(X::MappedAln, i)` returns a matrix/vector `X.dat[:, i]`.
# - `for s in X::MappedAln` iterates over sequences.
# - `eachsequence(X::MappedAln)` returns an iterator over sequences (`Vector{Int}`).
# - `eachsequence_weighted(X::MappedAln)` returns an iterator over sequences and weights.
# - `subsample(X::MappedAln, i)` constructs the subsample defined by index `i`.

# """
# @kwdef mutable struct MappedAln{T<:Integer}
#     dat::Matrix{T}
#     q::T = 21
#     mapping::Dict{T, Char} = default_mapping(q)
#     weights::Vector{Float64} = ones(size(dat,1))/size(dat,1)
#     names::Vector{String} = ["$i" for i in 1:size(dat, 1)]

#     function MappedAln(dat, q, mapping, weights, names)
#         # dat
#         @assert all(>(0), dat) "data must be a matrix/vector of strictly positive integers - got a zero or negative value"
#         @assert maximum(dat) <= q "values in `dat` must be smaller than `q`=$q - got $(maximum(dat))"
#         # mapping
#         @assert isempty(mapping) || q == length(mapping) "Inconsistent size for mapping $mapping and q=$q.
#         Use `mapping=Dict()` if you do not care about the mapping."
#         # weights
#         @assert length(weights) == size(dat, 1) "inconsistent number of weights: $(size(dat,1)) sequence and $(length(weights)) weights"
#         @assert all(>(0), weights) "Weights cannot be negative"
#         @assert isapprox(sum(weights), 1) "Weights must sum to 1"
#         @assert length(names) == size(dat, 1) "Must have as many names as sequences"

#         new(Matrix(dat'), q, mapping, weights, names)
#     end
# end
