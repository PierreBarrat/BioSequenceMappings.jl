using BioSequenceMappings
using DCAUtils # to compare weights
using Random
using Test

here = dirname(@__FILE__)

@testset "BioSequenceMappings.jl" begin
    @testset "Alphabet" begin
        println("# Alphabet")
        include(joinpath(here, "./alphabets/test.jl"))
    end
    @testset "Alignment" begin
        println("# Alignment")
        include(joinpath(here, "./alignments/test.jl"))
    end
    @testset "IO" begin
        println("# IO")
        include(joinpath(here, "./io/test.jl"))
    end
end
