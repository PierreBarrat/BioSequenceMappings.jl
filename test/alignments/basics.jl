@testset "Constructors without alphabet" begin
    @testset "Basic" begin
        data = Matrix{Int}([1 2 3; 2 3 4])
        alignment = Alignment(data; verbose=false)
        @test typeof(alignment) == Alignment{Char, Int} # should autofind nt alphabet
        @test size(alignment.data) == (2, 3)
        @test length(alignment.weights) == 3
        @test length(alignment.names) == 3
        @test Alphabet(alignment) == BioSequenceMappings.default_alphabet(4, Int)
    end

    @testset "With type" begin
        data = Matrix{Int8}([1 2 3; 2 3 4])
        alignment = Alignment(data; verbose=false)
        @test typeof(alignment) == Alignment{Char,Int8}
        @test Alphabet(alignment) == BioSequenceMappings.default_alphabet(4, Int8)
    end

    @testset "Without alphabet" begin
        data = Matrix{Int8}([1 2 3; 2 3 4])
        alignment = Alignment(data; alphabet = nothing)
        @test typeof(alignment) == Alignment{Nothing, Int8}
        @test isnothing(alignment.alphabet)
    end
end

@testset "Constructors with explicit alphabet" begin
    data = Matrix{Int16}([1 2 3; 2 3 4])
    @testset "Basic" begin
        alphabet = Alphabet("ACGT", Int16)
        alignment = Alignment(data, alphabet)
        @test typeof(alignment) == Alignment{Char, Int16}
        @test Alphabet(alignment) == alphabet
    end

    @testset "Mismatched types - Alphabet gets priority" begin
        alphabet = Alphabet("ACGT", Int8)
        alignment = Alignment(data, alphabet)
        @test typeof(alignment) == Alignment{Char, Int8}
        @test Alphabet(alignment) == alphabet
    end

    @testset "From alphabet constructor" begin
        alignment = Alignment(data, "BCDE")
        @test typeof(alignment) == Alignment{Char, Int16}
        @test Alphabet(alignment) == Alphabet("BCDE", Int16)
    end

    @testset "With keywords" begin
        weights = [1/4, 1/4, 1/2]
        names = ["A", "B", "CD"]
        alignment = Alignment(data, :dna; weights, names)
        @test Alphabet(alignment) == Alphabet(:dna, Int16)
        @test alignment.weights == weights
        @test alignment.names == names
    end

    @testset "From vector of vectors" begin
        vecdat = Matrix{Int}([1 2 3; 2 3 4]) |> eachcol |> collect
        alignment = Alignment(vecdat, :dna)
        @test typeof(alignment) == Alignment{Char, Int}
        @test size(alignment) == (2, 3) # three sequences of length 2
        @test Alphabet(alignment) == Alphabet(:dna, Int)
    end

    @testset "From vector of integers" begin
        sequence = Int16[1,2,3,4,5]
        alignment = Alignment(sequence, :dna)
        @test typeof(alignment) == Alignment{Char, Int16}
        @test size(alignment) == (5, 1)
        @test Alphabet(alignment) == Alphabet(:dna, Int16)
    end
end

@testset "Constructor with keyword alphabet" begin
    data = Matrix{Int}([1 2 3; 2 3 4])

    @testset "Explicit alphabet" begin
        alphabet = Alphabet("ACGT", Int)
        alignment = Alignment(data; alphabet)
        @test typeof(alignment) == Alignment{Char, Int}
        @test Alphabet(alignment) == alphabet
    end

    @testset "Symbol" begin
        alignment = Alignment(data; alphabet=:dna)
        @test Alphabet(alignment) == Alphabet(:dna)
    end

    @testset "Vector of vector data" begin
        data = [[1,2], [3,4]]
        alignment = Alignment(data; alphabet = :dna)
        @test size(alignment) == (2,2)
        @test alignment[1] == [1,2]
        @test Alphabet(alignment) == Alphabet(:dna)
    end

    @testset "Vector of integers data" begin
        data = [3,4]
        alignment = Alignment(data; alphabet = :dna)
        @test size(alignment) == (2,1)
        @test alignment[1] == [3,4]
        @test Alphabet(alignment) == Alphabet(:dna)
    end
end

@testset "Test for errors" begin
    data = Matrix{Int}([1 2 3; 2 3 21])
    @test Alignment(data; verbose=false).alphabet == default_alphabet(21)
    @test_throws AssertionError Alignment(data, :dna)
    @test_throws AssertionError Alignment(data; weights = [1/2, 1/4, 1/3], verbose=false)
    @test_throws AssertionError Alignment(data; weights = [1., 1., -1.], verbose=false)
    @test_throws AssertionError Alignment(data; weights = [1/2, 1/4, 1/4, 0.], verbose=false)
    @test_throws AssertionError Alignment(data; weights = [1/2, 1/4, 1/3], names = ["A"], verbose=false)

    data = Matrix{Int}([1 2 3; 2 3 22])
    @test_throws ErrorException Alignment(data; verbose=false)
    @test Alignment(data; alphabet=:none).alphabet == nothing
end
