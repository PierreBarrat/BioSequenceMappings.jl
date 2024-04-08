@testset "Iterating" begin
    data = Matrix{Int16}([1 2 3; 2 3 4])
    X = Alignment(data; verbose=false)

    @test size(X) == (2, 3)
    @test length(X) == 3
    @test sum(X) == Int16[6, 9] # sum of sequences, meaningless but it's for test
    @test map(sum, X) == Int16[3, 5, 7] # sum of each sequence, also meaningless
    @test findall(x -> x==[1,2], X) == [1]

    @test first(Iterators.reverse(X)) == [3, 4]
end

@testset "Indexing" begin
    data = Matrix{Int16}([1 2 3; 2 3 4])
    X = Alignment(data; verbose=false)

    @test X[1] == Int16[1, 2]
    @test X[1:2] == X.data[:, 1:2]
    @test X[:] == X.data
    @test X[end] == Int16[3, 4]
end

@testset "Subsample" begin
    data = Matrix{Int16}([1 2 3; 2 3 4])
    X = Alignment(data; verbose=false)

    Y = subsample(X, 1)
    @test typeof(Y) == typeof(X)
    @test size(Y) ==(2, 1)
    @test Y.alphabet == X.alphabet
    @test Y[1] == X[1]
    @test Y[1] !== X[1] # different in memory

    Y = subsample(X, 1:2)
    @test typeof(Y) == typeof(X)
    @test size(Y) == (2, 2)
    @test Y.alphabet == X.alphabet
    @test Y[1] == X[1] && Y[2] == X[2]
    @test Y[1] !== X[1] && Y[2] !== X[2]

    Y = subsample(X, 3:-1:2)
    @test typeof(Y) == typeof(X)
    @test size(Y) == (2, 2)
    @test Y[1] == X[3] && Y[2] == X[2]

    Y = subsample(X, :)
    @test all(i -> X[i] == Y[i], eachindex(Y))
end
