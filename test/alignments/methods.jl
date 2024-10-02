@testset "weights" begin
    Random.seed!(42)
    Y = Int8.(rand(1:2, 10, 100))
    X = Alignment(Y; verbose=false)
    for θ in 0:.2:1
        (w1, M1) = BioSequenceMappings.compute_weights(Y, θ; normalize=false)
        (w2, M2) = DCAUtils.compute_weights(Y, θ; verbose=false)
        @test w1 ≈ w2
        @test M1 ≈ M2
        @test sum(BioSequenceMappings.compute_weights(Y, .4; normalize=true)[1]) ≈ 1
        @test BioSequenceMappings.compute_weights(X, θ; normalize=false)[1] ≈ w1
    end
end

@testset "find sequence" begin
    @warn "I should write a test for `find_sequence` and `match_sequences`"
end

@testset "Concatenation" begin
    dat_1 = [[1,2], [3,4]]
    names_1 = ["A", "B"]
    aln1 = Alignment(dat_1; alphabet=:nt, names=names_1)

    dat_2 = [[1,2], [3,4]]
    names_2 = ["B", "C"]
    aln2 = Alignment(dat_2; alphabet=:nt, names=names_2)

    @testset "Two arguments" begin
        A = cat(aln1, aln2)
        @test size(A) == (2, 4)
        @test A.names == vcat(names_1, names_2)
        @test A[1] == dat_1[1]
        @test A[3] == dat_2[1]
    end

    @testset "Multiple arguments" begin
        A = cat(aln1, aln2, aln1)
        @test size(A) == (2, 6)
        @test A.names == vcat(names_1, names_2, names_1)
        @test all(≈(1/6; rtol = 1e-10), A.weights)
        @test A[end] == dat_1[end]
    end

    @testset "Errors" begin
        aln2 = Alignment(dat_2; alphabet=:aa)
        @test_throws ErrorException cat(aln1, aln2)

        dat_3 = [[1,2,3]]
        aln3 = Alignment(dat_3; alphabet=:nt)
        @test_throws ErrorException cat(aln1, aln3)
    end
end
