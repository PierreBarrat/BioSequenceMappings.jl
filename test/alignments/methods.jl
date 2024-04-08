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
