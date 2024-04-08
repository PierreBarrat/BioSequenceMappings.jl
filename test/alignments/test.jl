basedir = dirname(@__FILE__)

@testset "Basics" begin
    include(joinpath(basedir, "basics.jl"))
end

@testset "Interfaces" begin
    include(joinpath(basedir, "interfaces.jl"))
end

@testset "Methods" begin
    include(joinpath(basedir, "methods.jl"))
end
