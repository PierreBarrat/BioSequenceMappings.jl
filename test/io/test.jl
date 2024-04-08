basedir = dirname(@__FILE__)

@testset "Alphabet from fasta" begin
    auto_alph = BioSequenceMappings.auto_alphabet_from_fasta # convenience
    @test auto_alph(joinpath(basedir, "toy_fasta_aa.fasta")) == Alphabet(:aa)
    @test auto_alph(joinpath(basedir, "toy_fasta_dna.fasta")) == Alphabet(:dna)
    @test auto_alph(joinpath(basedir, "toy_fasta_bin.fasta")) == Alphabet(:binary)
    @test !in(
        auto_alph(joinpath(basedir, "toy_fasta_else.fasta")),
        [Alphabet(:dna), Alphabet(:aa), Alphabet(:binary)]
    )
end
