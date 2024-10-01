basedir = dirname(@__FILE__)

@testset "Alphabet from fasta" begin
    auto_alph = BioSequenceMappings.auto_alphabet_from_fasta # convenience
    @test auto_alph(joinpath(basedir, "toy_fasta_aa.fasta")) == Alphabet(:aa)
    @test auto_alph(joinpath(basedir, "toy_fasta_dna.fasta")) == Alphabet(:dna)
    @test auto_alph(joinpath(basedir, "toy_fasta_bin.fasta")) == Alphabet(:binary)
    @info "Tests: warning will be produced below"
    @test !in(
        auto_alph(joinpath(basedir, "toy_fasta_else.fasta")),
        [Alphabet(:dna), Alphabet(:aa), Alphabet(:binary)]
    )
end


@testset "Writing and reading" begin
    # test whether FASTX handles long names
    long_name = repeat("1", 250)
    short_name = "Hello"
    seq = [1,2,3,4]
    A = Alignment([seq, reverse(seq)], alphabet = :nt, names = [long_name, short_name])
    write(joinpath(basedir, "aln.fasta"), A)
    B = read_fasta(joinpath(basedir, "aln.fasta"), alphabet = :nt)
    @test find_sequence(long_name, B) == (1, seq)
    @test find_sequence(short_name, B) == (2, reverse(seq))
end
