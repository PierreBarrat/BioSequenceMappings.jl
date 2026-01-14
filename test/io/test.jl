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

    @testset "safe reading" begin
        # A normal amino acid alignment - the safe kwarg should not change anything
        fasta_normal = joinpath(basedir, "toy_fasta_aa.fasta")
        aln_1 = read_fasta(fasta_normal; safe=false)
        aln_2 = @test_nowarn read_fasta(fasta_normal; safe=true)
        @test aln_1 == aln_2

        # An alignment wih the X symbol - the safe kwarg should trigger a warning
        # note that if X appeared in the first sequences (default 5), the alphabet
        # would be auto-determined using these sequences and would contain X
        fasta_aa_X = joinpath(basedir, "toy_fasta_aa_X.fasta")
        @test_throws ArgumentError read_fasta(fasta_aa_X)
        @test_logs (:warn, r"Could not read")
    end
end
