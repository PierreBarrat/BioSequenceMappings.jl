# Very basic tests
@testset "Alphabet constructor tests" begin
    # Test valid input
    @testset begin
        alphabet = Alphabet("ACGT-", Dict{Char, Int}('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, '-'=>5))
        @test alphabet.string == "ACGT-"
        @test alphabet.char_to_index == Dict{Char, Int}('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, '-'=>5)
        @test alphabet.index_to_char == Dict{Int, Char}(1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'-')
        @test isa(alphabet, Alphabet{Int})
    end

    # Test constructor with AbstractDict input
    @test begin
        alphabet = Alphabet(Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4))
        alphabet.string == "ACGT" && alphabet.char_to_index == Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4)
    end

    # Test constructor with different lengths for string and mapping
    @test_throws AssertionError Alphabet("-ACGT", Dict('A'=>1, 'C'=>2, 'G'=>3))

    # Test constructor with incomplete mapping
    @test_throws AssertionError Alphabet("-ACGT", Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, 'X'=>5))

    # Test constructor with empty string and mapping
    @test begin
        alphabet = Alphabet("")
        alphabet.string == "" && alphabet.char_to_index == Dict()
    end

    # Test constructor from string alone
    @test begin
        alphabet = Alphabet("AB")
        typeof(alphabet) == Alphabet{Int} && alphabet.string == "AB" && alphabet.char_to_index == Dict{Char, Int}('A' => 1, 'B' => 2)
    end


    # Test constructor from string alone -- not integer types
    @test begin
        alphabet = Alphabet("AB", Int8)
        typeof(alphabet) == Alphabet{Int8} && alphabet.string == "AB" && alphabet.char_to_index == Dict{Char, Int8}('A' => 1, 'B' => 2)
    end

    # the test below wants to test overflow of Int8 for strings with diverse characters
    @test begin
        long_unicode_string = "≔⪉➰Ⅰ⚃⻘⧓⌏☷∉⻔⛴⿐ⱞ➒℃➷⮡⦊⻌⁆ⷢⅽ➫ⰹⓑⲂ⎛➻⮝ⷕ⮥⇑∇ⳣ⍹\u20c2⾣Ⲥⵀ➶⒰⨼⾧Ⓥⷵ〉⣀♈⽺⺆⦤⟹⁍⯕\u2fe4ⳘⱣ⡭⯐⪘℅⺈ⷽ⠴⇽⼝⿂⋅ⵥ␠⬯Ɀ⢒⢏⽹Ⓝⲵ⤢ⶕ⌲ↆ⑂⤻⚲\u2fe8▰⹋⦮⬄≢⣔⫿⅑Ⲏ❠\u2065▢⊺ⵡℙ⠰ⴈⱯ⎽☊ₜ⢳⧷⽌⯸⎾⡀⡚⍁⾱⨷⺊ℯ⼨⸱⸃⠷⛠⾢⤮ⱀ⒡⟽ⷃ⒋⑲⤆⑯⬊ⲁ☫\u2d29⬬⥄╟ⶮ ⢋♼⧂⧽⾰⫖⧼⡕⤦⎂Ⓨ◎⽩➛▍⮄␟♛ⵉ⺭▔╯ℂ✉⌃⤅K⽙⸐Ⓜ⣩⦘Ⱂℇ∞⽐ⴏ⮻☠┓⪡⾶ⓢ⻀≲┬⣅⦓⧈⤬≱⚀⫠⑾⅌"
        @test_throws InexactError Alphabet(long_unicode_string, Int8)

        alphabet = Alphabet(long_unicode_string, Int16)
        typeof(alphabet) == Alphabet{Int16}
    end

    # Test throw for abstract type
    @test_throws AssertionError Alphabet("ABC", Signed)

    # Test throw for non unique alphabet string
    @test_throws AssertionError Alphabet("AAA")

    # Test constructor with non-Char keys in mapping
    @test_throws MethodError Alphabet("-ACGT", Dict('A'=>1, 'C'=>2, 'G'=>3, 4=>4))
end

@testset "Equality and hash" begin
    A1 = Alphabet(:dna, Int16)
    A2 = Alphabet(:dna, Int16) # same data and same type
    A3 = Alphabet("-ACGT", Int16) # same data and same type
    B = Alphabet(:dna, Int8) # same data and different type
    C = Alphabet(:aa, Int16) # different data and same type

    @test A1 == A2 == A3
    @test hash(A1) == hash(A2) == hash(A3)
    @test A1 != B && hash(A1) != hash(B)
end

@testset "Copy, convert" begin
    # Testing the `copy` function
    @testset "copy" begin
        original_alphabet = Alphabet("ACGT")
        copied_alphabet = copy(original_alphabet)

        # Check that the copied alphabet is equal to the original (using data)
        @test copied_alphabet == original_alphabet

        # Check that the copied alphabet is not the same object as the original (in memory)
        @test copied_alphabet !== original_alphabet
    end

    # Testing the `convert` function
    @testset "convert function" begin
        original_alphabet = Alphabet("ACGT")

        # Convert to the same type should return the same object with === (Base.convert)
        @test convert(Alphabet{Int}, original_alphabet) === original_alphabet
        @test convert(Int, original_alphabet) === original_alphabet

        # Convert to a different type should return a new object with the same content
        converted_alphabet = convert(Alphabet{UInt8}, original_alphabet)
        @test converted_alphabet.string == original_alphabet.string
        @test converted_alphabet.char_to_index == original_alphabet.char_to_index
        @test converted_alphabet.index_to_char == original_alphabet.index_to_char

        # Check that the converted alphabet is not the same object as the original
        @test converted_alphabet !== original_alphabet

        # Check the second formulation of convert
        @test convert(UInt8, original_alphabet) == converted_alphabet
    end
end




@testset "Defaults" begin
    @test Alphabet(:dna).string == BioSequenceMappings._DEFAULT_NT_ALPHABET_STRING
    @test Alphabet(:nt).string == BioSequenceMappings._DEFAULT_NT_ALPHABET_STRING
    @test Alphabet(:nucleotide).string == BioSequenceMappings._DEFAULT_NT_ALPHABET_STRING

    @test Alphabet(:aa).string == BioSequenceMappings._DEFAULT_AA_ALPHABET_STRING
    @test Alphabet(:amino_acids).string == BioSequenceMappings._DEFAULT_AA_ALPHABET_STRING
    @test Alphabet(:aminoacids).string == BioSequenceMappings._DEFAULT_AA_ALPHABET_STRING
    @test Alphabet(:AA, Int8).string == BioSequenceMappings._DEFAULT_AA_ALPHABET_STRING

    @test Alphabet(:spin).string == BioSequenceMappings._DEFAULT_BINARY_ALPHABET_STRING
    @test Alphabet(:binary).string == BioSequenceMappings._DEFAULT_BINARY_ALPHABET_STRING
    @test typeof(Alphabet(:spin, Int8)) == Alphabet{Int8}
    @test_throws ErrorException Alphabet(:some_symbol)

    @test default_alphabet(4).string == BioSequenceMappings._DEFAULT_NT_ALPHABET_STRING_NOGAP
    @test default_alphabet(21).string == BioSequenceMappings._DEFAULT_AA_ALPHABET_STRING
    @test default_alphabet(5).string == BioSequenceMappings._DEFAULT_NT_ALPHABET_STRING
    @test default_alphabet(14).string == BioSequenceMappings._DEFAULT_AA_ALPHABET_STRING[1:14]
    @test_throws ErrorException default_alphabet(22)
end

@testset "Sequence to int" begin
    @test begin
        s = "TGCA"
        alphabet = Alphabet(:dna)
        X = alphabet(s)
        X == [5,4,3,2] && typeof(X) == Vector{Int}
    end

    @test begin
        s = "TGCA"
        alphabet = Alphabet(:dna, Int8)
        X = alphabet(s)
        X == [5,4,3,2] && typeof(X) == Vector{Int8}
    end

    @test begin
        s = "-TGCAB"
        alphabet = Alphabet(:dna, Int8)
        X = alphabet(s)
        X == [1,5,4,3,2,1] && typeof(X) == Vector{Int8}
    end

    @test begin
        s = "AB12"
        alphabet = Alphabet("CBA321", UInt8)
        X = alphabet(s)
        X == [3, 2, 6, 5] && typeof(X) == Vector{UInt8}
    end
end

@testset "Int to sequence" begin
    @test begin
        X = Int[1,2,3,4]
        alphabet = Alphabet(:dna)
        alphabet(X) == "-ACG"
    end

    @test begin
        X = Int8[5,4,3,2]
        alphabet = Alphabet(:dna, Int16) # Int16 should not matter here
        alphabet(X) == "TGCA"
    end

    @test_throws KeyError begin
        X = [1,2,3,4,5,6]
        alphabet = Alphabet(:dna)
        alphabet(X)
    end
    # Test with UTF8
    @test begin
        X = [1,2,3,4]
        alphabet = Alphabet("τ4ν2")
        alphabet(X) == "τ4ν2"
    end
end


