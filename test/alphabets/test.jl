# Very basic tests
@testset "Alphabet constructor tests" begin
    # Test valid input
    @testset "Constructor String + Dict" begin
        alphabet = Alphabet("ACGT-", Dict{Char, Int}('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, '-'=>5))
        @test alphabet.characters == collect("ACGT-")
        @test alphabet.char_to_index == Dict{Char, Int}('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, '-'=>5)
        @test alphabet.index_to_char == Dict{Int, Char}(1=>'A', 2=>'C', 3=>'G', 4=>'T', 5=>'-')
        @test isa(alphabet, Alphabet{Char,Int})
    end

    @testset "Constructor Vec + Dict - custom struct for characters" begin
        struct T c::Char end
        alphabet = Alphabet([T('A'), T('C')], Dict{T, Int}(T('A')=>1, T('C')=>2))
        @test alphabet.characters == [T('A'),T('C')]
        @test alphabet.char_to_index == Dict{T, Int}(T('A')=>1, T('C')=>2)
        @test alphabet.index_to_char == Dict{Int, T}(1=>T('A'), 2=>T('C'))
        @test isa(alphabet, Alphabet{T,Int})
    end

    # Test constructor with Dict input
    @testset "Constructor Dict" begin
        alphabet = Alphabet(Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4))
        @test alphabet.characters == collect("ACGT")
        @test alphabet.char_to_index == Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4)
    end

    # Test constructor with different lengths for string and mapping
    @test_throws AssertionError Alphabet("-ACGT", Dict('A'=>1, 'C'=>2, 'G'=>3))

    # Test constructor with incomplete mapping
    @test_throws AssertionError Alphabet("-ACGT", Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, 'X'=>5))

    # Test constructor with empty string
    @testset "Empty string" begin
        alphabet = Alphabet("")
        isempty(alphabet.characters) && alphabet.char_to_index == Dict()
    end

    # Test constructor from string alone
    @testset "String constructor" begin
        alphabet = Alphabet("AB")
        @test typeof(alphabet) == Alphabet{Char, Int}
        @test alphabet.characters == collect("AB")
        @test alphabet.char_to_index == Dict{Char, Int}('A' => 1, 'B' => 2)
    end


    # Test constructor from string alone -- not integer types
    @testset "String constructor - Int8" begin
        alphabet = Alphabet("AB", Int8)
        @test typeof(alphabet) == Alphabet{Char,Int8}
        @test alphabet.characters == collect("AB")
        @test alphabet.char_to_index == Dict{Char, Int8}('A' => 1, 'B' => 2)
    end

    # the test below wants to test overflow of Int8 for strings with diverse characters
    @testset "Unicode string + Int8" begin
        long_unicode_string = "≔⪉➰Ⅰ⚃⻘⧓⌏☷∉⻔⛴⿐ⱞ➒℃➷⮡⦊⻌⁆ⷢⅽ➫ⰹⓑⲂ⎛➻⮝ⷕ⮥⇑∇ⳣ⍹\u20c2⾣Ⲥⵀ➶⒰⨼⾧Ⓥⷵ〉⣀♈⽺⺆⦤⟹⁍⯕\u2fe4ⳘⱣ⡭⯐⪘℅⺈ⷽ⠴⇽⼝⿂⋅ⵥ␠⬯Ɀ⢒⢏⽹Ⓝⲵ⤢ⶕ⌲ↆ⑂⤻⚲\u2fe8▰⹋⦮⬄≢⣔⫿⅑Ⲏ❠\u2065▢⊺ⵡℙ⠰ⴈⱯ⎽☊ₜ⢳⧷⽌⯸⎾⡀⡚⍁⾱⨷⺊ℯ⼨⸱⸃⠷⛠⾢⤮ⱀ⒡⟽ⷃ⒋⑲⤆⑯⬊ⲁ☫\u2d29⬬⥄╟ⶮ ⢋♼⧂⧽⾰⫖⧼⡕⤦⎂Ⓨ◎⽩➛▍⮄␟♛ⵉ⺭▔╯ℂ✉⌃⤅K⽙⸐Ⓜ⣩⦘Ⱂℇ∞⽐ⴏ⮻☠┓⪡⾶ⓢ⻀≲┬⣅⦓⧈⤬≱⚀⫠⑾⅌"
        @test_throws InexactError Alphabet(long_unicode_string, Int8)

        alphabet = Alphabet(long_unicode_string, Int16)
        @test typeof(alphabet) == Alphabet{Char, Int16}
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
    B = Alphabet(:dna, Int8) # same data and different type
    C = Alphabet(:aa, Int16) # different data and same type

    @test A1 == A2
    @test hash(A1) == hash(A2)
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
        @test convert(Alphabet{Char, Int}, original_alphabet) === original_alphabet
        @test convert(Int, original_alphabet) === original_alphabet

        # Convert to a different type should return a new object with the same content
        converted_alphabet = convert(Alphabet{Char,UInt8}, original_alphabet)
        @test converted_alphabet.characters == original_alphabet.characters
        @test converted_alphabet.char_to_index == original_alphabet.char_to_index
        @test converted_alphabet.index_to_char == original_alphabet.index_to_char

        # Check that the converted alphabet is not the same object as the original
        @test converted_alphabet !== original_alphabet

        # Check the second formulation of convert
        @test convert(UInt8, original_alphabet) == converted_alphabet
    end
end




@testset "Defaults" begin
    local BSM = BioSequenceMappings
    @test Alphabet(:dna).characters == collect(BSM._DEFAULT_NT_ALPHABET_STRING)
    @test Alphabet(:nt).characters == collect(BSM._DEFAULT_NT_ALPHABET_STRING)
    @test Alphabet(:nucleotide).characters == collect(BSM._DEFAULT_NT_ALPHABET_STRING)

    @test Alphabet(:aa).characters == collect(BSM._DEFAULT_AA_ALPHABET_STRING)
    @test Alphabet(:amino_acids).characters == collect(BSM._DEFAULT_AA_ALPHABET_STRING)
    @test Alphabet(:aminoacids).characters == collect(BSM._DEFAULT_AA_ALPHABET_STRING)
    @test Alphabet(:AA, Int8).characters == collect(BSM._DEFAULT_AA_ALPHABET_STRING)

    @test Alphabet(:spin).characters == collect(BSM._DEFAULT_BINARY_ALPHABET_STRING)
    @test Alphabet(:binary).characters == collect(BSM._DEFAULT_BINARY_ALPHABET_STRING)
    @test typeof(Alphabet(:spin, Int8)) == Alphabet{Char, Int8}
    @test_throws ArgumentError Alphabet(:some_symbol)

    @test default_alphabet(4).characters == collect(BSM._DEFAULT_NT_ALPHABET_STRING_NOGAP)
    @test default_alphabet(21).characters == collect(BSM._DEFAULT_AA_ALPHABET_STRING)
    @test default_alphabet(5).characters == collect(BSM._DEFAULT_NT_ALPHABET_STRING)
    @test default_alphabet(14).characters == collect(BSM._DEFAULT_AA_ALPHABET_STRING)
    @test_throws ArgumentError default_alphabet(22)
    @test_throws ArgumentError default_alphabet(1)
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

    @test_throws ArgumentError begin
        s = "-TGCAB"
        alphabet = Alphabet(:dna, Int8)
        alphabet(s)
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

    @test_throws ErrorException begin
        X = [1,2,3,4,5,6]
        alphabet = Alphabet("-ACGT")
        alphabet(X)
    end
    # Test with UTF8
    @test begin
        X = [1,2,3,4]
        alphabet = Alphabet("τ4ν2")
        alphabet(X) == "τ4ν2"
    end
end


