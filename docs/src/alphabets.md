```@setup alpha1
using BioSequenceMappings
```
```@meta 
DocTestSetup = quote
	using BioSequenceMappings
	A = Alphabet(['A', 'C', 'G', 'T', '-'])
end
```

# The `Alphabet` type

## Basics

An `Alphabet` contains information necessary to map biological symbols to integers and inversely. 
The full type is `Alphabet{A,I}`, where `A` is the type of biological symbols (typically `Char`) and `I` is a subtype of `Integer`. 

The simplest way to create an alphabet is from a list of symbols: 
```@repl alpha1
A = Alphabet(['A', 'C', 'G', 'T', '-'])
```
The created alphabet `A` associates `Char` correspondings to nucleotides to `Int` according to the index at which they appear in the input vector: `'A' => 1`, `'C' => 2`, etc...
Note that we could have created the same alphabet from a string (since it's based on `Char`s) or from a dictionary: 
```jldoctest alpha2
julia> B = Alphabet("ACGT-");

julia> C = Alphabet(Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4, '-'=>5));

julia> A == B
true

julia> A == C
true
```

The alphabet is used to map symbols to integers and inversely. 
This is done by calling the object directly, as a function: 
```jldoctest alpha2
julia> A('A') # mapping Char to Int
1

julia> A(1) # mapping Int to Char
'A': ASCII/Unicode U+0041 (category Lu: Letter, uppercase)

julia> A("A-GT") # mapping a string to a vector 
4-element Vector{Int64}:
 1
 5
 3
 4

julia> A([1,2,3]) # mapping a vector to a string
"ACG"
```

If needed, the mapping is accessible using the `symbols`function: 
```jldoctest alpha2
julia> symbols(A)
5-element Vector{Char}:
 'A': ASCII/Unicode U+0041 (category Lu: Letter, uppercase)
 'C': ASCII/Unicode U+0043 (category Lu: Letter, uppercase)
 'G': ASCII/Unicode U+0047 (category Lu: Letter, uppercase)
 'T': ASCII/Unicode U+0054 (category Lu: Letter, uppercase)
 '-': ASCII/Unicode U+002D (category Pd: Punctuation, dash)

julia> A |> symbols |> prod # as a string
"ACGT-"
```

## Default alphabets

The package comes with three default alphabets: 
- an amino-acid alphabet `aa_alphabet` using the mapping `"-ACDEFGHIKLMNPQRSTVWY"`;
- a nucleotide alphabet `nt_alphabet` using the mapping `"-ACGT"`;
- a "binary" alphabet `BioSequenceMappings.binary_alphabet`, which I found useful for simulations, with the mapping: `"01"`. 

These can be also be accessed by calling `Alphabet(name)` where `name` is a symbol corresponding to any of the default alphabets. 
The symbolic names can be easily be found:
```jldoctest alpha2
julia> BioSequenceMappings.aa_alphabet_names # also works with nt and binary alphabets
(:aa, :AA, :aminoacids, :amino_acids)

julia> Alphabet(:aa) == aa_alphabet
true

julia> Alphabet(:amino_acids)([1,2,3])
"-AC"
```

Each default alphabet is also associated to a specific cardinality of biological symbols through the function `default_alphabet`. 
This means that an integer vector with elements ranging from 1 to `q` will be associated to the following alphabets: 
```juliadoctest alpha2
julia> default_alphabet(2) == Alphabet(:binary) # q == 2
true

julia> default_alphabet(5) == Alphabet(:nt) # q == 5
true

julia> default_alphabet(21) == Alphabet(:aa) # 5 < q <= 21
true

julia> default_alphabet(15) == Alphabet(:aa) # 5 < q <= 21
```
This association is useful to create `Alignment` objects from a matrix of integers without having to specify the alphabet manually. 

## Default characters

When reading biological sequences, it can be convenient to associate all unexpected characters to a default symbol, for instance the gap. 
This can be achieved by providing the `default_char` keyword argument when constructing the alphabet: 
```@repl alpha1
A_default = Alphabet("ACGT-"; default_char = '-')
A_default("ABCDEF") # 'unknown' chars are mapped to '-', in turn mapped to 5
A("ABCDEF") # if no defaults are provided, fails
```

This also works the other way around: integers that are not in the range of the alphabet are mapped to the default symbol: 
```@repl alpha1
A_default(1:10) # indices larger than 5 are mapped to the gap
A(1:10) # if no defaults are provided, fails
```

## Using specific integer types

When created as above, the alphabet will default to using `Int` as the integer type. 
If dealing with large amounts of data, it can be beneficial to use a more sober type. 
This is done by providing an extra argument of the desired type when constructing the alphabet: 
```@repl alpha1
B = Alphabet("ACGT-", UInt8)
B == A
B("A-")
```


## Translating between alphabets

It often happens to me that I have an integer vector `X` representing a sequence, but with a mapping different from the one I am used to. 
The [`translate`](@ref) function lets me convert it to another integer vector with the right mapping. 
```jldoctest alpha2
julia> strange_alphabet = Alphabet("TCGA-"); # the default is "-ACGT"

julia> X = Int[2, 2, 5, 4, 5]; # representing the sequence "CC-A-" according to the above

julia> strange_alphabet(X)
"CC-A-"

julia> nt_alphabet(X) # this is obviously wrong - nt_alphabet uses "-ACGT"
"AATGT"

julia> Y = translate(X, strange_alphabet, nt_alphabet)
5-element Vector{Int64}:
 3
 3
 1
 2
 1

julia> nt_alphabet(Y) # now this works
"CC-A-"
```