module BioTools

using FastaIO
using BioSequences

import BioSequences: @aa_str, @dna_str, @rna_str
export AbstractStrain, Strain

include("types.jl")
include("sequences.jl")
include("global.jl")
include("IO.jl")
include("filtering.jl")

end