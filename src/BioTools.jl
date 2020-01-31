module BioTools

using FastaIO
using BioSequences
using Dates

import BioSequences: @aa_str, @dna_str, @rna_str
export AbstractStrain, Strain

include("types.jl")
include("sequences.jl")
include("IO.jl")
include("filtering.jl")
include("global.jl")

include("Flu/Flu.jl")


end