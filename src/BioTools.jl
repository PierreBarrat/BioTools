module BioTools

using FastaIO
using BioSequences, BioSymbols
using Dates

# import BioSequences: @aa_str, @dna_str, @rna_str
import Base.getindex
import Base.isempty, Base.length, Base.enumerate, Base.iterate

export AbstractStrain, Strain
export Profile
export getindex
export Flu

include("types.jl")
include("sequences.jl")
include("IO.jl")
include("filtering.jl")
include("global.jl")

include("Flu/Flu.jl")


end