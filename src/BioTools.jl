module BioTools

using FastaIO
using BioSequences, BioSymbols
using Dates

# import BioSequences: @aa_str, @dna_str, @rna_str
import Base.getindex, Base.get, Base.setindex!, Base.keys, Base.values
import Base.isempty, Base.length, Base.enumerate, Base.iterate
import Base.show
import BioSequences.alphabet



export AbstractStrain, Strain
export SiteFrequency, Profile
export Flu

include("types.jl")
include("sequences.jl")
include("IO.jl")
include("filtering.jl")
include("global.jl")

include("Flu/Flu.jl")


end