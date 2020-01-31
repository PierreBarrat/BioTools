module Flu

using BioTools
using Dates

# Should I implement a flu strain that has segment / flustrain fields and a ntseq field? 
# For now, I'll do this for AAs only. 
# One way might also be to do codon alphabets in BioSequence 
mutable struct FluPop
	strains::Dict{String, Strain}
	datebin::Dict{Tuple{Date,Date},Array{Strain,1}}
end
function FluPop(; strains = Dict{String,Strain}(), 
				datebin = Dict{Tuple{Date,Date},Array{Strain,1}}())
	return FluPop(strains, datebin)
end



include("IO.jl")

end