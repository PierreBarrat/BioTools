module Flu

using BioTools
using TreeTools
using Dates, DelimitedFiles

export FluPop, AAFluPop
export bin_by_date!
# Should I implement a flu strain that has segment / flustrain fields and a ntseq field? 
# For now, I'll do this for AAs only. 
# One way might also be to do codon alphabets in BioSequence 
mutable struct FluPop{T<:AbstractStrain}
	strains::Dict{String, T}
	datebin::Dict{Tuple{Date,Date},Array{T,1}}
end
function FluPop(; strains = Dict{String,Strain}(), 
				datebin = Dict{Tuple{Date,Date},Array{eltype(values(strains)),1}}()) 
	return FluPop(strains, datebin)
end
"""
	add_strain_field!(fp::FluPop, field, default=0.)
"""
function add_strain_field!(fp::FluPop, field, default=0.)
	for s in values(fp.strains)
		s.data[field] = default
	end
	nothing
end

"""
	consensus(fp::FluPop) 
"""
consensus(fp::FluPop) = consensus([x for x in values(fp.strains)])



include("IO.jl")
include("global.jl")
include("filtering.jl")
include("lbi.jl")

end