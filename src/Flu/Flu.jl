module Flu

using BioTools
using TreeTools
using Dates, DelimitedFiles
using BioSequences

export FluPop, AAFluPop
export PosEvo
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
function FluPop(S::Array{<:AbstractStrain,1})
	return FluPop(strains = Dict(String(x[:strain])=>x for x in S))
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

mutable struct PosEvo{A}
	i::Int64
	alphabet::Array{A,1}
	data::Dict{Tuple{Date,Date}, SiteFrequency{A}}
end
function PosEvo(A::DataType, i::Int64) 
	return PosEvo(i, Array{A,1}(undef,0), Dict{Tuple{Date,Date}, SiteFrequency{A}}())
end


include("IO.jl")
include("global.jl")
include("filtering.jl")
include("frequencies.jl")
include("lbi.jl")

end