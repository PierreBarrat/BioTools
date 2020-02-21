

"""
	abstract type AbstractStrain end

The only requirement for concrete implementations of `AbstractStrain` is the existence of two fields: `seq` and `data::Dict`. `data` is often expected to contain a key `"strain"` referring to the name of the strain by default, but I should try not to enforce it. 
"""
abstract type AbstractStrain end

getindex(s::AbstractStrain, key::String) = get(s.data, key) do 
	haskey(s.data, Symbol(key)) ? s.data[Symbol(key)] : @error "Keys $key or $(Symbol(key)) not found."
end
getindex(s::AbstractStrain, key::Symbol) = get(s.data, key) do 
	haskey(s.data, String(key)) ? s.data[String(key)] : @error "Keys $key or $(String(key)) not found."
end

"""
	mutable struct Strain{A} <: AbstractStrain

Fields: 
- `seq::BioSequence{A}`
- `data::Dict`

	Strain(seq, dat, seqtype)

Values of `seqtype`: `aa`, `dna` or `rna`. 
"""
mutable struct Strain{A} <: AbstractStrain
	seq::LongSequence{A}
	data::Dict
end
function Strain(seqtype)
	if seqtype == :aa
		return Strain(LongAminoAcidSeq(), Dict())
	elseif seqtype == :dna
		Strain(LongDNASeq(), Dict())
	elseif seqtype == :rna
		Strain(LongRNASeq(), Dict())
	else
		unknown_seqtype()
	end

end
function Strain(seq, dat, seqtype ; verbose=false)
	if seqtype == :aa
		try 
			return Strain(LongAminoAcidSeq(seq), dat)
		catch
			cannot_read(seq, dat, seqtype)
			return Strain(seqtype)
		end
	elseif seqtype == :dna
		try
			return Strain(LongDNASeq(seq), dat)
		catch
			cannot_read(seq, dat, seqtype)
			return Strain(seqtype)
		end
	elseif seqtype == :rna
		return Strain(LongRNASeq(seq), dat)
	else
		unknown_seqtype()
	end
end
function cannot_read(seq, dat, seqtype)
	@warn "Could not read the following sequence:"
	println("$seq")
	println("Associated data:\n$dat")
	println("Called seq. type: $seqtype")
	println()
	return nothing
end

"""
	isempty(st::Strain)
"""
function isempty(st::Strain)
	return isempty(st.seq) && isempty(st.data)
end



# abstract type SequenceProfile end
# For now I don't see an efficient application of profile outside of this. 
# For DCA like stuff, it's a float vector, and it's highly unpractical to have something both general and efficient here. 
# So no abstract type

# One way would be to have `data` be an array of subtype `colProfile`.
# `colProfile` could then be designed independently
# As long as the getindex function is overloaded properly this could be practical.
# mutable struct SiteFrequency{A}
# 	i::Int64
# 	M::Int64
# 	alphabet::Array{A,1}
# 	freq::Dict{A,Float64}
# end
# getindex(F::SiteFrequency, a) = F.freq[a]
# counts(F::SiteFrequency, a) = round(Int64, F.freq[a] * F.M)
"""
	mutable struct Profile{A}
		data::Array{Dict{A, Float64},1}
		M::Int64
	end
Access `data` by indexing: `profile[i,a]` --> `profile.data[i][a]`.
"""
mutable struct Profile{A}
	data::Array{Dict{A, Float64},1}
	M::Int64 # Number of sequences the profile is based on
end
function Profile(T::DataType, L::Int64)
	return Profile([Dict{T, Float64}() for i in 1:L], 0)
end

getindex(P::Profile, i::Int64, a) = get(P.data[i], a, 0.)
getindex(P::Profile, i::Int64) = P.data[i]
length(P::Profile) = length(P.data)
iterate(P::Profile, n=1) = iterate(P.data, n)

"""
	Profile(S::Array{<:BioSequence{A}}) where A
	Profile(S::Array{<:Strain{A}}) where A
"""
function Profile(S::Array{<:BioSequence{A}}) where A
	if isempty(S)
		@error "Cannot build a profile  from empty alignment"
	else
		prof = Profile(eltype(A), length(S[1]))
		prof.M = length(S)
		for (m,s) in enumerate(S)
			for (i,a) in enumerate(s)
				prof.data[i][a] = prof[i,a] + 1. /prof.M
			end
		end
	end
	return prof
end
function Profile(S::Array{<:Strain{A}}) where A
	M = length(S)
	if isempty(S)
		@error "Cannot build a profile from empty alignment"
	else
		prof = Profile(eltype(A), length(S[1].seq))
		prof.M = length(S)
		for (m,s) in enumerate(S)
			for (i,a) in enumerate(s.seq)
				prof.data[i][a] = prof[i,a] + 1. /prof.M
			end
		end
	end
	return prof
end