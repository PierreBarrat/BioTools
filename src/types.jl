import Base.isempty
export AbstractStrain, Strain

"""
	abstract type AbstractStrain end

The only requirement for concrete implementations of `AbstractStrain` is the existence of two fields: `seq` and `data::Dict`. `data` is often expected to contain a key `"strain"` referring to the name of the strain by default, but I should try not to enforce it. 
"""
abstract type AbstractStrain end

getindex(s::AbstractStrain, key::String) = get(s.data, key, s.data[Symbol(key)])
getindex(s::AbstractStrain, key::Symbol) = get(s.data, key, s.data[String(key)])

"""
	mutable struct Strain{A} <: AbstractStrain

Fields: 
- `seq::BioSequence{A}`
- `data::Dict`

	Strain(seq, dat, seqtype)

Values of `seqtype`: `aa`, `dna` or `rna`. 
"""
mutable struct Strain{A} <: AbstractStrain
	seq::BioSequence{A}
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
		return Strain(LongAminoAcidSeq(seq), dat)
	elseif seqtype == :dna
		try
			return Strain(LongDNASeq(seq), dat)
		catch
			@goto cannot_read
		end
	elseif seqtype == :rna
		return Strain(LongRNASeq(seq), dat)
	else
		unknown_seqtype()
	end
	@label cannot_read 
	begin
		if verbose
			@warn "Could not read the following sequence:"
			println("$seq")
			println("Associated data:\n$dat")
			println()
		end
		return Strain(seqtype)
	end
end

"""
	isempty(st::Strain)
"""
function isempty(st::Strain)
	return isempty(st.seq) && isempty(st.data)
end