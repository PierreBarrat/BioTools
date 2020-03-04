## FASTA headers
const augur_header_reference = ["strain", "Strain", "STRAIN", "name", "Name", "label", "Label"]
const augur_header_date = ["date", "Date"]
const augur_header_virus = ["virus", "Virus"]
const augur_minimal_fields = [:strain, :date, :virus]
# Fieldnames below are IGNORED when reading a fasta file
const ignored_header_fields = ["?", "", '?',:?]
const special_fields = ["date", :date]
const parse_special_field = Dict("date" => parse_date, :date => parse_date)
## The following could be added to a BioSequences fork
const symbol_to_alphabet = Dict(BioSequences.AminoAcid => BioSequences.AminoAcidAlphabet,
								BioSequences.DNA => BioSequences.DNAAlphabet,
								BioSequences.RNA => BioSequences.RNAAlphabet)

ambiguous(::Type{AminoAcid}) = AA_X
ambiguous(::Type{DNA}) = DNA_N
ambiguous(::Type{RNA}) = RNA_N

# datatypes
function type(x::Symbol)
	if x == :aa return AminoAcidAlphabet
	elseif x == :rna return RNAAlphabet
	elseif x == :dna return DNAAlphabet
	elseif x == :int8 return Int8
	elseif (x == :int64 || x == :artificial) return Int64
	elseif x == :bool return Bool
	else
		@error "Possible symbols: $([:aa, :rna, :dna, :int8, :int64, :bool, :artificial])"
	end
end
# For artificial sequences  
isambiguous(::T) where T<:Real = false
# Sequence types
# const bioseqs = [:aa, :dna, :rna]

## 
## Errors
function unknown_seqtype()
	@error "Unknown sequence type `$seqtype`-- Possible values: (`:aa`, `:dna`, `:rna`)"
end

