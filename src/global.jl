# Declaration ... what is the clean way to do that? 
parse_date() = nothing
parse_float() = nothing

## FASTA headers
const augur_header_reference = ["strain", "Strain", "STRAIN", "name", "Name", "label", "Label"]
const augur_header_date = ["date", "Date"]
const augur_header_virus = ["virus", "Virus"]
const augur_minimal_fields = [:strain, :date, :virus]
# Fieldnames below are IGNORED when reading a fasta file
const ignored_header_fields = ["?", "", '?',:?]
const special_fields = ["date", :date, "Date", :date, :fitness]
const parse_special_field = Dict("date" => parse_date, :date => parse_date, :Date => parse_date, "Date" => parse_date, :fitness =>parse_float)
## The following could be added to a BioSequences fork
# const symbol_to_alphabet = Dict(BioSequences.AminoAcid => BioSequences.AminoAcidAlphabet,
								# BioSequences.DNA => BioSequences.DNAAlphabet,
								# BioSequences.RNA => BioSequences.RNAAlphabet)
symbol_to_alphabet(::Type{AminoAcid}) = BioSequences.AminoAcidAlphabet
symbol_to_alphabet(::Type{DNA}) = BioSequences.DNAAlphabet{4}
symbol_to_alphabet(::Type{RNA}) = BioSequences.RNAAlphabet{4}
ambiguous(::Type{AminoAcid}) = AA_X
ambiguous(::Type{AminoAcidAlphabet}) = AA_X
ambiguous(::Type{DNA}) = DNA_N
ambiguous(::Type{DNAAlphabet{4}}) = DNA_N
ambiguous(::Type{RNA}) = RNA_N
ambiguous(::Type{RNAAlphabet{4}}) = RNA_N
isambiguous(::T) where T<:Real = false # Artificial sequences are never ambiguous
#
alphabet_to_seqsymbol(::Type{AminoAcidAlphabet}) = :aa
alphabet_to_seqsymbol(::Type{DNAAlphabet{4}}) = :dna
alphabet_to_seqsymbol(::Type{RNAAlphabet{4}}) = :rna
alphabet_to_seqsymbol(::Integer) = :artificial
# datatypes
const sequenceymbols = (:aa, :rna, :dna, :int8, :int64, :bool, :artificial)
function type(x::Symbol)
	if x == :aa return AminoAcidAlphabet
	elseif x == :rna return RNAAlphabet{4}
	elseif x == :dna return DNAAlphabet{4}
	elseif x == :int8 return Int8
	elseif (x == :int64 || x == :artificial) return Int64
	elseif x == :bool return Bool
	else
		@error "Possible symbols: $(sequenceymbols)"
	end
end

#
const _DNAAlphabet = mapreduce(x->Char(x), *, BioSequences.alphabet(DNA))
const _RNAAlphabet = mapreduce(x->Char(x), *, BioSequences.alphabet(RNA))
const _AAAlphabet = mapreduce(x->Char(x), *, BioSequences.alphabet(AminoAcid))

## 
## Errors
function unknown_seqtype()
	@error "Unknown sequence type `$seqtype`-- Possible values: (`:aa`, `:dna`, `:rna`)"
end


const aa2int = Dict(AA_Gap => 1, AA_A=>2, AA_C=>3, AA_D=>4, AA_E=>5, AA_F=>6, AA_G=>7, AA_H=>8, AA_I=>9, AA_K=>10, AA_L=>11, AA_M=>12, AA_N=>13, AA_P=>14, AA_Q=>15, AA_R=>16, AA_S=>17, AA_T=>18, AA_V=>19, AA_W=>20, AA_Y=>21)
const int2aa = Dict(i => findfirst(==(i), BioTools.aa2int) for i in 1:21)
