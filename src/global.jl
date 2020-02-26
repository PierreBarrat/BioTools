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
## 
## Errors
function unknown_seqtype()
	@error "Unknown sequence type `$seqtype`-- Possible values: (`:aa`, `:dna`, `:rna`)"
end

