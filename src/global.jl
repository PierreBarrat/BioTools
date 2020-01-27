## FASTA headers
const augur_header_reference = ["strain", "Strain", "STRAIN", "name", "Name", "label", "Label"]
const augur_header_date = ["date", "Date"]
const augur_header_virus = ["virus", "Virus"]
const augur_minimal_fields = [:strain, :date, :virus]
# Fieldnames below are IGNORED when reading a fasta file
const ignored_header_fields = ["?", "", '?']
const flu_usual_header_fields = ["strain", "virus", "", "date", "region", "country", "", "", "", "segment"]
const augur_all_header_fields = ["strain", "virus", "isolate_id", "date", "region", "country", "division", "location", "passage", "authors", "age", "gender"]

## Errors
function unknown_seqtype()
	@error "Unknown sequence type `$seqtype`-- Possible values: (`:aa`, `:dna`, `:rna`)"
end

