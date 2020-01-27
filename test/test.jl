using Test
using BioSequences
using BioTools
using Dates

@testset "`writefasta` basics" begin 
	s = Strain(aa"ACGT", Dict(:strain=>"mystrain", :date=>today(), :somefield=>"someval")
		)
	t = Strain(aa"ACGT", Dict(:strain=>"myotherstrain", "secretfield"=>"secret"))
	@test writefasta(s, [:strain, :date]) == ("mystrain|2020-01-17|", "ACGT")
	@test writefasta(s, [:strain, :secretfield], fillvals = true) == ("mystrain|?|", "ACGT")
	@test writefasta(s, [:strain, "somefield"]) == ("mystrain|someval|", "ACGT")
	@test writefasta([s,t], [:strain, :secretfield], fillvals=true) == [("mystrain|?|", "ACGT"), ("myotherstrain|secret|", "ACGT")]
end

