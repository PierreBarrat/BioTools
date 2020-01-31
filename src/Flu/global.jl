outliers = Dict{Union{String,Missing},Array{String}}(missing=>String[])
for l in lineages
	p = dirname(pathof(BioTools)) * "/FluTools/config/outliers_$(l).txt"
	if isfile(p)
		outliers[l] = vec(readdlm(p, String))
	else
		outliers[l] = String[]
	end
	p = dirname(pathof(BioTools)) * "/FluTools/config/my_outliers_$(l).txt"
	if isfile(p)
		append!(outliers[l], vec(readdlm(p,String)))
	end
end