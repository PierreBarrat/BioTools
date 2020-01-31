function FluPop(f::Union{AbstractString,IO}, 
	sequencetype::Symbol, 
	headerfields; 
	flustrain=missing, 
	segment=missing, 
	strainfilter = [x->!is_flu_outlier(x, flustrain), hasdate],
	separator = '|')	
	
	strains = readfastastrains(f, sequencetype, headerfields, separator = separator, strainfilter = strainfilter)
	@warn "Flu outliers have not been set up properly. Need to move config file."
	return FluPop(strains = Dict(x[:strain]=>x for x in strains))
end

is_flu_outlier(x::AbstractStrain, flustrain) = in(x[strain], outliers[flulineage])