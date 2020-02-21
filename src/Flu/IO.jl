"""
	FluPop(f::Union{AbstractString,IO}, sequencetype::Symbol, headerfields; 
	flulineage=missing, segment=missing, strainfilters = [!is_flu_outlier(flulineage), BioTools.hasdate], separator = '|')

Call `readfastastrains` to read `f`. Store the result in a `FluPop` object. 
"""
function FluPop(f::Union{AbstractString,IO}, 
	sequencetype::Symbol, 
	headerfields; 
	flulineage=missing, 
	segment=missing, 
	strainfilters = [!is_flu_outlier(flulineage), BioTools.hasdate],
	separator = '|')	
	
	strains = readfastastrains(f, sequencetype, headerfields, separator = separator, strainfilters = strainfilters)
	return FluPop(strains = Dict(String(x[:strain])=>x for x in strains))
end
"""
	AAFluPop(f::Union{AbstractString,IO}, headerfields; [kwargs...])

Call `FluPop(f, :aa, headerfields, [kwargs...])
"""
function AAFluPop(f::Union{AbstractString,IO}, headerfields; 
	flulineage=missing, 
	segment=missing, 
	strainfilters = [!is_flu_outlier(flulineage), BioTools.hasdate],
	separator = '|')
	return FluPop(f, :aa, headerfields, flulineage=flulineage, segment=segment, strainfilters=strainfilters, separator=separator)
end

is_flu_outlier(x::AbstractStrain, flulineage) = in(x[:strain], outliers[flulineage])
function is_flu_outlier(flulineage)
	return x->is_flu_outlier(x,flulineage)
end
