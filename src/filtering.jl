"""
	remove(strains::Array{<:AbstractStrain}, names, [field=:strain]; verbose=false)

Remove strains for which `s.data[field]` is in `names`. 
"""
function remove(strains::Array{<:AbstractStrain}, names, field=:strain; verbose=false)
	out = Array{eltype(strains),1}(undef, length(strains))
	idx = 1
	for s in strains
		if !in(s[field], names)
			out[idx] = s
			idx += 1
		end 
	end
	verbose && println("Filtered $(length(strains) - idx + 1) strains")
	deleteat!(out, idx:length(strains))
end
"""
	remove!(strains::Array{<:AbstractStrain}, names, [field=:strain]; verbose=false)

Remove strains for which `s.data[field]` is in `names`. 
"""
function remove!(strains::Array{<:AbstractStrain}, names, field=:strain; verbose=false)
	indices = Int64[]
	for (i,s) in enumerate(strains)
		if in(s[field], names)
			push!(indices, i)
		end 
	end	
	verbose && println("Filtered $(length(indices)) strains")
	deleteat!(strains, indices)
end

"""
	gapfilter(s::Strain; threshold=0.1)

Return true if `s.seq` has less than `threshold * length(s.seq)` gaps.
"""
function gapfilter(s::Strain; threshold=0.1)
	return countgaps(s.seq) < threshold * length(s.seq)
end

hasdate(x::AbstractStrain) = (!ismissing(get(x.data, "date", missing)) || !ismissing(get(x.data, :date, missing)))

"""
	datebin_to_date(d)
"""
function datebin_to_date(d)
	if length(d) != 2 || d[2] < d[1]
		@error "Invalid date bin"
	end
	return d[1] + div(d[2] - d[1], 2)
end
"""
	date_to_datebin(date, binwidth)
"""
function date_to_datebin(date, binwidth)
	return (date-binwidth, date+binwidth)
end