"""
	bin_by_date!fsp::FluPop; start=:auto, last=:auto, binwidth=Day(121), binspacing=Day(121))

Bin `fp` by dates, from `start` to `last`. 
"""
function bin_by_date!(fp::FluPop; 
	start=:auto,
	stop=:auto,
	binwidth=Day(121),
	binspacing=Day(121))
	
	fp.datebin = Dict{Tuple{Date,Date}, Array{AbstractStrain,1}}()
	# Start and end dates
	if start == :auto
		startdate = findmin([x[:date] for x in values(fp.strains)])[1]
	else
		startdate = Date(start)
	end
	if stop == :auto
		stopdate = findmax([x[:date] for x in values(fp.strains)])[1]
	else
		stopdate = Date(stop)
	end
	# 
	now = startdate + binwidth
	while now < stopdate # If dlat - dstart != 0 mod[binwidth], the last bin will be smaller than others and is not considered here
		fp.datebin[now-binwidth, now] = Array{Strain,1}(undef, 0)
		now += binspacing
	end
	binlist = collect(keys(fp.datebin))

	for s in values(fp.strains)
		for b in binlist
			if s[:date] >= b[1] && s[:date] < b[2]
				push!(fp.datebin[b[1],b[2]], s)
			end
		end
	end

	nothing
	# for (d, s) in out
	# 	md = datebin_to_date(d)
	# 	sp.datepop[md] = s
	# end
end

"""
	filter_by_region!(fp::FluPop, r::Array{<:AbstractString,1})
	filter_by_region!(fp::FluPop, r::AbstractString)
"""
function filter_by_region!(fp::FluPop, r::Array{<:AbstractString,1})
	for S in values(fp.datebin)
		idx = findall(s-> !in(s[:region], r), S)
		for i in idx
			delete!(fp.strains, S[i][:strain])
		end
		deleteat!(S, idx)
	end
end
filter_by_region!(fp::FluPop, r::AbstractString) = filter_by_region!(fp, [r])
"""
	filter_by_region(fp::FluPop, r)
"""
filter_by_region(fp::FluPop, r) = begin out = deepcopy(fp); filter_by_region!(out, r); return out end

"""
	filter_by_country!(fp::FluPop, r::Array{<:AbstractString,1})
	filter_by_country(fp::FluPop, r::AbstractString)
"""
function filter_by_country!(fp::FluPop, r::Array{<:AbstractString,1})
	for S in values(fp.datebin)
		idx = findall(s-> !in(s[:country], r), S)
		for i in idx
			delete!(fp.strains, S[i][:strain])
		end
		deleteat!(S, idx)
	end
end
filter_by_country(fp::FluPop, r::AbstractString) = filter_by_country!(fp, [r])
"""
	filter_by_country(fp::FluPop, r)
"""
filter_by_country(fp::FluPop, r) = begin out = deepcopy(fp); filter_by_country!(out, r); return out end