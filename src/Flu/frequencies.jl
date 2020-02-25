"""
	PosEvo(fp::FluPop, i::Int64; 
			ambiguous = false, 
			threshold = 0.05)
"""
function PosEvo(fp::FluPop, i::Int64; 
			ambiguous = false, 
			threshold = 0.05)
	A = eltype(first(fp.strains)[2].seq)
	out = PosEvo(A,i)
	#
	if isempty(fp.datebin)
		@error "`fp.datebin` is empty"
	end
	#
	for (d, S) in fp.datebin
		# Frequencies for that site and timebin
		f = SiteFrequency(A, i, 0)
		for s in S
			if ambiguous || !isambiguous(s.seq[i])
				f[s.seq[i]] = get(f, s.seq[i], 0.) + 1.
				f.M += 1
			end
		end
		for a in BioTools.alphabet(f)
			f[a] /= f.M
			if !in(a, out.alphabet)
				push!(out.alphabet, a)
			end
		end
		# 
		out.data[d] = f
	end
	remove_rare_symbols!(out, threshold)
	return out
end
function PosEvo(fp::FluPop; 
				ambiguous=false,
				threshold=0.05)
	ph = []
	for i in 1:length(first(fp.strains)[2].seq)
	    print("$i       \r")
	    push!(ph, PosEvo(fp, i, ambiguous=ambiguous, threshold=threshold))		
	end
	return ph
end
# Regional weights could be implemented at the PosEvo level
# For instance, a function `PosEvo(fp, i, weights)`
# It's the only place where we compute frequencies, and we also know the strains since we have `fp` 

# Why is this useful? I'll leave it empty for now... 
function remove_rare_symbols!(ph::PosEvo, threshold)
	
end

function frequency_series(ph::PosEvo)
	X = Array{Date,1}(undef, 0)
	Y = zeros(Float64, length(ph.data), length(ph.alphabet))
	pop = Array{Int64,1}(undef, 0)
	for (i, (d,f)) in enumerate(sort(ph.data))
		push!(pop, f.M)
		push!(X, BioTools.datebin_to_date(d))
		for (k,a) in enumerate(ph.alphabet)
			Y[i,k] = f[a]
		end
	end
	return (X,Y,pop)
end