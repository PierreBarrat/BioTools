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
# Why is this useful? I'll leave it empty for now... 
function remove_rare_symbols!(ph::PosEvo, threshold)
	
end