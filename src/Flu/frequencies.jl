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
	A = eltype(first(fp.strains)[2].seq)
	ph = PosEvo{A}[]
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

"""
	frequency_series(ph::PosEvo)

Return arrays `X`, `Y` and `Z`, with `X` containing dates, `Y` a matrix with columns containing frequencies of symbols, and `Z` a vector containing population sizes frequencies are based on. 
"""
function frequency_series(ph::PosEvo; freq_threshold = 0.05)
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
	# Delete positions with frequencies always smaller than `freq_threshold`
	if freq_threshold > 0.
		idx = Int64[]
		for i in 1:size(Y,2)
			for t in 1:size(Y,1)
				if Y[t,i] > freq_threshold
					push!(idx, i)
					break
				end
			end
		end
	else
		idx = 1:length(Y,2)
	end
	return (X,Y[:,idx],pop)
end

"""
	entropy(Z::PosEvo)

Return a dictionary `Date => Float` giving the entropy as a function of time for position history `Z`. 
"""
function entropy(Z::PosEvo; freq_threshold = 0.05)
	X,Y,tmp = frequency_series(Z, freq_threshold=freq_threshold)
	S = zeros(Float64, length(X))
	for i in 1:length(S)
		f = Y[i,:] ./ sum(Y[i,:])
		S[i] = isnan(StatsBase.entropy(f)) ? 0. : StatsBase.entropy(f)
	end
	return S,X
end