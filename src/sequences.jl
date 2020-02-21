export hamming, profile, consensus

"""
	countgaps(s::BioSequence)
"""
function countgaps(s::BioSequence)
	return count(x->isgap(x)||isambiguous(x), s)
end
countgaps(s::AbstractStrain) = countgaps(s.seq)

function hamming_ignore_ambiguous(a::BioSequence, b::BioSequence)
	h = 0
	for (x,y) in zip(a,b)
		if x != y && !isambiguous(x) && !isambiguous(y)
			h += 1
		end
	end
	return h
end
function hamming_ambiguous_mismatch(a::BioSequence, b::BioSequence)
	h = 0
	for (x,y) in zip(a,b)
		if x != y 
			h += 1
		end
	end
	return h	
end
"""
	hamming(a::BioSequence, b::BioSequence, [normalization=false]; ambiguous=:ignore)
	hamming(a::AbstractStrain, b::AbstractStrain [, normalization=false]; ambiguous=:ignore)

`ambiguous` keyword:
- `:ignore` : ambiguous symbols are considered equal to anything .
- `:mismatch` : ambiguous symbols behave like regular symbols.
"""
function hamming(a::BioSequence, b::BioSequence, 
	normalization=false; 
	ambiguous=:ignore)
	if ambiguous == :ignore && normalization
		return hamming_ignore_ambiguous(a,b)/length(a)
	elseif ambiguous == :ignore && !normalization
		return Float64(hamming_ignore_ambiguous(a,b))
	elseif ambiguous == :mismatch && normalization
		return hamming_ambiguous_mismatch(a,b)/length(a)
	elseif ambiguous == :mismatch && !normalization
		return Float64(hamming_ambiguous_mismatch(a,b))
	else
		@error "Incorrect arguments."
	end
end
hamming(a::AbstractStrain, b::AbstractStrain, normalization=false; ambiguous=:ignore) = hamming(a.seq, b.seq, normalization, ambiguous=ambiguous)
"""
	hamming(P::Profile{A}, Q::Profile{A}, normalization=false; 
				ambiguous=:ignore) where A <: BioSymbol
"""
function hamming(P::Profile{A}, Q::Profile{A}, normalization=false; 
				ambiguous=:ignore) where A <: BioSymbol
	length(P) != length(Q) && @error "Profiles of different lengths"
	out = 0.
	for (i,(p,q)) in enumerate(zip(P,Q))
		val = 1.
		for a in BioSequences.alphabet(A)
			if ambiguous == :ignore && isambiguous(a) 
				val -= get(p,a,0.) + get(q,a,0.) - get(p,a,0.)*get(q,a,0.)
			else
				val -= get(p,a,0.)*get(q,a,0.)
			end
		end
		out += val
	end
	if normalization
		out /= length(P)
	end
	return out
end
function hamming(P::Profile, s::BioSequence, normalization = false; ambiguous = :ignore)
	length(P) != length(s) && @error "Sequence and profile of different lengths"
	out = 0.
	for (i,(a,p)) in enumerate(zip(s,P))
		out += 1.
		for b in keys(p)
			if b == a || (ambiguous==:ignore && (isambiguous(a)||isambiguous(b)))
				out -= p[b]
			end
		end
	end
	if normalization 
		out /= length(P)
	end
	return out
end
hamming(s::BioSequence, P::Profile) = hamming(P,s)
hamming(s::Strain, P::Profile) = hamming(s.seq, P)
hamming(P::Profile, s::Strain) = hamming(s.seq, P)


"""
	consensus(P::Profile{A}) where A<:BioSymbol
	consensus(X::Array{<:BioSequence,1})

Return the consensus *sequence*, of type `LongSequence{alphabet}` where `eltype(alphabet)==A`.

	consensus(X::Array{Strain} [, data = Dict(:strain=>"consensus")])

Return consensus `Strain`. 
"""
function consensus(P::Profile{A}) where A<:BioSymbol
	seq = LongSequence{symbol_to_alphabet[A]}(length(P))
	for (i,p) in enumerate(P)
		seq[i] = findmax(p)[2]
	end
	return seq
end
consensus(X::Array{<:BioSequence,1}) = consensus(Profile(X))
function consensus(X::Array{<:Strain,1}, data = Dict(:strain=>"consensus"))
	p = Profile(X)
	seq = consensus(p)
	return Strain(seq, data)
end

### OLD
# """
# 	consensus(Q::Array{LongSequence{A},1}) where A
# 	consensus(Q::Array{<:Strain,1})
# """
# function consensus(Q::Array{LongSequence{A},1}) where A
# 	profQ = profile(Q)
# 	seq = []
# 	for p in profQ
# 		push!(seq, findmax(p)[2])
# 	end
# 	return LongSequence{A}(seq)
# end
# consensus(Q::Array{<:Strain,1}) = Strain(consensus([x.seq for x in Q]), Dict(:strain=>"consensus"))


# """
# 	profile(A::Array{BioSequence,1})
# 	profile(Q::Array{<:AbstractStrain,1})
# """
# function profile(A::Array{<:BioSequence,1})
# 	out = [Dict() for i in 1:length(A[1])]
# 	for (m,s) in enumerate(A)
# 		for (i,a) in enumerate(s)
# 			out[i][a] = get(out[i], a, 0) + 1
# 		end
# 	end
# 	return out
# end
# function profile(Q::Array{<:AbstractStrain,1})
# 	return profile([x.seq for x in Q])
# end



