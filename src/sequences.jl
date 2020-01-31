export hamming

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
	hamming(a::AbstractStrain, b::AbstractStrain, [normalization=false]; ambiguous=:ignore)

`ambiguous` keyword:
- `:ignore` : ambiguous symbols are considered equal to anything .
- `:mismatch` : ambiguous symbols are different from regular symbols.
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
