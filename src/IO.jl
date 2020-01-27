import FastaIO.writefasta
export writefasta, readfastastrains

# Closure for assigning labels to strains
let strain_number = 0
	global new_strain_number() = strain_number+=1
	global reset_strain_number() = strain_number = 0
end
function new_strain_label()
	return "strain_$(new_strain_number())"
end

"""
	writefasta(s::AbstractStrain, fields; fillvals = false)
	writefasta(f::IO, s::AbstractStrain, fields; fillvals=false)
	writefasta(S::Array{<:AbstractStrain}, fields; fillvals = false)
	writefasta(f::String, S::Array{<:AbstractStrain}, fields; fillvals=false, mode="w")

Write strain `s` to a fasta format. Header is built using `fields` from `s.data`. 
"""
function writefasta(s::AbstractStrain, fields; fillvals = false)
	# Make the header
	header = mapreduce(*, fields, init = "") do x
		if haskey(s.data, String(x))
			"$(get(s.data, String(x), 0))|"
		elseif haskey(s.data, Symbol(x))
			"$(get(s.data, Symbol(x), 0))|"
		elseif fillvals
			"?|"
		else
			@error "Field `$x` not in strain"
		end
		# "$(get(s.data, x) do ; fillvals ? "?|" : @error "Field $x not in strain"; end)|"  -- One liner, but does not handle the symbol/string 
	end[1:end-1] # Removing trailing '|'

	return (header, String(s.seq))
end
writefasta(f::IO, s::AbstractStrain, fields; fillvals=false) = writefasta(f, [writefasta(s, fields, fillvals=fillvals)])
function writefasta(f::String, S::Array{<:AbstractStrain}, fields; fillvals=false, mode="w")
	open(f, mode) do io
		for s in S
			writefasta(io, s, fields, fillvals=fillvals)
		end
	end
end
function writefasta(S::Array{<:AbstractStrain}, fields; fillvals = false)
	out = Array{Tuple{String, String},1}(undef, length(S))
	for (i,s) in enumerate(S)
		out[i] = writefasta(s, fields, fillvals = fillvals)
	end
	return out
end

"""
	readfastastrains(f::Union{AbstractString,IO}, sequence_type::Symbol, headerfields; separator = '|', strainfilter=x->true)

Read strains in `f` to an array of `Strain`. Info about arguments: 
- `sequence_type`: `:dna`, `:aa` or `:rna`
- `headerfields`: Array of field names for parsing the fasta headers. `"?"` will be ignored. 
- `separator`: Typically `|` 
- `strainfilter`: Ignore strains `s` for which strainfilter(s) == false`. 
"""
function readfastastrains(f::Union{AbstractString,IO}, sequence_type::Symbol, headerfields; separator = '|', strainfilter=x->true)
	strains = Array{Strain,1}(undef, 0)
	nfiltered = 0
	nunread = 0
	ntot = 0
	typeof(f) <: AbstractString ? println("Reading $f...") : println("Reading alignment...")
	for (n,s) in FastaReader(f)
		dat = parse_header(n, headerfields, separator)
		st = Strain(s, dat, sequence_type)
		if BioTools.isempty(st)
			nunread += 1
		elseif strainfilter(st)
			push!(strains, st)
		else
			nfiltered += 1
		end
		ntot += 1
	end
	println("Read $(length(strains)) strains out of $ntot. Filtered $nfiltered. Could not read $nunread")
	return strains
end

function parse_header(h, headerfields, separator)
	sh = split(h, separator)
	dat = Dict()
	for (i,f) in enumerate(headerfields)
		if !in(f, ignored_header_fields)
			dat[f] = sh[i]
		end
	end
	return dat
end
function gapfilter(s::Strain; threshold=0.1)
	return countgaps(s.seq) < threshold * length(s.seq)
end