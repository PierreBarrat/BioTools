"""
	all_trajectories(posh::PosEvo; fixed_thr = 0.95, lost_thr = 0.05, keep_unfinished=false)

Find all active trajectories for given position `posh`. 
An active trajectory is defined by being above fixing threshold (resp. below absence threshold) for 2 timebins in a row, and below (resp. above) at the next timebin. 
Those can then be filtered for given characteristics, such as being seen in some frequency bin, being absent or fixed in the past, etc...  

## Technical note
`start_indices`: `Y[i-1,a]>=fixed_thr && Y[i,a]>=fixed_thr` (*i.e.* mutation was fixed) AND `Y[i+1,a]<fixed_thr, *i.e.* it is now not fixed. This and the reverse for `lost_thr`.   
`end_indices`: `Y[i,a] >= fixed_thr && Y[i+1,a] >= fixed_thr` (*i.e.* mutation is now fixed) AND `Y[i-1,a] < fixed_thr` (*i.e.* it was not fixed before). 
"""
function all_trajectories(posh::PosEvo{A}; fixed_thr = 0.95, lost_thr = 0.05, keep_unfinished=false) where A
	X,Y,N = frequency_series(posh)
	out = Array{FrequencyTraj{A}, 1}(undef, 0)
	for a in 1:size(Y,2)
		#
		start_indices = findall(collect(2:length(X)-1)) do i 
			if Y[i-1,a] >= fixed_thr && Y[i,a] >= fixed_thr # Was fixed
				if Y[i+1,a] < fixed_thr # Now not fixed
					return true
				else
					return false
				end
			elseif Y[i-1,a] <= lost_thr && Y[i,a] <= lost_thr # Was lost
				if Y[i+1,a] > lost_thr  # Now exists
					return true
				else
					return false
				end
			else
				return false
			end
		end .+ 1
		#
		end_indices = findall(collect(2:length(X)-1)) do i 
			if Y[i,a] >= fixed_thr && Y[i+1,a] >= fixed_thr # Now fixed
				if Y[i-1,a] < fixed_thr # Was not fixed
					return true
				else
					return false
				end
			elseif Y[i,a] <= lost_thr && Y[i+1,a] <= lost_thr # Now lost
				if Y[i-1,a] > lost_thr # Was not lost
					return true
				else
					return false
				end
			else
				return false
			end
		end .+ 1


		# Fixing problems with boundaries : trajectories that we did not see start 
		while !isempty(end_indices) && !isempty(start_indices) && end_indices[1] < start_indices[1] 
			splice!(end_indices,1)
		end
		# ... or that we will not see end.
		if !isempty(end_indices) && !isempty(start_indices) && end_indices[end] < start_indices[end]
			if keep_unfinished
				push!(end_indices, length(X)-1)
			elseif !isempty(end_indices) && !isempty(start_indices)
				splice!(start_indices, length(start_indices))
			end
		elseif isempty(end_indices) && !isempty(start_indices)
			if length(start_indices)>1
				error("Empty `end_indices` and more than one in `start_indices`")
			elseif keep_unfinished
				push!(end_indices, length(X)-1)
			else
				start_indices=Int64[]
			end
		end



		# If nothing found, break
		if isempty(start_indices) || isempty(end_indices)
			continue
		end
		# Storing found trajectory
		for (is,ie) in zip(start_indices, end_indices)
			val = posh.alphabet[a]
			date = X[is]
			t = X[is:ie] .- X[is]
			freq = Y[is:ie,a]
			pop = N[is:ie]
			index = Dict(:start=>1, :end=>ie-is+1	, :active=>missing)
			fixation = :poly
			if Y[ie,a] >= fixed_thr
				fixation = :fixed
			elseif Y[ie,a] <= lost_thr
				fixation = :lost
			end
			push!(out, FrequencyTraj(posh.i, val, date, t, freq, pop, index, fixation, Dict()))
		end
	end
	return out
end