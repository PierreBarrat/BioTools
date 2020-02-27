"""
	compute_strains_fitness!(traj::FrequencyTraj, fp::FluPop, strainfield=:lbi, trajfield=strainfield; shift = :mean)
	compute_strains_fitness!(traj::Array{<:FrequencyTraj,1}, fp::FluPop, strainfield=:lbi, trajfield=strainfield; shift = :mean)
"""
function compute_strains_fitness!(traj::FrequencyTraj, fp::FluPop, strainfield=:lbi, trajfield=strainfield; shift = :mean)
	traj.data[trajfield] = zeros(Float64, length(traj.strains))
	for (i, strains) in enumerate(traj.strains)
		current_date = traj.date + traj.t[i]
		db = find_datebin(current_date, fp)
		# Strains in trajectory
		Ntraj = count(x->!ismissing(fp.strains[x][strainfield]), strains)
		Ftraj = mapreduce(x->ismissing(fp.strains[x][strainfield]) ? 0. : fp.strains[x][strainfield], +, strains, init=0.) 
		# All strains in the datebin
		N = count(x->!ismissing(x[strainfield]), fp.datebin[db])
		F = mapreduce(x->ismissing(x[strainfield]) ? 0. : x[strainfield], +, fp.datebin[db], init=0.) 
		if ismissing(F) || ismissing(Ftraj)
			traj.data[trajfield][i] = missing # This would be the proper way to do it. 
		elseif Ntraj == 0 || Ntraj == N
			traj.data[trajfield][i] = 0.
		elseif N > Ntraj
			if shift==:exclusive_mean
				traj.data[trajfield][i] = Ftraj/Ntraj  -  (F - Ftraj)/(N - Ntraj)
			elseif shift==:mean
				traj.data[trajfield][i] = Ftraj/Ntraj - F/N
			elseif shift==:none
				traj.data[trajfield][i] = Ftraj/Ntraj
			else
				@error "Possible values of `shift` are `:none`, `:mean`, `:exclusive_mean`."
			end
		else
			@error "Ntraj $(Ntraj) > N $N - More strains strains in trajectory than in datebin."
		end
	end
end
function compute_strains_fitness!(traj::Array{<:FrequencyTraj,1}, fp::FluPop, strainfield=:lbi, trajfield=strainfield; shift = :mean)
	for t in traj
		compute_strains_fitness!(t, fp, strainfield, trajfield, shift=shift)
	end
end