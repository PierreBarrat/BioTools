function compute_strains_fitness!(traj::FrequencyTraj, fp::FluPop, field=:lbi)
	for (i, strains) in enumerate(traj.strains)
		current_date = traj.date + traj.t[i]
		datebin = find_datebin(current_date, sp)
		
		findall(x -> fp.strains[x][field]>=0., strains)

		Ftraj = mapreduce(x->fp.strains[x][field], +, strains, init=0.)

	end

end