
function FindLocations(W, T, Tstd, micPos, globIter)
	
	#globIter = 100;
	ns = size(W,2);

	tic()
	Try = Parallel_Tri(globIter, T, Tstd, W, micPos, 1e-8);
	toc()

	if isdir("./LocTrials") == false
		mkdir("./LocTrials");
	end

	f=open("./LocTrials/LocSols.json", "w");
	JSON.print(f,Try);
	close(f);

	allVec= Array(Float64, 3, ns, globIter);

	for i = 1:globIter
		allVec[:,:,i] = [Try[i][1]'; repmat([Try[i][2]], 1, size(T,2))];
	end

	#Idx = cluster_solutions(allVec, 100);

	orderedPos = allVec[:,:,:];

	positions = Array(Float64, size(orderedPos, 2), 2, size(orderedPos, 3));
	for i=1:globIter
		positions[:,:,i] = orderedPos[1:2,:,i]';
		
	end

	
	avgPos = ParseLoc(micPos);
	return avgPos
end
