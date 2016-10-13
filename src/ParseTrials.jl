
function ParseTrials(k)


	#Try = JSON.parsefile(string("TrialX$k.json"));
	Try = JSON.parsefile(string("./Trials/TrialX$k.json"));

	allW = Array(Float64, size(Try[1][1][1])[1], size(Try[1][1])[1], size(Try)[1]);
	allT = Array(Float64, size(Try[1][1][1])[1], size(Try[1][1])[1], size(Try)[1]);

	allH = Array(Float64, size(Try[1][2][1])[1], size(Try[1][2])[1], size(Try)[1]);

	allCost = Array(Float64, size(Try)[1]);

	for i=1:size(Try)[1]

		allCost[i] = Try[i][5];
		for s=1:size(Try[1][1])[1]

			allW[:,s,i] =  Try[i][1][s];
			allT[:,s,i] =  Try[i][3][s];
		end

		for g=1:size(Try[1][2])[1]

			allH[:,g,i] =  Try[i][2][g];
		end

	end

	#///////////////////  Normalize H and W  ///////////////////////////
	for c=1:size(allH,3)	
		for i=1:size(allH,1)
		    sumH=sum(allH[i,:,c]);
		    allH[i,:,c] = allH[i,:,c]./sumH;
		    allW[:,i,c] = allW[:,i,c].*sumH;
		end
	end

	Vars = Array(Float64,size(allH,3));
	for i=1:length(Vars)
		Vars[i] = sqrt( 1-Try[i][4]);
	end


	return allW, allH, allT, allCost , Vars

end
