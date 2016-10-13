function ParseLoc(micPos)

	Try = JSON.parsefile("./LocTrials/LocSols.json")

	ns = length(Try[1][1][1]);
	allPositions = Array(Float64, length(Try[1][1][1]), 2, length(Try));

	allSpeed = Array(Float64, length(Try));
	AllInit = Array(Float64, length(Try[1][3]), length(Try));
	allMinf = Array(Float64, length(Try));
	Center = Array(Any, ns);
	polPos = Array(Any, ns);


	for i = 1:length(Try)
		if length(Try[i]) == 4
			println(i)
			allPositions[:,:,i] = [Try[i][1][1] Try[i][1][2]];
			allSpeed[i] = Try[i][2];
			AllInit[:,i] = Try[i][3];
			allMinf[i] = Try[i][4];
		end


	end

	#to 75% quantile
	idx = find(allMinf.<=quantile(allMinf, .50));
	#allPositions = allPositions[:,:,idx];

	avgPos = Array(Float64, ns, 2);
	posStd = Array(Float64, ns, 2);
	rad = Array(Float64, ns);

	for i=1:ns
		Center[i]=[median(allPositions[i,1,idx])  median(allPositions[i,2,idx])];
		polPos[i] = ((allPositions[i,1,idx] - Center[i][1,1]).^2 + (allPositions[i,2,idx] - Center[i][1,2]).^2).^(1/2);
		
	end

	AllPol = zeros(size(polPos[1]));
	for i=1:ns

		AllPol = AllPol + polPos[i];
	end


	#id = find(quantile(AllPol[1,1,:][:], .25) .< AllPol[1,1,:].<quantile(AllPol[1,1,:][:], .75));
	id = find(AllPol[1,1,:].<=quantile(AllPol[1,1,:][:], .50));
	idx = idx[id];


	v = sqrt(allPositions[:,1,idx].^2 + allPositions[:,2,idx].^2);
	Vstd = std(v,3);
	for i=1:ns
		avgPos[i,:] = [median(allPositions[i,1,idx])  median(allPositions[i,2,idx])];
		posStd[i,:] = [std(allPositions[i,1,idx])  std(allPositions[i,2,idx])];
		#rad[i] = sqrt(posStd[i,1].^2 + posStd[i,2].^2);
	end
	rad[:]=Vstd[:];


	# ///////////////////////////// Plotting /////////////////////////////////////////////////////////////////////////////////////
	#SoPos = readcsv("SourcePosition.csv");
	#micPos = readcsv("MicPosition.csv")

	#=TFD = plot(	layer(x = avgPos[:,1], y = avgPos[:,2], Geom.point, Theme(default_color = colorant"orange")),
				#layer(x = SoPos[:,1], y = SoPos[:,2], Geom.point, Theme(default_color = colorant"black")),
				layer(x = micPos[:,1], y = micPos[:,2], Geom.point, Theme(default_color = colorant"magenta")),
				
				 Guide.manual_color_key("Legend", [ "Real Location of Sources", "Location of Sensors", "FOund Sources"], ["black", "magenta", "orange"]),
				 
		);


	draw(SVG("SourceLocations.svg", 20cm, 20cm), TFD);=#
	writecsv("FoundSourcePosition.csv", avgPos);
	return avgPos;

end
