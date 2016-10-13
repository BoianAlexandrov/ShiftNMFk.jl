function Plot(X, maxSource)

	#X = readcsv("./Input/Observation.csv");

	if isdir("./Results") == false
		println("ERROR: Results directory does not exist!")
	else

		Trials = collect(1:maxSource);

		Norm = Array(Float64, size(Trials));
		Silhouette = Array(Float64, size(Trials));

		for z = 1:length(Trials)

			numberOfProcesses = Trials[z];

			#///////// This code computes the reconstruction from our averaged answers ///////

			#=Rec = readcsv("./Results/Reconstruction$numberOfProcesses.csv");
			#normalize X & Rec
			for i=1:size(X,1)
			    sumX=sum(X[i,:]);
			    sumRec=sum(Rec[i,:]);
			    X[i,:] = X[i,:]./sumX;
			    Rec[i,:] = Rec[i,:]./sumRec;
		  	end

			cost=0.5*vecnorm(X-Rec,2)^2;
			=#

			cost = readcsv("./Results/Cost$numberOfProcesses.csv");
			Norm[z] = mean(cost);

			avgStability = readcsv("./Results/avgStability$numberOfProcesses.csv");
			Silhouette[z] = mean(avgStability);
			println(Norm[z]);

		end

		writecsv("Norm.csv", Norm);
		writecsv("Silhouette.csv", Silhouette);

		x = Trials;

		ResultPlot = plot(
		  layer(x=x, y=Norm, Geom.point, Geom.line, Theme(default_color=colorant"green")),
		  layer(x=x, y=Silhouette, Geom.point,  Geom.line, Theme(default_color=colorant"red")),
		  Guide.xlabel("Number of Sources"), Guide.ylabel("Fro. Norm and Silhouette Value"), Guide.title("Results"), Guide.manual_color_key("Legend", 
		  	["Silhouette Value", "Norm Minimization"], ["red", "green",])
		);

		draw(SVG("ResultPlot.svg", 25cm, 20cm), ResultPlot); 
	end

	return Silhouette, Norm

end