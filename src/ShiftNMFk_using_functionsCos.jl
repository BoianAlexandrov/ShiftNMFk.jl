
function shiftNMFk(X, maxSource, globalIter, nmfIter)

	#Trials =  [1, 2, 3, 4, 5];
	Trials = collect(1:maxSource);

	elapsedTime = cell(length(Trials),4);

	#X = readcsv("./Input/Observation.csv");
	inputMatrix = X;
	inputMatrix=inputMatrix';
	numberOfPoints = size(inputMatrix, 1);
	numberOfSamples = size(inputMatrix, 2);

	All_NMFanswers =0;
	for z = 1:length(Trials)

				t1 = time_ns();
				#InInfo = readdlm("Input.txt");		


				#globalIter = InInfo[1,2];
				#nmfIter = InInfo[2,2];
				numberOfProcesses = Trials[z];

				# reading input file and initilizing arrays
				
				numberOfPoints   = size(inputMatrix, 1);
				numberOfSamples = size(inputMatrix, 2);

				allProcesses = zeros( numberOfPoints, numberOfProcesses, globalIter );
				allMixtures  = zeros( numberOfProcesses, numberOfSamples, globalIter );
				allDelays = zeros(size(allMixtures)); 
				allCost = Array(Float64, globalIter);
				Vars = Array(Float64,globalIter);
				Steps = Array(Float64,globalIter);

			
				opts=Dict();
				opts = Dict("runit"=>0, "convcrit"=>1e-8, "auto_corr"=>1, "maxiter"=>nmfIter, "dispiter"=>0);

				# DON'T FORGET the input matrix needs to be transposed as well as the output W and H matricies 
				All_NMFanswers = Parallel_ShiftNMF2(globalIter, nmfIter, inputMatrix' ,numberOfProcesses, opts);
			
			# //////////////////////////////////////////////////////////	
				
				
				if isdir("./Trials") == false
					mkdir("./Trials");
				end

				f=open("./Trials/TrialX$numberOfProcesses.json", "w");
				JSON.print(f,All_NMFanswers);
				close(f)

				t2 = time_ns();
				elapsedTime[z,1] = "Elapsed time for $z sources in seconds: "
				elapsedTime[z,2] = (t2 - t1)/1.0e9;
				elapsedTime[z,3] = "The average iterations per NMF for $z are: ";
				elapsedTime[z,4] = mean(Steps);

	end

	writedlm("Output.txt", elapsedTime);
	ShiftNMFk_CosCluster(X, maxSource);
	return All_NMFanswers;

end


