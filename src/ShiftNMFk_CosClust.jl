
function ResultsForNumSources(N);

	N = N[1];
	if isdir("./Results") == false
			println("ERROR: Results directory does not exist!")
			throw(ArgumentError("ERROR: Results directory does not exist!"))
	else
		W = readcsv("./Results/W$N.csv");
		H = readcsv("./Results/H$N.csv");
		T = readcsv("./Results/T$N.csv");
		Tstd = readcsv("./Results/Tstd$N.csv");
	end
	return W, H, T, Tstd;
end

function ShiftNMFk_CosCluster(X, maxSource)

	#Trials =  [1, 2, 3, 4, 5];

	Trials = collect(1:maxSource);

	inputMatrix = X;
	inputMatrix=inputMatrix';
	numberOfPoints = size(inputMatrix, 1);
	numberOfSamples = size(inputMatrix, 2);

	VarChange = cell(length(Trials), 6);
	ALLCOST = cell(length(Trials));


	for z = Trials

				

				# DON'T FORGET the input matrix needs to be transposed as well as the output W and H matricies 
				W, H, T, cost, Vars = ParseTrials(z);
				globalIter = size(H,3);
				nmfIter = 100000;
				numberOfProcesses = size(H,1);

				numberOfPoints   = size(inputMatrix, 1);
				numberOfSamples = size(inputMatrix, 2);

				allProcesses = zeros( numberOfPoints, numberOfProcesses, globalIter );
				allMixtures  = zeros( numberOfProcesses, numberOfSamples, globalIter );
				allDelays = zeros(size(allMixtures));
				allProcessesOG = zeros( numberOfPoints, numberOfProcesses, globalIter );
				allMixturesOG  = zeros( numberOfProcesses, numberOfSamples, globalIter );
				allDelaysOG = zeros(size(allMixtures));
				allCostOG = Array(Float64, globalIter);




				for i =1:globalIter
					allProcessesOG[:, :, i] = H[:,:,i]';					# Again transposed
					allMixturesOG[:, :, i] = W[:,:,i]';
					allDelaysOG[:,:,i] = T[:,:,i]';
					allCostOG[i] = cost[i];
				end


				# Getting Rid of the worst solutions
				
				GoodTrials = ones(globalIter);
				maxVar = 0.08

				for f = 1:10

					for i=1:globalIter
						
						if Vars[i] > maxVar

							GoodTrials[i] = 0;
						end
						if minimum(W[:,:,i]) < .001
							GoodTrials[i] = 0;
						end
					end

					
					if length(findin(GoodTrials,1)) > 10
						break;
					else
						maxVar = maxVar + 0.1;
						GoodTrials = ones(globalIter);
					end
					println("The Norm cut off for $z is now $maxVar");
				end

				sizeGood = size(find(GoodTrials.==1));
				println("Good Trials after norm cut: $sizeGood");

				if length(findin(GoodTrials,1)) < 1    # So that the script doesnt stop if GoodInd=0
					GoodTrials = ones(globalIter);
					
				end

				sizeGood = size(find(GoodTrials.==1));
				println("Good Trials after regen cus goodtrials  = 0: $sizeGood");



				GoodInd = find(GoodTrials.== 1);

				allProcesses = allProcessesOG[:,:,GoodInd];
				allMixtures = allMixturesOG[:,:,GoodInd]
				allDelays = allDelaysOG[:,:,GoodInd];

				allCost = allCostOG[GoodInd];


				allProcessesForK = reshape(allProcesses, size(allProcesses,1), size(allProcesses,2)*size(allProcesses,3));
				allDelaysForK = Array(Float64, size(allDelays,2),0);
				allMixturesForK = Array(Float64, size(allMixtures,2),0);


				for i = 1:size(allDelays,3)

					allDelaysForK = [allDelaysForK allDelays[:,:,i]'];
					allMixturesForK = [allMixturesForK allMixtures[:,:,i]'];
				end



				# clustering extracted processes
				
				if z !=1
					Rtoy = kmeans(allProcessesForK, z; maxiter=1000);

					#idx = reshape(assignments(Rtoy), z, Int(length(assignments(Rtoy))/z));
					idc = assignments(Rtoy);
				else
					idc = ones(Int64, 1, size(allProcesses,3));
				end

				idx = cell(z);
				idxOLD = cell(z);

				# ///////////////// Clustering T ////////////////////
				
				orderedT = cell(z);

				for d =1:z
		
						orderedT[d] = allDelaysForK[:,findin(idc, d)];
						idx[d] = findin(idc, d);
						idxOLD[d] = findin(idc, d);
			
				end

				maxSpM = 0.8;

				for jk = 1:10

					tLeft = 0;

					for f = 1:length(orderedT)

						kg = 1;
						for i=1:size(orderedT[f],2)

								Tpos = abs(orderedT[f][:,i]);
								Tmean = mean(Tpos); 
								StdTrue = std(Tpos);

							if StdTrue./Tmean >= maxSpM
								idx[f] = [idx[f][1:kg-1]; idx[f][kg+1:end]];
								kg = kg-1;
							end
							kg = kg+1;

						end

						tLeft = tLeft + length(idx[f]);
						
					end

				
					if tLeft > 10*z 
						break;
					else
						maxSpM = maxSpM + 0.1;

						idx = cell(z);
						for k = 1:length(idx)
							idx[k] = idxOLD[k][:];
						end
		
					end
					println("The SpM cut off for $z is now $maxSpM");

				end


			# //////////////////////////////////////////////////////////	


			BIGidx = [];
			for i = 1:length(idx)
				BIGidx = [BIGidx; idx[i]];
			end
			sort!(BIGidx);

			allMixturesForK = allMixturesForK[:,BIGidx];
			allDelaysForK = allDelaysForK[:,BIGidx];

			mixtures  = Array(Float64, size(allMixturesForK,1), z);
			delays  = Array(Float64, size(allDelaysForK,1), z);
			Tstd = Array(Float64, size(allDelaysForK,1), z);

			if z !=1

				Dt = allProcessesForK[:,BIGidx];

				R = kmeans( Dt, z; maxiter=1000);
								 	
				Dist = pairwise(CosineDist(), Dt);
				avgStabilityProcesses = mean(silhouettes(R, Dist));
				processes = R.centers;

				for i=1:z
					mixtures[:,i] = mean(allMixturesForK[:,findin(R.assignments, i)], 2);
					delays[:,i] = mean(allDelaysForK[:,findin(R.assignments, i)], 2);
					Tstd[:,i] = std(allDelaysForK[:,findin(R.assignments, i)], 2);
				end

			else
				avgStabilityProcesses =1;
				processes = mean(allProcessesForK[:,BIGidx],2);
				mixtures[:,1] = mean(allMixturesForK, 2);
				delays[:,1] = mean(allDelaysForK, 2);
				Tstd[:,1] = std(allDelaysForK, 2);
			end



			sizeGood = length(BIGidx)./z;
			println("Good Trials after maxSpM cut: $sizeGood");

			ALLCOST[z] = allCost[:];


				# Generate Observation Matrix by taking into account Shift

				usecomp = 1:numberOfProcesses;
				H =  processes';   
				W =  mixtures;  
				T=  delays; 
				X=Array(Float64,size(W,1),size(H,2));
						Hf=fft(H,2);                         
						Hf=Hf[:,1:Int(floor(size(Hf,2))/2)+1];      
						N=size(H,2);                            
						#f=im*2*pi*[0:N-1]'/N;  		# depricated syntax 
						f=im*2*pi*collect(0:N-1)'/N;                  
						f=f[1:size(Hf,2)]'*(-1);                      
						for i=1:size(W,1)                       
						   Hft=Hf[usecomp,:].*exp(T[i,usecomp]'*f); 
						   Hft=[Hft conj(Hft[:,end-1:-1:2])];       
						   Ht=real(ifft(Hft,2));                
						   X[i,:]=W[i,usecomp]*Ht;						
						end
				X[X.<0]=0;



				dataRecon = X'; #processes * mixtures and delays;

				dataReconCorr = zeros(numberOfSamples, 1);

				for i = 1 : numberOfSamples
					dataReconCorr[i] = cor( inputMatrix[:,i], dataRecon[:, i] );
				end


				if isdir("./Results") == false
						mkdir("./Results");
				end


				writecsv("./Results/dataReconCorr$numberOfProcesses.csv",dataReconCorr);
				writecsv("./Results/T$numberOfProcesses.csv",T);
				writecsv("./Results/W$numberOfProcesses.csv",W);
				writecsv("./Results/H$numberOfProcesses.csv",H);
				writecsv("./Results/Tstd$numberOfProcesses.csv",Tstd);
				writecsv("./Results/Reconstruction$numberOfProcesses.csv",X);
				writecsv("./Results/avgStability$numberOfProcesses.csv",avgStabilityProcesses);
				writecsv("./Results/Cost$numberOfProcesses.csv",allCost);

				VarChange[z,1] = "NMF runs for $z Sources are taken below this % of reconstruction Error: ";
				VarChange[z,2] = maxVar*100;
				VarChange[z,3] = "For $z Sources MaxSpM: ";
				VarChange[z,4] = maxSpM;
				VarChange[z,5] = "For $z Sources we have this many trials after SpM cut: ";
				VarChange[z,6] = sizeGood;

	end
			
			writedlm("ErrorLog.txt", VarChange);
			Plot(inputMatrix', maxSource);

end