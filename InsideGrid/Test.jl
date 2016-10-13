@everywhere function setdir(dir)	
	if isdir(dir)
		cd(dir)
	end
end

@everywhere function setdir()
	dir = remotecall_fetch(1, ()->pwd())
	setdir(dir)
end


@everywhere setdir()

@everywhere include("../src/ShiftNMFkMD.jl"); # including the module on all cores
@everywhere using ShiftNMFk    				# using the module ShiftNMFk on all cores

X = readcsv("./InputInGrid/Observation.csv");		#Inputing the observation matrix of the desired example
micPos = readcsv("./InputInGrid/MicPosition.csv"); # The coordinates of the detectors in the grid
maxSource = 5;										#5 max number of sources
globalIter = 	10;								# 1000 NMF runs for each guess of a sources								
nmfIter = 3000;									# 80,000 max number of iterations for each source.
locIter = 10;										# 1000 minimizations are performed to find the location
numT = 180; 										#The numbes of sample points that make up our signals (time samples) 
nd = 16;											# The number of detectors in our grid



shiftNMFk(X, maxSource, globalIter, nmfIter);					 

Sil, Norm = Plot(X, maxSource);						#Plots the Norm and Silhouette Value graph and returns both 

aic_min, nopt = AIC_final( Norm, Sil, numT, nd)		# nopt returs the correct number of sources in the experiment. 

W,H,T,Tstd = ResultsForNumSources(nopt);			#returns the mixing(W), signal(H), and delay(T) matricies for the number
													#of sources nopt. Also returs sigma of the delays(Tstd)


SourcePositions = FindLocations(W, T, Tstd, micPos, locIter);	#Returns the coordinates of the sources. 
