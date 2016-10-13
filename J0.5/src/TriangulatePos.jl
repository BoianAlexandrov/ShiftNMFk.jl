

#@everywhere 
function myfunc(x::Vector, grad::Vector, T, Tstd, micPos)
    if length(grad) > 0
        
    end

    s=0;
  	k = 1;
    for i = 2:2:2*size(T,2)+1

    		r1 = sqrt( (x[i] - micPos[1,1])^2 + (x[i+1] - micPos[1,2])^2 );

    		for h = 2:size(T,1)
    			r2 = sqrt( (x[i] - micPos[h,1])^2 + (x[i+1] - micPos[h,2])^2 );

    			s = s + ((T[1,k] - T[h,k])/sqrt(Tstd[1,k]^2 + Tstd[h,k]^2) - (( (r1 - r2)/x[1] )/sqrt(Tstd[1,k]^2 + Tstd[h,k]^2)))^2

    		end

    	
    	k=k+1;	 
    end
    return s
end


#@everywhere 
function Triangulate(T, Tstd,  W, micPos, xtol) 

	ns = size(T,2);
	
	positions =  zeros(size(T,2) ,2);

	# //// Generating initial conditions range. /////////////////////

	Rmax = Array(Float64, ns);
	searchSquare = Array(Float64, ns, 4);
	farMics = Array(Float64, ns, 2);
	for i = 1:ns
		idClose = findin(T[:,i], minimum(T[:,i]));
		idFar = findin(T[:,i], maximum(T[:,i]));
		Tcf = sqrt( (micPos[idClose, 1]-micPos[idFar, 1]).^2 + (micPos[idClose, 2]-micPos[idFar, 2]).^2)
		Wc = W[idClose,i];
		Wf = W[idFar,i];
		R = Tcf./((1-Wf/Wc));		# This is the only thing I have changed so far for the ne chi-minimization, I added the 2* denomonator
		Rmax[i] = abs(R[1]);

		searchSquare[i,:] = [micPos[idFar, 1]+Rmax[i] micPos[idFar, 1]-Rmax[i] micPos[idFar, 2]+Rmax[i] micPos[idFar, 2]-Rmax[i]];
		farMics[i,:] = [micPos[idFar, 1] micPos[idFar, 2]];

	end

	lb = [0.00001];
	ub = [Inf];
	init = [10*rand()];
	for i = 1:ns
		lb = [lb; [searchSquare[i,2] ; searchSquare[i,4] ] ];
		ub = [ub; [searchSquare[i,1] ; searchSquare[i,3] ] ];
		init =[init ; farMics[i,1]+rand([-1,1])*Rmax[i]*rand() ; farMics[i,2]+rand([-1,1])*Rmax[i]*rand()];
	end


	#//////////////////////////////
	

				#init = 20*rand(1+2*ns)
				opt = Opt(:LN_COBYLA, 2*ns+1)
				lower_bounds!(opt, lb);
				upper_bounds!(opt, ub);
				xtol_rel!(opt,xtol)
				maxtime!(opt, 120)

				min_objective!(opt, (x,g) -> myfunc(x, g, T, Tstd, micPos))

				(minf,minx,ret) = optimize(opt, init)
				#println("got $minf at $minx after $count iterations (returned $ret)")

				k=1
				for i = 2:2:2*size(T,2)+1
					positions[k,1] = minx[i];
					positions[k,2] = minx[i+1];
					
					k=k+1;
				end
				speed = minx[1]

				return positions, speed, init, minf;
end



		