#using ProfileView

#@everywhere 
function mgetopt(opts, varname, default, varargin...)
  if haskey(opts, varname)
    var = opts[varname];
  else
    var = default;
  end
  return var;
end

#@everywhere 
function estTimeAutCor(Xf, W, Hf, T, f, smoothtype)
  noc=size(W,2);
  sX=2*size(Hf,2)-2;
  Xft=Array(Complex{Float64}, size(Xf));                   
  Xft[:]=Xf[:];
  t1=randperm(size(Xf,1));                                           
  t2=randperm(noc);

  
  for k=t1
    for d=t2
        nocm=setdiff(1:noc,d);

          Xft[k,:]=Xf[[k],:]-W[[k],nocm]*(Hf[nocm,:].*exp(T[[k],nocm]'*f[[1],:]));
          C=(conj(Xft[[k],:]).*Hf[[d],:]);
          C=[C conj(C[:,end-1:-1:2])];
          C=real(ifft(C,2));

          y=maximum(C);
          ind = 1;
          ind = find(C.==y);
          ind = ind[1];

          T[k,d] =(ind-sX)-1;
          if smoothtype =="none"
            W[k,d]=C[ind]/((sum(2*(Hf[[d],:].*conj(Hf[[d],:])))-sum(Hf[[d],[1,end]].*conj(Hf[[d],[1,end]])))/sX);
          end
          if abs(T[k,d])>(sX/2)
              if T[k,d]>0
                  T[k,d]=T[k,d]-sX;
              else
                   T[k,d]=T[k,d]+sX;
              end
          end
    end
  end
  #Wc = real(Wc);

  return T, W;
end

#@everywhere 
function update_TmH(T,P,Recfd,Hall,Told)
  cost=100;

  nyT = P["nyT"];
  Hf = P["Hf"];
  W = P["W"];
  sizeX2 = P["sizeX2"];
  Xf=P["Xf"];
  f=P["f"];
  w=P["w"];


  Rep1 = Array(Float64,size(W,1),length(f));  
  Rep2 = Array(Complex{Float64}, size(W,1), size(Hf,2)); 

  for d=1:size(W,2)
    Exp = exp(T[:,d]*f);                                                       
    
    for i=1:length(f)
      Rep1[:,i] = W[:,d];
    end
    Rep2 = Array(Complex{Float64}, size(W,1), size(Hf,2)); 
    for i=1:size(W,1)
      Rep2[i,:] = Hf[d,:];
    end
    D=(Rep1.*Exp).*Rep2;
    Recfd[:,:,d] = D;
  end

  Recf=sum(Recfd,3);

  Xnew = Xf-Recf;                                                
  RepX= Array(Complex{Float64}, size(Xnew,1), size(Xnew,2), size(Recfd,3));  
  for i=1:size(Recfd,3)
      RepX[:,:,i] = conj(Xnew[:,:]);
  end
  Q=Recfd.*RepX;                                                  

  WnF=w.*(f.^2);                                                  
  RP1= Array(Complex{Float64}, size(Q,1), length(WnF), size(Q,3));      
  for i=1:size(Q,3)
    for k=1:size(Q,1)
    RP1[k,:,i] = WnF; 
  end
  end
  RP = RP1.*real(Q);

  SM=sum(RP,2);
  Hdiag= Array(Float64,size(SM,1),size(SM,3)); 

  for i=1:size(RP,3)
    Hdiag[:,i] =2*SM[:,:,i];
  end


  for d=1:size(W,2)

    SM2 = sum(RP1.*real(Recfd.*conj(Recfd[:,:,d])),2);
    Hall[:,:,d]=2*squeeze(SM2,2);
  end

  WnF = w.*f;

  for i=1:size(Q,3)
    for k=1:size(Q,1)
    RP1[k,:,i] = WnF; 
    end
  end

  grad=squeeze(sum(RP1.*(conj(Q)-Q),2),2);

  for i=1:size(grad,1)
    ind=find(abs(grad[i,:]).>1e-6);
    grad[i,ind]=grad[[i],ind]/(-diagm(conj(Hdiag[[i],ind])'[:,1])-squeeze(Hall[[i],ind,ind],1));
  end

  ind1=find(w.==2); # Areas used twice
  ind2=find(w.==1); # Areas used once
  cost_old=vecnorm(Xf[:,ind1]-Recf[:,ind1],2)^2;
  cost_old=cost_old+0.5*vecnorm(Xf[:,ind2]-Recf[:,ind2],2)^2;

  keepgoing=1;
  Told[:]=T[:];

  kill = 0;
  while keepgoing == 1
        T=Told-nyT*grad;
        for d=1:size(W,2)

            Recfd[:,:,d]=(repmat(W[:,d],1 ,length(f)).*exp(T[:,d]*f)).*repmat(Hf[[d],:],size(W,1),1);

        end
        Recf=sum(Recfd,3);
        
        cost=vecnorm(Xf[:,ind1]-Recf[:,ind1],2)^2;
        cost=cost+0.5*vecnorm(Xf[:,ind2]-Recf[:,ind2],2)^2;

        if cost<=cost_old || kill >300
            keepgoing=0;
            nyT=nyT*1.2;
        else
            keepgoing=1;
            nyT=nyT/2;
            kill=kill+1;
        end
  end
   T = mod(real(T),sizeX2);
   ind=find(T.>floor(sizeX2/2));
    T[ind]=T[ind]-sizeX2;


  return T, nyT, cost;
end

#@everywhere 
function ShiftNMF(X, noc, opts, varargin...)
  # can get the return values by calling:   W, H, T, varexpl, cost = ShiftNMF(X, noc, argin)


  #opts = Dict();                        #opts is a dictionary in julia instead of a struct in MATLAB
  if length(varargin)>0
    for  i = 1:2:length(varargin)-1
      #d1 = [varargin[i] => varargin[i+1]]; depricatd syntax
      d1 = Dict(varargin[i] => varargin[i+1]);
      opts = merge(opts, d1);
    end
  end
  varexplold=0;
  cost=1;
  runit = mgetopt(opts,"runit",0);

  if !haskey(opts,"H") && runit != 1
    println("Finding the best out of", 10, "initial solutions");
    for k=1:10
      if k==1
        println("Now Estimating ", k, " st solution");
      elseif k==2
        println("Now Estimating ", k, " nd solution");
      elseif k==3
        println("Now Estimating ", k, " rd solution");
      else
        println("Now Estimating ", k, " th solution");
      end

      optsn = opts;
      
      optsn = merge(optsn,Dict("dispiter"=>0));
      optsn = merge(optsn, Dict("runit"=>1));
      optsn = merge(optsn, Dict("maxiter"=>25));

     Wq, Hq, Tq, varexpl, cost= ShiftNMF(X,noc,optsn);
      #cost=Int64;

      if varexpl>varexplold
        println("Best variation explained ", varexpl);
        
        varexplold=varexpl;
        W=Wq;
        H=Hq;
        T=Tq;
      end
    end
  else
    mx = maximum(abs(X[:]));
    W=mgetopt(opts,"W",mx*rand(size(X,1),noc));
    H=mgetopt(opts,"H",mx*rand(noc,size(X,2)));
    T=mgetopt(opts,"T",zeros(size(W)));
    
  end


  SumofRowX = Array(Float64,size(X,1));

  # New normalization of X matrix
  for i=1:size(X,1)
    sumX=sum(X[i,:]);
    X[i,:] = X[i,:]./sumX;
    SumofRowX[i]=sumX;
  end

	#Normalizing W from original ShiftNMF

  #W=W./repmat(sum(W,2),1,size(W,2));
  W=W./repmat(sqrt(sum(W.^2,1)),size(W,1),1);

  nyT=mgetopt(opts,"nyT",1);
  maxiter=mgetopt(opts,"maxiter",1000);
  conv_crit=mgetopt(opts,"convcrit",1e-6);
  SST=vecnorm(X,2)^2;                        # This is the Frobenius norm, sqrt(sum(diag(X'*X))) of the oservation matrix X.
  constH=mgetopt(opts,"constH",0);
  constW=mgetopt(opts,"constW",0);
  constT=mgetopt(opts,"constT",0);
  lambda=mgetopt(opts,"lambda",0);
  auto_corr=mgetopt(opts,"auto_corr",1);
  smoothtype= mgetopt(opts,"smoothtype","none");
  smoothnorm= mgetopt(opts,"smoothnorm","L2");
  alpha=mgetopt(opts,"alpha",1);
  dispiter=mgetopt(opts,"dispiter",1);


  if smoothtype == "curv"          #Smoothtype
    L=-2*eye(size(H,2));
    L[2:end,1:end-1]=L[2:end,1:end-1]+eye(size(H,2)-1);
    L[1:end-1,2:end]=L[1:end-1,2:end]+eye(size(H,2)-1);
    L[:,1]=0;
    L[:,end]=0;
  elseif smoothtype == "grad"
    L=-eye(size(H,2));
    L[2:end,1:end-1]=L[2:end,1:end-1]+eye(size(H,2)-1);
    L[:,end]=0;
  elseif smoothtype == "regularization"
    L=eye(size(H,2));
  else
    L = 0;
  end

  if L != 0
    L=sparse(L);                  #Done to save memory, may not be needed in julia
  end
  N=size(X,2);
  #f=im*2*pi*[0:N-1]'/N;
  f=im*2*pi*collect(0:N-1)'/N;
  Xf=fft(X,2);                         # Applies fft across 2nd dimension (columns)
  #Xf=Xf[:,1:floor(size(Xf,2)/2)+1];    # depricated syntax
  Xf=Xf[:,1:Int(floor(size(Xf,2)/2))+1];       #Setting Xf and Hf to half of the fft() (columns) because fft repeats (symetric)
  Hf=fft(H,2);
  #Hf=Hf[:,1:floor(size(Hf,2)/2)+1];    # depricated syntax
  Hf=Hf[:,1:Int(floor(size(Hf,2)/2))+1];
  f=f[1:size(Xf,2)]'*(-1);


  # Initial Cost
  Rec=Array(Float64,size(X));                       

  for i=1:size(W,1)
     Hft = Hf.*exp(T[[i],:]'*f);
     if mod(size(X,2),2)==0
         Hft=[Hft conj(Hft[:,end-1:-1:2])];
     else
         Hft=[Hft conj(Hft[:,end:-1:2])];
     end
     Ht=real(ifft(Hft,2));
     Rec[i,:]=W[[i],:]*Ht;
  end
  cost=0.5*vecnorm(X-Rec,2)^2;          # This is the Frobenius norm, sqrt(sum(diag(X'*X))) of the oservation matrix X.
  varexpl=(SST-2*cost)/SST;             # The Fro norm of the observations X is SST.

  if !isempty(smoothtype)
    HL=H*L;
    if smoothnorm == "L1"
      smoothcost = lambda*sum(sum(abs(HL)));
    elseif !isempty(smoothtype)
      smoothcost = lambda*0.5*vecnorm(HL,2)^2;
    end
  end

  cost_oldt=cost+smoothcost;
  dcost=Inf16;
  told=Int(time_ns())/1e9;
  iter=0;

  #Display Algorithm
  if dispiter == 1
    println(" ");
    println("Shifted Non-negative matrix factorization");
    println("A ", noc, " component model will be fitted");
    println("To stop algorithm press control C");
    println("Smoothness by ", smoothtype, "using ", smoothnorm, "-norm, imposed with strenght ", lambda);
    dheader = "Iteration,    Expl.var.,    Cost func.    ,Delta costf.,   Time(s)   ,    H-stepsize    ,   Plato";
    dline = "-------------+--------------+--------------+--------------+--------------+--------------+--------------+";
    println(dline);
    println(dheader);
    println(dline);
  end

  plato = 0;

  #costiter = Float64[];
  #varexpl = [varexpl];

    # Declarations for SPEED `
  gradnH=Array(Complex{Float64}, size(Hf));
  gradpH=Array(Complex{Float64}, size(Hf));

  gradconp=Array(Float64,size(H));                            
  gradconn=Array(Float64,size(H));                            

  Hold=Array(Float64,size(H));
  Recf=Array(Complex{Float64}, size(Xf));

  gradnW = Array(Float64,size(W));                                 
  gradpW = Array(Float64,size(W));

  P = Dict();

  small_w = ones(1,length(f))
  if mod(size(X,2),2)==0
    small_w[2:end-1]=2;               
  else
    small_w[2:end]=2;
  end    

  #P={"Xf"=>Xf, "sizeX2"=>size(X,2), "w"=>small_w, "f"=>f};  depricated syntax
  P=Dict("Xf"=>Xf, "sizeX2"=>size(X,2), "w"=>small_w, "f"=>f);
  
  Recfd=Array(Complex{Float64},size(W,1),size(Hf,2),size(W,2));
  
  Hall=Array(Float64,size(W,1),size(W,2),size(W,2));
  
  Told=Array(Float64,size(T));       

  while iter < maxiter  && (dcost >= cost*conv_crit || mod(iter,20) == 0)

    iter = iter+1;
    if dcost <= cost*conv_crit
      plato = plato+1;
    else
      plato = 0;
    end

    if mod(iter,10)==0 && dispiter == 1
      println(dline);
      println(dheader);
      println(dline);
    end

    #Update H
    if constH == 0
         
          for i=1:size(Hf,2)
              Wf=(W.*exp(T*f[i]));
              gradnH[:,i]=Wf'*Xf[:,i];
              gradpH[:,i]=Wf'*Wf*Hf[:,i];
          end

          if mod(size(X,2),2)==0
              gradnH2=[gradnH conj(gradnH[:,end-1:-1:2])];
              gradnHt=real(ifft(gradnH2,2));
              gradpH2=[gradpH conj(gradpH[:,end-1:-1:2])];
              gradpHt=real(ifft(gradpH2,2));
          else
                  gradnH2=[gradnH conj(gradnH[:,end:-1:2])];
                  gradnHt=real(ifft(gradnH2,2));
                  gradpH2=[gradpH conj(gradpH[:,end:-1:2])];
                  gradpHt=real(ifft(gradpH2,2));
          end
          gradnHt[gradnHt.<0]=0;
          gradpHt[gradpHt.<0]=0;

          if smoothtype != "none"
                if smoothnorm =="L1"
                    gradconstr=lambda*sign(HL).*repmat(sum(L,2)',size(H,1),1);
                else
                    gradconstr=lambda*HL*L';
                end

                indp = find(x->(x > 0),gradconstr);
                indp = find(x->(x < 0),gradconstr);
                
                gradconp[indp]=gradconstr[indp];
                gradconn[indn]=-gradconstr[indn];
                grad=(gradnHt+gradconn)./(gradpHt+gradconp+eps());
          else
                grad=(gradnHt)./(gradpHt+eps());
            end
          
          Hold[:]=H[:];
          
          keepgoing=1;
          kill =0;
          while keepgoing == 1
                 H=Hold.*(grad.^alpha);
                 Hf=fft(H,2);
                 Hf=Hf[:,1:Int(floor(size(Hf,2)/2))+1];
                 for i=1:size(W,1)
                    Hft=Hf.*exp(T[[i],:]'*f);
                    Recf[i,:]=W[[i],:]*Hft;
                 end
                 if mod(size(X,2),2)==0
                      cost=1/size(H,2)*(0.5*vecnorm(Xf[:,[1,end]]-Recf[:,[1,end]],2)^2+vecnorm(Xf[:,2:end-1]-Recf[:,2:end-1],2)^2);
                 else
                      cost=1/size(H,2)*(0.5*vecnorm(Xf[:,1]-Recf[:,1],2)^2+vecnorm(Xf[:,2:end]-Recf[:,2:end],2)^2);
                 end
                 if ~isempty(smoothtype)
                     HL=H*L;
                     if smoothnorm == "L1"
                        smoothcost=lambda*sum(sum(abs(HL)));
                     else
                        smoothcost=lambda*0.5*vecnorm(HL,2)^2;
                     end
                     cost=cost+smoothcost;
                 end
                 if cost>cost_oldt && kill <300
                     alpha=alpha/2;
                     keepgoing=1;
                     kill = kill+1;
                 else
                     alpha=alpha*1.2;
                     keepgoing=0;
                 end
          end
    end
    if ~isempty(smoothtype)
        cost=cost-smoothcost;
    end


    # Update W and calculate cost function
     if constW == 0 || (constH == 0 && constT == 1)

        
        for i=1:size(W,1)
           Hft=Hf.*exp(T[[i],:]'*f);
           if mod(size(X,2),2)==0
               Hft=[Hft conj(Hft[:,end-1:-1:2])];
           else
               Hft=[Hft conj(Hft[:,end:-1:2])];
           end
           Ht=real(ifft(Hft,2));
           Ht[Ht.<0]=0;
           gradnW[i,:]=X[[i],:]*Ht';
           gradpW[i,:]=W[[i],:]*(Ht*Ht');
        end
        if constW == 0 && smoothtype != "none"
            tx = sum(gradpW.*W,1);
            ty = sum(gradnW.*W,1);
            gradnW = gradnW + repmat(tx,size(W,1),1).*W;
            gradpW = gradpW + repmat(ty,size(W,1),1).*W;
            W=W.*(gradnW)./(gradpW+eps());
            W=W./repmat(sqrt(sum(W.^2,1)),size(W,1),1);
        else
            W=W.*(gradnW)./(gradpW+eps());
        end
        if  (constH ==0 && constT == 1)              # Calculate Reconstruction if cost function is needed
           for i=1:size(W,1)
               Hft=Hf.*exp(T[[i],:]'*f);
               if mod(size(X,2),2)==0
                   Hft=[Hft conj(Hft[:,end-1:-1:2])];
               else
                    Hft=[Hft conj(Hft[:,end:-1:2])];
               end
               Ht=real(ifft(Hft,2));
               Rec[i,:]=W[[i],:]*Ht;
           end
        end
      end

    #Update T and calculate cost

      if constT != 1

          if mod(iter,20)==0 && iter<maxiter-20 && auto_corr == 1

              println("Reestimating shifts by cross-correlation");
              T, W= estTimeAutCor(Xf,W,Hf,T,f,smoothtype);                         #Auto Correlation Method
          end

          P["W"] = W;
          P["Hf"] = Hf;
          P["nyT"] = nyT;

          (T,nyT,cost)=update_TmH(T,P,Recfd,Hall,Told);
          cost=cost/size(H,2);                                # Use Parseval identity for cost function

          sse=2*cost;
      elseif constW == 0
          sse=vecnorm(X-Rec,2)^2;
          cost=0.5*sse;
      else
          sse=2*cost;
      end
      if !isempty(smoothtype)
          cost=cost+smoothcost;
      end

      #Display output of iterations
      dcost=cost_oldt-cost;
      cost_oldt = cost;
      #push!(costiter,cost);
      varexpl=((SST-sse)/SST);

      if rem(iter,1)==0 && dispiter == 1
          t=Int(time_ns())/1e9;
          tim=t-told;
          told=t;
          println("  ",iter,"    " ,varexpl," ",cost, "   ",dcost,"        ",tim,"      ",alpha ,"     ", nyT, "    ", plato);
      end
  end

  if(iter >= maxiter)
      println("Stopped because max # of iterations reached");
  elseif(dcost<cost*conv_crit)
      println("Stopped because dcost is smaller than cost*conc_criteria");
      println(dcost);
  else
      println("I dont know why it stopped");
  end


  #cost=costiter;

  #Align H and T
  tmean=mean(T,1);
  T=T-repmat(tmean,size(T,1),1);
  Hf=Hf.*exp(tmean'*f);
  if mod(size(X,2),2)==0
      Hf=[Hf conj(Hf[:,end-1:-1:2])];
  else
      Hf=[Hf conj(Hf[:,end:-1:2])];
  end
  H=real(ifft(Hf,2));
  H[H.<0]=0;


  # Changing W to reverse X normalization in the beginnning
  #=for i=1:size(W,1)
    W[i,:]= W[i,:].*SumofRowX[i];
  end=#

  #///////////////////  Normalize H and W  ///////////////////////////
   
    for i=1:size(H,1)
        sumH=sum(H[i,:]);
        H[i,:] = H[i,:]./sumH;
        W[:,i] = W[:,i].*sumH;
    end


  
  return W ,H, T, varexpl, cost, iter;
end



