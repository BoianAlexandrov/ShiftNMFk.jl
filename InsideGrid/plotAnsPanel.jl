using Gadfly
#using Cairo

Htrue = readcsv("./InputInGrid/Htrue.csv");
Ttrue = readcsv("./InputInGrid/Ttrue.csv");
Wtrue = readcsv("./InputInGrid/Wtrue.csv");

x=1:size(Htrue,2)

Ha = readcsv("./Results/H3.csv");
for i=1:size(Htrue,1)
    sumH=sum(Htrue[i,:]);
    Htrue[i,:] = Htrue[i,:]./sumH;
    Wtrue[:,i] = Wtrue[:,i].*sumH;
end

Ta = readcsv("./Results/T3.csv");

Wa = readcsv("./Results/W3.csv");

# ////////////// Creating the error bars for Tfound

sd = 1*rand(size(Ta));
#sd = readcsv("./Results/Tstd3.csv");
sdmin=Ta-sd;
sdmax = Ta+sd;

# normalization of Wtrue and Wfound so that the contributions of each row add to 100%

for i = 1:size(Wa,1)
    Wa[i,:] = (1/sum(Wa[i,:]))*Wa[i,:];
    Wtrue[i,:] = (1/sum(Wtrue[i,:]))*Wtrue[i,:];
end



WaWide = Array(Float64, size(Wa, 1), 6, size(Wa, 2));
WtrueWide = Array(Float64, size(Wa, 1), 6, size(Wa, 2));
for i = 1:size(Ta,2)
  WaWide[:,:,i] = [Wa[:,i] Wa[:,i] Wa[:,i] Wa[:,i] Wa[:,i] Wa[:,i]];
  WtrueWide[:,:,i] = [Wtrue[:,i] Wtrue[:,i] Wtrue[:,i] Wtrue[:,i] Wtrue[:,i] Wtrue[:,i]];
end

ns = size(Ttrue,2);
tmean = mean(Ttrue,1);
Ttrue = Ttrue - repmat(tmean, size(Ttrue,1), 1);
#Ta = Ta + repmat([tmean[2] tmean[3] tmean[1]], size(Ttrue,1), 1);
shiftby = ceil(tmean);
shiftby =[8 8 7];
shiftby =[8 7 8];
Hnew = Array(Any,ns);
for i=1:ns
  Hnew[i] = [zeros(1,shiftby[i]) Ha[i,:]];

end

# //////////////////////// Results Plot /////////////////////

Wtp = vstack(hstack( spy(WtrueWide[:,:,1], Guide.title(""), Guide.xlabel(""), Guide.xticks(ticks=nothing),Theme(minor_label_font_size = 22pt), Scale.color_continuous(minvalue=0.0, maxvalue=0.7 )),
           spy(WaWide[:,:,1], Guide.title(""), Guide.xlabel(""), Guide.xticks(ticks=nothing),Theme(minor_label_font_size = 22pt), Scale.color_continuous(minvalue=0.0, maxvalue=0.7 )))
  , hstack( spy(WtrueWide[:,:,2], Guide.title(""), Guide.xlabel(""), Guide.xticks(ticks=nothing),Theme(minor_label_font_size = 22pt), Scale.color_continuous(minvalue=0.0, maxvalue=0.7 )),
           spy(WaWide[:,:,3], Guide.title(""), Guide.xlabel(""), Guide.xticks(ticks=nothing),Theme(minor_label_font_size = 22pt), Scale.color_continuous(minvalue=0.0, maxvalue=0.7 )))
  , hstack( spy(WtrueWide[:,:,3], Guide.title(""),  Guide.xlabel(""), Guide.xticks(ticks=nothing),Theme(minor_label_font_size = 22pt), Scale.color_continuous(minvalue=0.0, maxvalue=0.7 )),
           spy(WaWide[:,:,2], Guide.title(""), Guide.xlabel(""), Guide.xticks(ticks=nothing),Theme(minor_label_font_size = 22pt), Scale.color_continuous(minvalue=0.0, maxvalue=0.7)))
  );

#xdis = [1,2,3,4,5,6,7,8,9];
xdis = 1:size(Ttrue,1);

#=Ttp = vstack( hstack( plot(x=xdis, y=Ttrue[:,1], Geom.bar, Scale.x_discrete, Guide.xlabel("Delay to each Observer"), Guide.ylabel("Shift"), Guide.title("Ttrue"))
      , plot(x=xdis, y=Ta[:,3], Geom.bar, Scale.x_discrete, Guide.xlabel("Delay to each Observer"), Guide.ylabel("Shift"), Guide.title("Tfound")))
, hstack( plot(x=xdis, y=Ttrue[:,2], Geom.bar, Scale.x_discrete, Guide.xlabel("Delay to each Observer"), Guide.ylabel("Shift"), Guide.title("Ttrue"))
      , plot(x=xdis, y=Ta[:,2], Geom.bar, Scale.x_discrete, Guide.xlabel("Delay to each Observer"), Guide.ylabel("Shift"), Guide.title("Tfound")))
, hstack( plot(x=xdis, y=Ttrue[:,3], Geom.bar, Scale.x_discrete, Guide.xlabel("Delay to each Observer"), Guide.ylabel("Shift"), Guide.title("Ttrue"))
      , plot(x=xdis, y=Ta[:,1], Geom.bar, Scale.x_discrete, Guide.xlabel("Delay to each Observer"), Guide.ylabel("Shift"), Guide.title("Tfound"))));
=#  

Ttp = vstack( plot( 
      layer(x=xdis, y=Ta[:,1], #=ymin=sdmin[:,1], ymax=sdmax[:,1], Geom.point, Geom.errorbar=# Geom.bar, Theme(default_color=color("red"), bar_spacing=6mm))
      ,layer(x=xdis, y=Ttrue[:,1], Geom.bar)
      , Scale.x_discrete, Guide.xlabel(""), Guide.ylabel("Shift"), Guide.title("Panel D: Delays")
      , Guide.manual_color_key("Legend", ["True Value", "Found Value"], ["blue", "red"])
      ,Theme(minor_label_font_size = 22pt), Guide.yticks(ticks = [-10, 0, 10]), Guide.xticks(ticks = [1, 4, 8, 12, 16]))
,  plot(
      layer(x=xdis, y=Ta[:,3], #=ymin=sdmin[:,2], ymax=sdmax[:,2], Geom.point, Geom.errorbar,=# Geom.bar, Theme(default_color=color("red"), bar_spacing=6mm))
      ,layer(x=xdis, y=Ttrue[:,2], Geom.bar)
      , Scale.x_discrete, Guide.xlabel(""), Guide.ylabel("Shift")
      , Guide.manual_color_key("Legend", ["True Value", "Found Value"], ["blue", "red"])
      ,Theme(minor_label_font_size = 22pt), Guide.yticks(ticks = [-10, 0, 10]), Guide.xticks(ticks = [1, 4, 8, 12, 16]))
, plot(
      layer(x=xdis, y=Ta[:,2], #=ymin=sdmin[:,3], ymax=sdmax[:,3], Geom.point, Geom.errorbar,=# Geom.bar, Theme(default_color=color("red"), bar_spacing=6mm)) 
      , layer(x=xdis, y=Ttrue[:,3], Geom.bar)
      , Scale.x_discrete, Guide.xlabel(""), Guide.ylabel("Shift")
      , Guide.manual_color_key("Legend", ["True Value", "Found Value"], ["blue", "red"])
      ,Theme(minor_label_font_size = 22pt), Guide.yticks(ticks = [-10, 0, 10]), Guide.xticks(ticks = [1, 4, 8, 12, 16]))
);
      



Htp = vstack(plot( layer(x=x, y=Htrue[1,:], Geom.line, Theme(default_color=colorant"green"))
    ,layer(x=1:length(Hnew[1]), y=Hnew[1], Geom.point, Theme(default_color=color("blue")))
    , Guide.xlabel(""), Guide.ylabel("Mag"), Guide.title("Panel C: Signals"), Guide.manual_color_key("Legend", ["True Value", "Found Value"], ["green", "blue",])
    ,Theme(minor_label_font_size = 20pt))
, plot( layer(x=x, y=Htrue[2,:], Geom.line, Theme(default_color=color("green")))
    ,layer(x=1:length(Hnew[3]), y=Hnew[3], Geom.point, Theme(default_color=color("blue")))
    , Guide.xlabel(""), Guide.ylabel("Mag"),  Guide.manual_color_key("Legend", ["True Value", "Found Value"], ["green", "blue",])
    ,Theme(minor_label_font_size = 20pt))
, plot( layer(x=x, y=Htrue[3,:], Geom.line, Theme(default_color=color("green")))
    ,layer(x=1:length(Hnew[2]), y=Hnew[2], Geom.point, Theme(default_color=colorant"blue"))
    , Guide.xlabel(""), Guide.ylabel("Mag"),  Guide.manual_color_key("Legend", ["True Value", "Found Value"], ["green", "blue",])
    ,Theme(minor_label_font_size = 20pt))
);


# ////////////////// Norm and Silhouette Plot

X = readcsv("./InputInGrid/Observation.csv");

Trials = [1, 2, 3, 4, 5];

Norm = Array(Float64, size(Trials));
Silhouette = Array(Float64, size(Trials));

for z = 1:length(Trials)

  numberOfProcesses = Trials[z];
  #Rec = readcsv("./Results/Reconstruction$numberOfProcesses.csv");
  Rec = readcsv("./Results/Reconstruction$numberOfProcesses.csv");
  #normalize X & Rec
  for i=1:size(X,1)
      sumX=sum(X[i,:]);
      sumRec=sum(Rec[i,:]);
      X[i,:] = X[i,:]./sumX;
      Rec[i,:] = Rec[i,:]./sumRec;
    end

  cost=0.5*vecnorm(X-Rec,2)^2;

  Norm[z] = cost;

  #avgStability = readcsv("./Results/avgStability$numberOfProcesses.csv");
  avgStability = readcsv("./Results/avgStability$numberOfProcesses.csv");
  Silhouette[z] = mean(avgStability);
  println(Norm[z]);
  #println(Silhouette[z]);

end

x = Trials;

ResultPlot = plot(
  layer(x=x, y=Norm, Geom.point, Geom.line, Theme(default_color=color("green"))),
  layer(x=x, y=Silhouette, Geom.point,  Geom.line, Theme(default_color=color("red"))),
  Guide.xlabel("Number of Sources"), Guide.ylabel("Fro. Norm and Silhouette Value"), Guide.title("Results"), Guide.manual_color_key("Legend", 
    ["Silhouette Value", "Norm Minimization"], ["red", "green",])
);


#//////////////////////// Putting together plots in panels //////////////////


Left = vstack(ResultPlot, Htp);
Right = vstack(Wtp, Ttp);

Allt = hstack(Left ,Right);
#Allt = hstack(Wtp, Htp, Ttp);

#draw(SVG("ResCompare2.svg", 50cm, 25cm), Allt);

draw(SVG("ResComparePanel.svg", 50cm, 50cm), Allt);


