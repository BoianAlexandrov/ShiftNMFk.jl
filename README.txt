The module ShiftNMFk is contained within the ShiftNMFkMD.jl file in src.
It uses the rest of the files in src to export several functions used to perform the entire ShiftNMFk procedure and find the locations of the sources. 

Two examples are provided. Both feature a grid of 16 detectors and 3 sources. In one example the sources are inside the grid while in the other they are outside. These examples are the same found in our ShiftNMFk paper.  

To run and example:

First run FirstLoad.jl which adds the required packages. (In the future this will be done by the REQUIRE file).

Start julia with the number of extra cores you desire IN the folder of the example : julia -p X

Then simply run one of the Test.jl files in OutsideGrid and InsideGrid.
The Test.jl files are commented in a way that explains the ShifNMFk functions and their inputs. 
Initially only a small number of runs and iterations are performed to see if the code runs alright. These variables can be changed within the Test.jl files for running the substantial runs needed to reproduce the examples in our paper. 
Running the example creates several folders with results and plots the Silhouette and Norm plot as well as the locations of the sources.

Run plotAnsPanel.jl for a graph that compares the ShiftNMFk results with the original signals used to create the example.  

Unfortunately my functions still through up warnings regarding deprecated syntax. I will fix this soon, I just wanted to get the module up an running first. 


