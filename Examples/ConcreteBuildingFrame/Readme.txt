Contents under \ConcreteBuildingFrame:

"ConcreteBuildingFrame.m" 
* This is the main function to perform FE model updating of concrete building frame from multiple starting points;
* Before running the code, please add "_Shared" folder at the upper level into MATLAB path. 
* The code first loads the data for updating the frame model by calling "LoadStructure.m". Then the code calls the functions in the "_Shared" folder to perform FE model updating.


"LoadStructure.m"
* This is a function called by "ConcreteBuildingFrame.m" to load the influence matrices of updating variables, mass matrix from "ConcreteBuildingFrame.mat".  The function then builds the nominal as well as actual stiffness matrices.


"ConcreteBuildingFrame.mat" 
* This is a file including the influence matrices of updating variables, mass matrix, and measured DOFs.


"\SAP_Model"
* This is a folder including the SAP2000 model files of the concrete building frame.


"\Results" 
* This is a folder including the model updating results of the concrete building frame and the code to plot the figures.

