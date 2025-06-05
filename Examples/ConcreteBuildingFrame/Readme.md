
Contents under \ConcreteBuildingFrame:

The concrete building frame example is a full-scale test structure in the Structural Engineering and Materials Laboratory on the Georgia Tech campus. Model updating is performed entirely through simulation models. Detailed descriptions of the structure and model updating studies can be found in the companion monograph, ["Formulation and application of SMU - an open-source MATLAB package for structural model updating."](https://github.com/ywang-structures/Structural-Model-Updating/blob/master/Formulation%20and%20application%20of%20SMU%20%E2%80%93%20an%20open-source%20MATLAB%20package%20for%20structural%20model%20updating.pdf)


"ConcreteBuildingFrame.m" 
* This is the main function to perform FE model updating of concrete building frame from multiple starting points;
* Before running the code, please add "SMU" folder at the upper level into MATLAB path. 
* The code first loads the data for updating the frame model by calling "LoadStructure.m". Then the code calls the functions in the "SMU" folder to perform FE model updating.


"LoadStructure.m"
* This is a function called by "ConcreteBuildingFrame.m" to load the influence matrices of updating variables, mass matrix from "ConcreteBuildingFrame.mat".  The function then builds the nominal as well as actual stiffness matrices.


"ConcreteBuildingFrame.mat" 
* This is a file including the influence matrices of updating variables, mass matrix, and measured DOFs.


"\SAP_Model"
* This is a folder including the SAP2000 model files of the concrete building frame.


"\Results" 
* This is a folder including the model updating results of the concrete building frame and the code to plot the figures.
