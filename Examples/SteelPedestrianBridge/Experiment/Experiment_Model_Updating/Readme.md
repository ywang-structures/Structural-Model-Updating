Contents under \SteelPedExp:

"SteelPedestrianBridge_Exp.m" 
* This is the main function to perform FE model updating of steel pedestrian bridge from 100 starting points.

"ExpModeInfo_FEM.mat"
* This is the data file including the experiment modal properties of the steel pedestrian bridge.
	'lambdaExp' variable - a 5x1 double array including the identified eigenvalues of the steel pedestrian bridge
	'numTest' variable - a 5x1 double array including the number of modal analysis results for each identified eigenvalue
	'psiExp_m' variable - a 58x5 double array including the identified mode shapes of the steel pedestrain bridge
	'numVldCH' variable - a 58 x 5 double array including the number of modal analysis results for each entry of the identified mode shapes
	'std_lambdaExp' variable - a 5x1 double array including the standard deviation of the identified eigenvalues
	'std_psi' variable - a 58 x 5 double array including the standard deviation for each entry of the identified mode shapes

"SteelPedestrianBridgeExp.mat" 
* This is a file including the influence matrices of updating variables, mass matrix, and measured DOFs.

"LoadStructure.m"
* This is a function called by "SteelPedestrianBridgeExp.m" to load the influence matrices of updating variables and mass matrix from "SteelPedestrianBridgeExp.mat". The function then builds the nominal stiffness matrix.


"\Results" 
* This is a folder including the model updating results of the steel pedestrian bridge and the code to plot the result figures in the paper
* Please run "PostProcess.m" to plot the best solution result figures in the paper
* Please run "loc_Min.m" to plot the local minimum solution result figure in the paper
* The three model updating result files are:
  1. SteelPedBridg_form1_Jacon_interior-point.mat  - Case A: applying interior-point algorithm on MAC value formulation 
  2. SteelPedBridg_form2_Jacon_interior-point.mat  - Case B: applying interior-point algorithm on eigenvector difference formulation 
  3. SteelPedBridg_form3_Jacon_interior-point.mat  - Case C: applying interior-point algorithm on modal dynamic residual formulation
  Each *.mat file contains following variables:
  	'alpha' variable - a 17x100 double array including the updated alpha calculated from 100 starting points
	'alpha0' variable - a 17x100 double array including the 100 starting points for alpha
	'exit_flag' variable - a 100x1 double array including the exit flags for the 100 optimization processes
	'fval' variable - a 100x1 double array including the objective function values at the optimal points of the 100 optimization processes
	'gradient' variable - a 17x100 double array including the gradient of the objective function at the optimal points
	't' variable - a 1x100 double array including the time consumed by each of the 100 successful optimization processes
	't_waste' variable - a 1x1 double array including the time consumed by the failed optimization processes
	'expModes' variable - a 1x1 structure including the experimental results. Refer relevant code for details
	'structModel' variable - a 1x1 structure including the structure information. Refer relevant code for details
	'optimzOpts' variable - a 1x1 structure including the optimization options. Refer relevant code for details
	'updatingOpts' variable - a 1x1 structure including the model updating options. Refer relevant code for details