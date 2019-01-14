Contents under \FourStoryStructureExperiment:

"FourStoryStructureExperiment.m" 
* This is the main function to perform FE model updating of four-story structure from 1000 start points;

"FourStoryStructure.mat"
* This is the data file including the experiment data on the four-story structure
	'test_frequency' variable - a 4x1 double array including the identified natural frequencies of the four-story structure
	'test_mode' variable - a 4x4 double array including the identified mode shapes of the four-story structure

"\StructrualDynamics" 
* This is a folder including necessary functions to build the four-story structure model

"\Results" 
* This is a folder including the model updating results of the four-story structure and the code to plot the result figures in the paper
* Please run "PostProcess_Exp.m" to plot the result figures in the paper
* Followings are the descriptions of each result file:
	1. Exp_form1_Jacon_interior-point.mat  - Case A: applying interior-point algorithm on MAC value formulation 
	2. Exp_form2_Jacon_interior-point.mat  - Case B: applying interior-point algorithm on eigenvector difference formulation 
	3. Exp_form3_Jacon_interior-point.mat  - Case C: applying interior-point algorithm on modal dynamic residual formulation
	Each *.mat file contains following variables:
		'alpha' variable - a 4x1000 double array including the updated alpha calculated from 1000 starting points
		'alpha0' variable - a 4x1000 double array including the 1000 starting points for alpha
		'exit_flag' variable - a 1000x1 double array including the exit flags for the 1000 optimization processes
		'fval' variable - a 1000x1 double array including the objective function values at the optimal points of the 1000 optimization processes
		'gradient' variable - a 4x1000 double array including the gradient of the objective function at the optimal points
		't' variable - a 1x1000 double array including the time consumed by the 1000 successful optimization processes
		't_waste' variable - a 1x1 double array including the time consumed by the faild optimization processes
		'expModes' variable - a 1x1 structure including the experimental results. Refer relevant code for details
		'structModel' variable - a 1x1 structure including the structure information. Refer relevant code for details
		'optimzOpts' variable - a 1x1 structure including the optimization options. Refer relevant code for details
		'updatingOpts' variable - a 1x1 structure including the model updating options. Refer relevant code for details