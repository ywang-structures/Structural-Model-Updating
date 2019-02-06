Contents under \SAP_Model

“ConcreteBuildingFrame.sdb”  SAP2000 model file of the concrete building frame.


“DOFInfo.mat”   MATLAB data file with DOF information from the SAP2000 model.  
	'DOFlist' variable – a 67x1 cell array including the node labels of the SAP 2000 FE model;  
	'DOF' variable – a 67x6 double array listing the DOFs for each node in the stiffness and mass matrices. The i-th row in DOF variable corresponds to the i-th cell in DOFlist. For example, the 2nd row in 'DOF' contains the six DOFs for node 'A2' (the 2nd cell in DOFlist).  

“seltMeasDOFs.m”    MATLAB code that selects the measured DOFs from all the DOFs. The code loads DOFInfo.mat, and generates measDOFs.mat.


“measDOFs.mat”    MATLAB data file containing the measured DOFs in the SAP2000 model.

