% This function reads the DOFInfo output from SAP2000 and select the
% measured DOFS
clc;clear;

% “DOFInfo.mat”   MATLAB data file with DOF information from the SAP2000
% model. 
%   'DOFlist' variable – a 67x1 cell array including the node labels of
%      the SAP 2000 FE model;  
%   'DOF' variable – a 67x6 double array listing the DOFs for each node in
%      the stiffness and mass matrices. The i-th row in DOF variable
%      corresponds to the i-th cell in DOFlist. For example, the 2nd row in
%      'DOF' contains the six DOFs for node 'A2' (the 2nd cell in
%      DOFlist).  
load DOFInfo

%% Pick the node in FE model where sensor is instrumented
% Nodes whose longitudinal vibration is instrumented
logtdNode = [1:9 15:18 10 19 25:27 30:32 21 33 23 34 66 13 29 14 35];  
% Nodes whose vertical vibration is instrumented
vertNode = [5 11 24 10 28 29 36 33 9 12 14 37 27 20 35 66];         

for i = 1 : length(logtdNode)
    measDOFs(i) = DOF(logtdNode(i),1);
end
for i = 1 : length(vertNode)
    measDOFs(i + length(logtdNode)) = DOF(vertNode(i),3);
end

% Remove DOFs associated with A1, E1, I1, the three fixed nodes at the column bases.
measDOFs = measDOFs( measDOFs ~= 0 );

save measDOFs.mat measDOFs 
