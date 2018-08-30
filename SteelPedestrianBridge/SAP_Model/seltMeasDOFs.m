% This function reads the DOFInfo output from SAP2000 and selects the
% measured DOFS
clc;clear;

% “DOF.txt”   text file with DOF information from the SAP2000 model.  
%      Column 1 lists the node labels of the SAP 2000 FE model; 
%      Column 2~7 lists the DOFs for each node in the stiffness and mass
%      matrices, i.e. the 2nd row in 'DOF.txt' contains the six DOFs for
%      node '2'.
dof = load('DOF.txt');

%% Pick the node in FE model where sensor is instrumented
% Nodes where biaxial accelerometers are instrumented
biaxlNode = [22 8 26 12 30 20 33] ;
% Nodes where uniaxial accelerometers are instrumented
unaxlNode = [36 7 40 11 44 19 35] ;

numdir = 2;
for i = 1 : length(biaxlNode)
    for j = 1 : numdir
        measDOFs_biaxle( (i - 1) * numdir + j, 1) = dof(biaxlNode(i), 2 + j);
    end
end

for i = 1 : length(unaxlNode)
    measDOFsOrig_unaxle(i,1) = dof(unaxlNode(i),4);
end
measDOFs = [measDOFs_biaxle; measDOFsOrig_unaxle];

measDOFs = measDOFs(measDOFs ~= 0) ;

save measDOFs measDOFs
