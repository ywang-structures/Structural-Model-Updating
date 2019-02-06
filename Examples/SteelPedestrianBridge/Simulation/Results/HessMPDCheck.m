% Check positive/negative definiteness of the Hessian matrix at local
% minima
clc;clear;close all

fileID = 1;
updResults = dir('SteelPedBrdg*.mat');
filename = updResults(fileID).name
load(filename)

n_run = find(alpha(1,:) ~= 0);
for i = n_run
    [hess] = HessMPD(alpha(:,i), structModel, expModes, updatingOpts);
    eigMax(i,1) = max(eig(hess));
    eigMin(i,1) = min(eig(hess));
end