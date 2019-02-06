function [M] = assemblyM(M_i, N)

% function [M] = assemblyM(M_j, N)
%
%   (c) Yang Wang, Xinjun Dong, Dan Li (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% This function assemblies the mass matrix for shear frame structures.
%
% Input:
%   M_i (N x 1)  - a vector with values of floor mass 
%   N - the number of stories in the shear frame structure
%
% Output:
%	M -  the mass matrix of the shear frame structure

[r, c] = size(M_i);
N_ = max(r, c);

if N_ ~= N
    error('\nWrong input for creating the stiffness matrix!');
else
    M = diag(M_i);
end

end


