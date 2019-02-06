function [K] = assemblyK(K_i, N)

% function [K] = assemblyK(K_j, N)
%
%   (c) Yang Wang, Xinjun Dong, Dan Li (all rights reserved)
%       School of Civil and Environmental Engineering
%       Georgia Institute of Technology
%       2018
%
% Revision: 1.0
%
% This function assemblies the stiffness matrix for shear frame structures.
%
% Input:
%   K_i (N x 1)  - a vector with values of inter-story stiffness 
%   N - the number of stories in the shear frame structure
%
% Output:
%	K -  the stiffness matrix of the shear frame structure

[r, c] = size(K_i);
N_ = max(r, c);

if N_ ~= N
    error('\nWrong input for creating the stiffness matrix!');
else
    K = zeros(N, N);
    K(1,1) = K_i(1);
    for i = 2:N
        K(i-1, i-1) = K(i-1, i-1) + K_i(i);
        K(i-1, i) = K(i-1, i) - K_i(i);
        K(i, i-1) = K(i, i-1) - K_i(i);
        K(i, i) = K(i, i) + K_i(i);
    end
end

end


