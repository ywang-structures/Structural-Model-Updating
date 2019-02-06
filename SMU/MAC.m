function macs = MAC(vec, vecs)

% function mc = MAC(vec, vects) 
%
%   Yang Wang, Xinjun Dong, Dan Li
%   School of Civil and Environmental Engineering
%   Georgia Institute of Technology
%   2018
%
% Revision: 1.0
%
% This function calculates the modal assurance criterion between vect1
% and each column in vects (Reference: R. J. Allemang and D. L. Brown, "A
% correlation coefficient for modal vector analysis," in Proceedings of the
% 1st international modal analysis conference, pp. 110-116, 1982.)
%

if size(vec, 1) ~= size(vecs, 1)
    error('\nMAC: the input vectors must have the same length.');
end

numVecs = size(vecs, 2);
macs = zeros(numVecs, 1);
for i = 1 : numVecs
    vec_i = vecs(:, i);
    macs(i) = (vec' * vec_i) ^ 2 / ( norm(vec)^2 * norm(vec_i)^2 );
end
