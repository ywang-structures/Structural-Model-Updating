
load SteelPedestrianBridgeExp

% Number of DOFs
N = size(K_Frame,1);

K_j{1} = sparse(K_Frame);
K_j{2} = sparse(K_Truss);
K_j{3} = sparse(K_Conc1);
K_j{4} = sparse(K_Conc2);
K_j{5} = sparse(K_Conc3);
K_j{6} = sparse(K_Conc4);
K_j{7} = sparse(K_Conc5);
K_j{8} = sparse(K_Conc6);
K_j{9} = sparse(K_Conc7);
K_j{10} = sparse(K_Conc8);
K_j{11} = sparse(K_Conc9);
K_j{12} = sparse(K_Conc10);
K_j{13} = sparse(K_Conc11);
K_j{14} = sparse(K_y1);
K_j{15} = sparse(K_z1);
K_j{16} = sparse(K_y2);
K_j{17} = sparse(K_z2);

n_alpha = length(K_j);

M0 = sparse(M0);
K0 = sparse( zeros( N, N));
for i = 1:n_alpha
    K0 = K0 + K_j{i};
end





