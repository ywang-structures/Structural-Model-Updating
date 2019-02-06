
load SteelPedestrianBridge
   
%% Nominal structure
N = size(E_1,1);
K_j{1} = sparse(E_1);
K_j{2} = sparse(E_2);
K_j{3} = sparse(E_3);
K_j{4} = sparse(E_4);
K_j{5} = sparse(E_5);
K_j{6} = sparse(E_6);
K_j{7} = sparse(E_t2);
K_j{8} = sparse(E_t3);
K_j{9} = sparse(E_t4);
K_j{10} = sparse(E_t5);
K_j{11} = sparse(E_t6);
K_j{12} = sparse(k_y1);
K_j{13} = sparse(k_z1);
K_j{14} = sparse(k_y2);
K_j{15} = sparse(k_z2);

clear E_* k_*;

n_alpha = length(alpha_act);
K0 = sparse( zeros( N, N));
for i = 1:n_alpha
    K0 = K0 + K_j{i};
end

% Build the actual stiffness matrix starting from K0; then add the
% influence of each parameter change
K_act = K0;
for i = 1:n_alpha
    K_act = K_act + K_j{i} * alpha_act(i);
end


