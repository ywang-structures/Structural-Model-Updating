
load ConcreteBuildingFrame

%% Nominal structure
N = size(E_1,1);
K_j{1} = sparse(E_1);
K_j{2} = sparse(E_2);
K_j{3} = sparse(E_3);
K_j{4} = sparse(E_4);
K_j{5} = sparse(E_5);
K_j{6} = sparse(E_6);
K_j{7} = sparse(E_7);
K_j{8} = sparse(E_8);
K_j{9} = sparse(E_9);
K_j{10} = sparse(E_10);
K_j{11} = sparse(E_11);
K_j{12} = sparse(E_12);

clear E_*;

n_alpha = length(alpha_act);

% Rebar stiffness is not updated, but contributes to K0;
%
% In the SAP model, rigid link elements are usually used to connect nodes
% of a concrete element and nodes of a rebar element at the same segment.
% When a rigid link element is used to connect two nodes, SAP provides some
% DOFs of these two nodes as negative numbers for identifying the fixated
% relationship with the link. This prevents a measured DOF to explicitly
% exist in the DOF list.  To bypass the difficulty, when a measured DOF is
% associated with a node connected with a link, a non-rigid link element
% with very high stiffness is allocated, so that the measured DOF explicitly
% exists as a positive DOF in the list (instead of as a negative number).
% Four such links are allocated in total, at nodes C2, C3, G2, and G3, in
% order to have their longitudinal (x-direction) DOF properly represented.
K0 = Influence_Link + Influence_Rebar;
for i = 1:n_alpha
    K0 = K0 + K_j{i};
end
% Build the actual stiffness matrix starting from K0; then add the
% influence of each parameter change
K_act = K0;
for i = 1:n_alpha
    K_act = K_act + K_j{i} * alpha_act(i);
end


