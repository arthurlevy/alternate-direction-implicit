function TzFinal = Pz1D(Tinit,hinf, hsup, Tinf, Tsup, dt)

%% number of dofs through thickness
% nbnodes = length(Tinit);
% 
% nb_plies = 8;
% 
% nbnodes_per_ply = nbnodes / nb_plies;
% 
% %array of integer indices givinfg ply number
% ply_number = fix((0:nbnodes-1)'/nbnodes_per_ply);
% 
% %array of boolean giving ply orientation 0°
% ply_orientation = mod(ply_number,2);
% 
% % through thcknes properties
% PzProperties.k_a = 4*ply_orientation + 0.42*(1-ply_orientation);
% PzProperties.k_b = 7.5e-3*ply_orientation + 1e-3*(1-ply_orientation);
% PzProperties.rho_a = 1600*ones(nbnodes,1);
% PzProperties.rho_b = -0.2*ones(nbnodes,1);
% PzProperties.cp_a = 800*ones(nbnodes,1);
% PzProperties.cp_b = 2.25*ones(nbnodes,1);


PzProperties.k_a = (4+0.42) / 2;
PzProperties.k_b = (7.5e-3+1e-3) / 2;
PzProperties.rho_a = 1600;
PzProperties.rho_b = -0.2;
PzProperties.cp_a = 800;
PzProperties.cp_b = 2.25;

L = 2e-3;

TzFinal = ...
    OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, PzProperties);


end