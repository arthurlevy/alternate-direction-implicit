function TzFinal = Pz1D(Tinit,hinf, hsup, Tinf, Tsup, dt)

%% number of dofs through thickness
nbnodes = length(Tinit);

nb_plies = 8;

nbnodes_per_ply = nbnodes / nb_plies;

%array of integer indices givinfg ply number
ply_number = fix((1:nbnodes+1-1e-3)'/nbnodes_per_ply);

%array of boolean giving ply orientation 0°
ply_orientation = mod(ply_number,2);

% through thcknes properties
PzProperties.k_a = 0.5*ply_orientation + 2*(1-ply_orientation);
PzProperties.k_b = 0*ply_orientation + 0*(1-ply_orientation);
PzProperties.rho_a = 1000*ply_orientation + 1000*(1-ply_orientation);
PzProperties.rho_b = 0*ply_orientation + 0*(1-ply_orientation);
PzProperties.cp_a = 1000*ply_orientation + 1000*(1-ply_orientation);
PzProperties.cp_b = 0*ply_orientation + 0*(1-ply_orientation);



L = 2e-3;

TzFinal = ...
    OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, PzProperties);


end