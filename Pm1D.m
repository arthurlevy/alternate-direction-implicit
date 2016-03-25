function TmFinal = Pm1D(Tinit,dt, L)
%solves the problem Pm over one single time step dt 

%properties are supposed to depend linearly on temperature such that:
% k = k_a + T*k_b
%rho = rho_a + T*rho_b ...

%the in plane conductivity is the average over all plies.

PmProperties.k_a = (4+0.42) / 2;
PmProperties.k_b = (7.5e-3+1e-3) / 2;
PmProperties.rho_a = 1600;
PmProperties.rho_b = -0.2;
PmProperties.cp_a = 800;
PmProperties.cp_b = 2.25;

hinf = 0;
Tinf = 0;
hsup = 0;
Tsup = 0;

TmFinal = ...
    OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, PmProperties);

end