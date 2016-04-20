function TzFinal = Pz1D(Tinit,hinf, hsup, Tinf, Tsup, dt)
%solves one problem Pz over one single time step dt 

%properties are supposed to depend linearly on temperature such that:
% k = k_a + T*k_b
%rho = rho_a + T*rho_b ...

PzProperties.k_a = 0.42;
PzProperties.k_b = 0*1e-3;
PzProperties.rho_a = 1600;
PzProperties.rho_b = -0*0.2;
PzProperties.cp_a = 800;
PzProperties.cp_b = 0*2.25;

L = 2e-3;

TzFinal = ...
    OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, PzProperties);


end