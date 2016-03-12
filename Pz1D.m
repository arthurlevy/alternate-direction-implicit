function TzFinal = Pz1D(Tinit,hinf, hsup, Tinf, Tsup, dt)


PzProperties.k_a = 0.5;
PzProperties.k_b = 0;
PzProperties.rho_a = 1000;
PzProperties.rho_b = 0;
PzProperties.cp_a = 1000;
PzProperties.cp_b = 0;

L = 2e-3;

TzFinal = ...
    OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, PzProperties);


end