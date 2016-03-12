function TmFinal = Pm1D(Tinit,dt, L)

Np = length(Tinit);

PmProperties.k_a = 0.5;
PmProperties.k_b = 0;
PmProperties.rho_a = 1000;
PmProperties.rho_b = 0;
PmProperties.cp_a = 1000;
PmProperties.cp_b = 0;

hinf = 0;
Tinf = 0;
hsup = 0;
Tsup = 0;

TmFinal = ...
    OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, PmProperties);

end