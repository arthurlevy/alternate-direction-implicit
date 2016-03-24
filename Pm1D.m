function TmFinal = Pm1D(Tinit,dt, L)

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