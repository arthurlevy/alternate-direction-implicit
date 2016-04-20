
function Tfinal =  OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, oneStepProperties)
%solves the 1D heat transfer problem on one single step dt

%% Discretisation
Np = length(Tinit);

Nel=Np-1; Le = L/Nel;

%% define residual function
residual = @(currT) ...
    capacity_matrix(currT,dt,Le,oneStepProperties) * (currT - Tinit) ...
    +  conductivity_matrix(currT,Le,oneStepProperties) * currT ...
    + flux_BC(currT,hinf, hsup, Tinf, Tsup)' ;
    
%% solve
toleranceT = 1e-3; %degC characteristic tolerance on temperatures
toleranceResidual = 1e-3;% characetristic tolerance for residual

options = optimset('Display','iter',...
    'TolFun', toleranceResidual, 'TolX', toleranceT,...
    'Algorithm', 'Levenberg-marquardt');

Tfinal = fsolve (residual, Tinit, options);
end

%mass matrix assembling
function M = capacity_matrix(T,dt,Le,properties)

Np = length(T);
rho = properties.rho_a + properties.rho_b .* T;
cp = properties.cp_a + properties.cp_b .* T;
rhocp_dt = rho .* cp / dt;

% set values
idxcolumn =    [1:Np,    1:Np-1,    2:Np];

idxline   =    [1:Np,   2:Np,    1:Np-1];

values    = Le/6 *...
    [2*rhocp_dt(1) ; 4*rhocp_dt(2:end-1); 2*rhocp_dt(end) ;...
    rhocp_dt(1:end-1);...
    rhocp_dt(2:end)];

%assemble
M = sparse ( idxcolumn  , idxline  , values  ,     Np,Np);
end

%stifness matrix assembling
function K = conductivity_matrix(T,Le,properties)

Np = length(T);

k_ = properties.k_a + properties.k_b .* T;

% set values
idxcolumn = [1:Np, 1:Np-1, 2:Np];

idxline   = [1:Np, 2:Np, 1:Np-1];

values    = 1/Le *...
    [k_(1); 2*k_(2:end-1); k_(end) ;...
    -k_(1:end-1);...
    -k_(2:end)];

%assemble
K = sparse ( idxcolumn  , idxline  , values  ,     Np,Np);
end

%right hand side assembling boundary conditions
function RHS = flux_BC(T,hinf, hsup, Tinf, Tsup)
Np = length(T);

RHS(1) = hinf * (T(1) - Tinf);
RHS(Np) = hsup * (T(end) - Tsup);
end
