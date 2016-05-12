
function Tfinal =  OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, oneStepProperties)
%solves the 1D heat transfer problem on one single step dt

%% Discretisation
Np = length(Tinit);

Nel=Np-1; Le = L/Nel;

%% solver options
toleranceT = 1e-9; %degC characteristic tolerance on temperatures
toleranceResidual = 1e-12;% characetristic tolerance for residual

options = optimset('Display','none',...
    'TolFun', toleranceResidual, 'TolX', toleranceT,...
    'Algorithm', 'Levenberg-marquardt');%,...
    %'Jacobian', 'on');

%% solve
Tfinal = fsolve (...
    @(T) residual(T, dt, Le, oneStepProperties, Tinit, hinf,hsup,Tinf,Tsup),...
    Tinit, options);
end


%residual function to be zeroed
function [residual, jacobian] = residual(T, dt, Le, oneStepProperties, Tinit,...
    hinf,hsup,Tinf,Tsup)

M_dt = capacity_matrix(T,dt,Le,oneStepProperties);
K = conductivity_matrix(T,Le,oneStepProperties);
[RHS, dRHS_dT] = flux_BC(T,hinf, hsup, Tinf, Tsup);

residual = M_dt * (T-Tinit)  +  K * T  +  RHS';

jacobian = M_dt + K + dRHS_dT;% ...
 %   + ;
end


%mass matrix assembling (M/dt)
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
function [RHS, dRHS_dT] = flux_BC(T,hinf, hsup, Tinf, Tsup)
Np = length(T);

RHS(1)  =  hinf * (T(1) - Tinf);
RHS(Np) =  hsup * (T(end) - Tsup);

dRHS_dT = sparse([1,Np],[1,Np],[hinf, hsup]);
end
