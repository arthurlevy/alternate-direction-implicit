
function Tfinal =  OneDFEMSolverOneStep(Tinit,hinf, hsup, Tinf, Tsup, dt,L, oneStepProperties)
%solves the 1D heat transfer problem on one single step dt

% Discretisation
Np = length(Tinit);

Nel=Np-1; Le = L/Nel;

%% define residual function
residual = @(currT) ...
    capacity_matrix(currT,dt,Le,oneStepProperties) * (currT - Tinit) ...
    +  conductivity_matrix(currT,Le,oneStepProperties) * currT ...
    - flux_BC(currT,hinf, hsup, Tinf, Tsup)' ; 
    

%% solve
%Jacob PAttern
%idxcolumn = [1:Np, 1:Np-1, 2:Np];
%idxline   = [1:Np, 2:Np, 1:Np-1];
%JacobPattern = sparse(idxcolumn, idxline, ones(3*Np-2,1), Np,Np);
options = optimoptions('fsolve','Display','none');%,...
    %'TolFun', tolfun, 'TolX', 1e-2);%,...
%    'JacobPattern', JacobPattern);
Tfinal = fsolve (residual, Tinit, options);

%Linear case
% M = capacity_matrix(Tinit,dt,Le,oneStepProperties);
% K =    conductivity_matrix(Tinit,Le,oneStepProperties);
% RHS = flux_BC(Tinit,hinf, hsup, Tinf, Tsup)';
% Tfinal = (M+K) \ (RHS + M*Tinit); 
end


function M = capacity_matrix(T,dt,Le,properties)

Np = length(T);

rho = properties.rho_a + properties.rho_b * T;
cp = properties.cp_a + properties.cp_b * T;
rhocp_dt = rho .* cp / dt;

%% Mass matrix
% set values
idxcolumn = [1:Np, 1:Np-1, 2:Np];
idxline   = [1:Np, 2:Np, 1:Np-1];
values    = [(2*Le/3)*rhocp_dt ; Le/6*rhocp_dt(1:end-1); Le/6*rhocp_dt(2:end)];

%assemble
M = sparse ( idxcolumn  , idxline  , values  ,     Np,Np);
M(1,1)=rhocp_dt(1)*Le/3;
M(Np,Np)=rhocp_dt(1)*Le/3;
end

function K = conductivity_matrix(T,Le,properties)

Np = length(T);

k_ = properties.k_a + properties.k_b .* T;

%% stiffness
% set values
idxcolumn = [1:Np, 1:Np-1, 2:Np];
idxline   = [1:Np, 2:Np, 1:Np-1];
values    = [(2/Le)*k_ ; -k_(1:end-1)/Le; -k_(2:end)/Le];

%assemble
K = sparse ( idxcolumn  , idxline  , values  ,     Np,Np);
K(1,1)=k_(1)/Le;       K(Np,Np)=k_(end)/Le;
end

function RHS = flux_BC(T,hinf, hsup, Tinf, Tsup)
Np = length(T);

RHS(1) = -hinf * (T(1) - Tinf);
RHS(Np) = hsup * (T(end) - Tsup);
end
