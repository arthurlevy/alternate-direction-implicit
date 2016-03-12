%solves the heat transfer problem using the additive decomposition method
%[Hoang et al 2016] which is better called locally one dimensional. 
function T_history = LocallyOneDHeatTransfer

%% time discretization
t_final = 1;
dt = 0.1;
tspan = 0:dt:t_final;
nb_time_step = length(tspan);

%% space discretization
nbnodePm = 25;
nbnodePz = 15;

%% initial temperature
T = 250 * ones(nbnodePm, nbnodePz);

%% preallocation
T_half = zeros(nbnodePm, nbnodePz);
%temperature history
T_history = zeros(nbnodePm, nbnodePz, nb_time_step);

%% geometry
Lf = 5e-3; %flange length
R = 2*2e-3 ; %midplane radius
L = 2*Lf + pi/2 * R ; %L shape midplane length.

xspan = linspace(0,L,nbnodePm);

%% set Pm boundary values
RTC = 1e-4; % thermal contact resistance
hmould = 1/RTC;% exchange coeff with mould
hair = 20; %exchange coeff with air
Tmould = 200; %mould temperature
Tair = 20; %air temeprature

hinf = hmould * ones(nbnodePm, 1);
hsup(xspan>Lf) = hair ;
hsup(xspan<=Lf) = hmould;

Tinf = Tmould * ones(nbnodePm, 1);
Tsup(xspan>Lf) = Tair ;
Tsup(xspan<=Lf) = Tmould;

for i_time = 1:nb_time_step
    display ( [num2str(i_time), ' : time = ', num2str(tspan(i_time))])
    
    %% for each in plane position
    for iter_position_m = 1:nbnodePm  
        %solve Pz : T_half is T_{n+1/2}
        T_half(iter_position_m,:) =...
            Pz1D(T(iter_position_m,:)',...
            hinf(iter_position_m), hsup(iter_position_m),...
            Tinf(iter_position_m), Tsup(iter_position_m),...
            dt);
    end%solved all Pz problems
     
    %% averaging : Tm_half is <T_{n+1/2}>
    Tm_half = mean(T_half,2);
    
    %% in plane Pm problem Tm_one is <T_{n+1}>
    Tm_one = Pm1D(Tm_half,dt,L);
    
    %% recomposition
    T = bsxfun(@plus, Tm_one - Tm_half, T_half);
     
    %% store solution
    T_history(:,:,i_time) = T;
    
end% time loop

surf(T_history(:,:,end));

end