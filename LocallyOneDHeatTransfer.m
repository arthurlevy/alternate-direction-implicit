%solves the heat transfer problem using the additive decomposition method
%[Hoang et al 2016] which is better called locally one dimensional or
%alternate direction implicit method. 
function T_history = LocallyOneDHeatTransfer

%% time discretization
t_final = 20;
nb_time_step = 70;
%non constant time mesh
tspan = linspace(0,t_final,nb_time_step).^2 / t_final;

%% space discretization
nbnodePm = 50;
nbnodePz = 20;

%% initial temperature
T = 400 * ones(nbnodePm, nbnodePz);

%% preallocation
T_half = zeros(nbnodePm, nbnodePz);
%temperature history
T_history = zeros(nbnodePm, nbnodePz, nb_time_step);
T_history (:,:,1) = T;

%% geometry
L_f = 20e-2; %flange length
L_poincon = L_f/2; %dye width
shear_edge = 0.004;

xspan = linspace(0,L_f/2,nbnodePm);

%% set Pm boundary values
hmould = 5000;% exchange coeff with mould
hair = 10; %exchange coeff with air
Tmould = 200; %mould temperature
Tair = 70; %air temeprature

hinf(xspan > L_poincon/2+shear_edge) = hmould;
hinf(xspan <= L_poincon/2+shear_edge) = hair;

hsup(xspan > L_poincon/2) = hair ;
hsup(xspan <= L_poincon/2) = hmould;

Tinf(xspan > L_poincon/2+shear_edge) = Tmould;
Tinf(xspan <= L_poincon/2+shear_edge) = Tair;

Tsup(xspan >  L_poincon/2) = Tair ;
Tsup(xspan <=  L_poincon/2) = Tmould;


for i_time = 2:nb_time_step
    display ( [num2str(i_time), ' : time = ', num2str(tspan(i_time))])
    %update dt
    dt = tspan(i_time)-tspan(i_time-1);
    
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
    Tm_one = Pm1D(Tm_half,dt,L_f/2);
    
    %% recomposition
    T = bsxfun(@plus, Tm_one - Tm_half, T_half);
     
    %% store solution
    T_history(:,:,i_time) = T;
    
end% time loop

surf(T_history(:,:,end));

end