clear data
global data
data.dt = 1;
%Engine
% roll 0.6/(0.3*0.3)/1500 % JW/Rw^2 / MV
data.e      = [0.4];
data.pme0   = [0.1*10^6];
data.Vd     = [1.497*10^-3] ;
data.Je     = [0.2];
%...
%Fuel
data.Hlhv   = [44.6*10^6];
data.rho_f  = [737.2];

%Battery
data.Ri     = [0.65];
data.Uoc    = [300];
data.Q      = [6.5];%*3600 %Ah battery capacity
%....
%Electric motor
data.n_em   = [0.9];
%...
%Limits
data.Temmax = [400]; %T_em
data.Pemmax = [50*10^3];%P_em
data.wdotmax = [300];% acceleration max
data.Ticemax = [115]; % engine torque
data.Imax = [200]; % max discharge
data.w_max = [5000*2*pi/60]; % max speed

% NM*rad/s - U/U*s
