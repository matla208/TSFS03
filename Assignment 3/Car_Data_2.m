% Car_data_2

%Constants

g = 9.81; % m/s2
Pa = 1.18; % kg/m3
Hl = 44.6*10^6 ;%J/kg
pl = 737.2; % kg/m3

%Car Data 
Je = 0.2 ;% kgm^2
Temax = 115;%Nm
Vd = 1.497*10^-3 ;%m^3
%Vd = 2.2579*10^-3;% more power 
nr = 2 ;% rev per stroke
me = 1.2*10^-3 ;% Kg/W

% Willians approx
e = 0.4; 
pme0 = 0.1*10^6; % Mpa

%
Q0 = 6.5 *3600; %Ah battery capacity
Uoc = 300;
Imax = 200; % +- <--
Ri = 0.65;
mbatt = 45;
%
n = 0.9;
Temmax = 400;
Pemmax = 50*10^3;%kom ihåg +-
mem = 1.5*10^-3; % kg / W

Pptmax = 90.8*10^3;

cd = 0.32;
cr = 0.015;
Af = 2.31;
mv = 1500;%-120 (90.8*1.2-72.2) ;
rw = 0.3;
Jw = 0.6;
% Tabell 1.3
egb = 0.98;
CL = 0.3*10^3; %kW
Gear_ratio = [13.0529 8.1595 5.6651 4.2555 3.2623];

% Tabell 2.2

Ne = [0 800:1:5000];    % Allowed engine speed range [rpm]
wdot = 300;          % Maximum acceleration in the engine [rad/s2]

