% Task 2'
clear all
close all
clc
format shortg
%load('EUDC_MAN_DDP')
load('City_MAN_DDP')
run('Car_Data_2')

global V_z
global Gear
Gear = zeros(length(G_z),1);
Gear(G_z==0) = 1;
Gear(G_z==1) = Gear_ratio(1);
Gear(G_z==2) = Gear_ratio(2);
Gear(G_z==3) = Gear_ratio(3);
Gear(G_z==4) = Gear_ratio(4);
Gear(G_z==5) = Gear_ratio(5);

%% Vehicle calc
Xtot = trapz(V_z);
V_avg =[(V_z(1:end - 1) + V_z(2:end))/2;0];
A_z = [(V_z(2:end)-V_z(1:end - 1));0];
Fa = 1/2*cd*Pa*Af*V_avg.^2;
Fr = cr*mv.*g;
Ft = mv.*A_z + Fa + Fr;

FG = Ft + Jw/(rw^2).*A_z;


Tw = FG * rw; 

%% Gear_box

Gear = zeros(length(G_z),1);
Gear(G_z==0) = inf;
Gear(G_z==1) = Gear_ratio(1);
Gear(G_z==2) = Gear_ratio(2);
Gear(G_z==3) = Gear_ratio(3);
Gear(G_z==4) = Gear_ratio(4);
Gear(G_z==5) = Gear_ratio(5);

w_g = V_z/(2*rw*pi)*2*pi; % avg v ?
w_e = Gear.*w_g;
w_e(isnan(w_e))=0;

Tee = Tw./Gear.*(egb.^(-sign(Tw))) ;
Te = Tee.*sign(Tw);
% no clutch consideration
%% Engine
Te(Te>Temax) = inf;

% räkna ut PC behöver ne = pme/pmf ne = e enligt mahdi
ne = e;
Pc = w_e.*Te/ne;
mdot_f = Pc/Hl;

%% Hybrid motor

P2 = Tee.*w_e;

P1 = P2.*n.^(-sign(P2));

%% Power converter

% for i 1:length(p2)
% % 
% %Uoc(t)?Ri(t)·I2(t) = U2(t)
% % Ibatery fpr U2 and U2 for the power feedback
% U2 = P2/I2;
% I2= Uoc - U2 /Ri
% I2 = 

Soc = 0.5;
t = 1;
I21 = Soc*Q0/t
U2 = Uoc - Ri*I21;
I2 = P1/U2
h =  1
Q = I2*h
Soc = Q/q0
delta_Soc = I2*Q0/h
% SoC(k+1) = SoC(k) + delta(SoC) <-- aka soc from the current I2 value


%% test
% cost1 =parallelHybrid([15 16], 0.5, [0.49 0.498 0.50 0.501 0.51])*10^4
% cost2 =parallelHybrid([4 5], 0.5, [0.49 0.498 0.50 0.501 0.51])*10^4
% cost3 =parallelHybrid([59 60], 0.5, [0.49 0.498 0.50 0.501 0.51])*10^4

%% Hybrid
t = [15 16; 4 5; 59 60;158 159];

for i= 1:4
    costp1(i,:) = parallelHybrid([t(i,:)], 0.5, [0.49 0.498 0.50 0.501 0.51]);% uppgift 1 parallel hybrid
end

costp12 = parallelHybrid([59 60], 0.5, [0.49 0.498 0.50 0.501 0.51])
%% series
seris =seriesHybrid([4 5], 0.5, [49 49.8 50 50.2 51]*1e-2, 3e3, 3e3)
seris2 =seriesHybrid([4 5], 0.5, [0.499 0.5], 0, [0 8 20]*1e2)
seris2 =seriesHybridtemp([4 5], 0.5, [0.499 0.5], 0, [0 8 20]*1e2)
%% costs series
t = [15 16; 4 5; 59 60;158 159]

    costs1 = [seriesHybrid([t(1,:)], 0.5, [49 49.8 50 50.2 51]*1e-2, 3e3, 3e3);
                seriesHybrid([t(2,:)], 0.5, [49 49.8 50 50.2 51]*1e-2, 3e3, 3e3);
                seriesHybrid([t(3,:)], 0.5, [49 49.8 50 50.2 51]*1e-2, 3e3, 3e3);
                seriesHybrid([t(4,:)], 0.5, [49 49.8 50 50.2 51]*1e-2, 3e3, 3e3)]
    costs2 = [seriesHybrid([t(1,:)], 0.5, [0.5], 3e3, [0 2 3 5]*1e3);
                seriesHybrid([t(2,:)], 0.5, [0.5], 3e3, [0 2 3 5]*1e3);
                seriesHybrid([t(3,:)], 0.5, [0.5], 3e3, [0 2 3 5]*1e3);
                seriesHybrid([t(4,:)], 0.5, [0.5], 3e3, [0 2 3 5]*1e3)];

    costs3 = [seriesHybrid([t(1,:)], 0.5, [0.499 0.5], 0, [0 8 20]*1e2);
                seriesHybrid([t(2,:)], 0.5, [0.499 0.5], 0, [0 8 20]*1e2);
                seriesHybrid([t(3,:)], 0.5, [0.499 0.5], 0, [0 8 20]*1e2);
                seriesHybrid([t(4,:)], 0.5, [0.499 0.5], 0, [0 8 20]*1e2)];


%% 1.1 e
% battery + electrical motor
% maximum contribution = 
electrical_weight = mem * Pptmax + 45