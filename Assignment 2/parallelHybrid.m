function costVector=parallelHybrid(t_vec,SOC_start,SOC_final);
global V_z
global Gear
run('Car_Data_2')
% costVector=parallelHybrid(t_vec,SOC_start,SOC_final);
%
% Function for calculating the all arc-costs from one node to all other possible nodes.
% 
% This function is called from dynProg1D.m
%  
% Inputs:
%  t_vec     - 1x2 matrix, with the start and stop time for the interval.
%  SOC_start - A single start value for SOC during the interval.
%  SOC_final - Vector (from the dicretization) with all possible final values for the
%             interval.
%
% Output:
%  costVector - Vector with the costs for all arcs from SOC_start. 
  
% Version 1.0,  2008-06-30 Lars Eriksson


% Implement your parallel hybrid model and calculate the arc costs below.
% temp test
% t_vec = [60 61];
% SOC_start = 0.5;
% SOC_end = [0.49 0.498 0.50 0.501 0.51];

h = t_vec(2)-t_vec(1); %<--- should be 1 everywhere
V_avg = (V_z(t_vec(1)) +V_z(t_vec(2)))/2; 
A_z = (V_z(t_vec(2))-V_z(t_vec(1)))/h;% for each time step h = 1
if V_avg == 0
    mdot_f(SOC_final==SOC_start)= 0;
    mdot_f(SOC_final~=SOC_start) = inf;
    costVector = mdot_f;
    return
end
%Fa = 1/2*cd*Pa*Af*V_avg^2;
Fa = (1/2)*0.32*Af*Pa*V_avg^2;
%display(Fa)
Fr = cr*mv.*g;
Ft = mv.*A_z + Fa + Fr;
FG = Ft + Jw/(rw^2).*A_z;

gear = Gear(t_vec(1));
w_g = V_avg/(2*rw*pi)*2*pi; % avg v ?
w_e = gear.*w_g;
w_gdot = A_z/rw;
w_edot =gear.* w_gdot;
Tw = FG * rw; 
Tg = Tw./gear.*(egb.^(-sign(Tw)))  ;

%Te = Tg.*sign(Tw);
PG = Tg.*w_e ;

%Pc = w_e.*Te/e; % <-- fel jao Pe beror på battery
%mdot_f = Pc/Hl;

% Battery stuff

DSoC = SOC_final-SOC_start;
I2 = -DSoC*Q0/h;
Pm = (Uoc^2-(Uoc-2*I2*Ri).^2)/(4*Ri); %eq 4.60 book 
Pb = Pm.*n.^(sign(Pm));
Tb = Pb./w_e; % neccesary?
%power split for cost function

Pc = PG-Pb; % whats left to reach required effect at the gear
Te = Pc/w_e;
mdot_f = w_e/(e*Hl)*(Te+ pme0*Vd/(pi*4)+Je*w_edot);
% check against limits
mdot_f(Te>Temax) = inf;
mdot_f(abs(I2)>Imax) = inf; 
mdot_f(abs(Pb)>Pemmax) = inf;
mdot_f(Tb>Temmax) = inf;
mdot_f(mdot_f<0) = 0;
costVector = mdot_f;
