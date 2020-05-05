function costVector=seriesHybrid(t_vec,SOC_start,SOC_final,Ne_start,Ne_final);
global V_z
run('Car_Data_2')
% costVector=seriesHybrid(t_vec,SOC_start,SOC_final,Ne_start,Ne_final);
% 
% Function for calculating all the arc-costs during one time interval, from one node to
% all other possible nodes.
% 
% This function is called from dynProg2D.m with the following argument list
%    funName([tVec(t) tVec(t+1)],discX(ii),discX,discY(jj),discY)
%  
% Inputs:
%  t_vec      - 1x2 matrix, with the start and stop time for the interval.
%  SOC_start  - A single start value for SOC during the interval.
%  SOC_final  - Vector (from the dicretization) with all possible final values for the
%               interval.
%  Ne_start   - A single start value for engine speed during the interval.
%  Ne_final   - Vector (from the dicretization) with all possible final values for the
%               interval.
%
% Output:
%  costMatrix - Matrix with the costs for all arcs from SOC_start. 
%               The size of the returned cost-matrix should have the following
%               size:  length(SOC_final) x length(Ne_final)
  
% Version 1.0,  2008-06-30 Lars Eriksson
h = t_vec(2)-t_vec(1); %<--- should be 1 everywhere
V_avg = (V_z(t_vec(1)) +V_z(t_vec(2)))/2; 
A_z = (V_z(t_vec(2))-V_z(t_vec(1)))/h;% for each time step h = 1

Fa = (1/2)*0.32*Af*Pa*V_avg^2;
Fr = cr*mv.*g;
Ft = mv.*A_z + Fa + Fr;
FG = Ft + Jw/(rw^2).*A_z;

%gear = Gear(t_vec(1));
w_e = (Ne_start+Ne_final)/2 *2*pi/60;
%w_e(w_e==0)= 1
w_edot = (Ne_final-Ne_start/h)*2*pi/60;

Tw = FG * rw; 
Tg = Tw.*(n.^(-sign(Tw)))  ; % no gears ty electric motor only 

%Te = Tg.*sign(Tw);
PG = Tg.*V_avg/rw ; % w_e not correlated to wheel

% Battery stuff

DSoC = SOC_final-SOC_start;
I2 = -DSoC*Q0/h;
Pm = (Uoc^2-(Uoc-2*I2*Ri).^2)/(4*Ri); %eq 4.60 book 
% ic motor generates electricity to battery
Pgen = PG-Pm ;% we now only need to know how much ic need to contribute to
% the motor and therefore dont need the battery contribution.
Pc =Pgen/n;% engine effect * n = generated effect in motor
Te = Pc'./(w_e+(w_e==0));
Te(Te<0) = 0;
mdot_f = w_e./(e*Hl).*(Te+ pme0.*Vd/(pi*4)+Je.*w_edot);

%check limits
We_min = 800*2*pi/60 ;% lowest allwoed enigine speed
mdot_f(abs(Pgen)>Pemmax)=inf; 
mdot_f(abs(w_edot)>wdot)=inf;

mdot_f(abs(I2)>Imax) = inf; 
mdot_f(Te>Temax) = inf;
We_min= repmat(We_min,1,size(Te,1))';

mdot_f(Ne_start==0 & Ne_final*2*pi/60>We_min) = inf;

mdot_f(Te>0 & w_e<We_min) = inf; % check if w_e is too low at positve torqu

mdot_f(mdot_f<0) = 0 ;% cant have negative flow
costVector = mdot_f;
%