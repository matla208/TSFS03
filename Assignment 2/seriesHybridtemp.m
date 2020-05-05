function costVector=seriesHybridtemp(t_vec,SOC_start,SOC_final,Ne_start,Ne_final);

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
  
run('Car_Data_2')
%format shortG

global drivecycle;
drivecycle = {'EUDC_MAN_DDP.mat','City_MAN_DDP.mat'};

load(drivecycle{2})


We_lim = 800*2*pi/60;


Ne1 = Ne_start*2*pi/60; 
Ne2 = Ne_final*2*pi/60;

t1 = t_vec(1);
t2 = t_vec(2);
h = t2 - t1;
v = (V_z(t1) + V_z(t2))/2;
a = (V_z(t2) - V_z(t1))/h;

Fa = 0.5*Pa*Af*0.32*v^2;
Fg = mv*g*sin(0);
Fr = cr*mv*g;
Fi = mv*a;
Fw = Jw/(rw^2)*a;

Ft = Fa + Fg + Fr + Fi + Fw;


Pw = Ft*v;
P_em = Pw*n^(-sign(Pw));

delta_soc = SOC_final-SOC_start;
I2 = -delta_soc*Q0/h;
Pm = (Uoc^2-(Uoc-2*I2*Ri).^2)/(4*Ri);

P_gen = P_em-Pm;


P_ice = P_gen/n;
We_dot = Ne2 - Ne1;

Pe = P_ice/h;
We = (Ne1+Ne2)/2;
Te = Pe'./(We+(We==0));
Te(Te<0) = 0;

mf = We/e/Hl.*(Te+pme0*Vd/(4*pi)+Je*We_dot);

mf(abs(P_gen)>Pemmax) = inf;
mf(abs(We_dot)>wdot) = inf;
mf(abs(I2)>Imax) = inf;
mf(Te>Temax) = inf;


We_lim = repmat(We_lim,1,size(Te,1))';
mf(Te>0 & We<We_lim) = inf;


mf(mf<0) = 0;


mf(Ne1==0 & Ne2>We_lim) = inf;

costVector = mf;