% Drive Cycles 
%Task 1
clear all
close all
clc
%run('Car_data')
load('NEDC_MAN');
%load('FTP_75')

% constants
g = 9.81; % m/s2
Pa = 1.18; % kg/m3
qLHV = 44.6 ;%MJ/kg
pf = 737.2;%*1000 % kg/m3


Ti = T_z;
Vi = V_z;
Xtot = Vi(1)*Ti(1);


for i = 2:length(T_z);
    
  Xtot = Xtot + Vi(i)*(Ti(i)-Ti(i-1)); 
    
end

Ai = D_z;
G = G_z;
H = length(T_z);
V3 = 0;
V = 0;
AV = 0;
%T_trac = 0;
%T_idle = 0;
h = 1; %step size
% configs 
mv = [1000 2100];
Af = [2 1 ];
cd = [0.4 0.1];
cr = [0.02 0.005];
for j = 1:2;
    V3 = 0;
    V = 0;
    AV = 0;
    T_trac = 0;
    T_idle = 0;
    Xtrac = 0;
    
    mvj = mv(j);
    Afj = Af(j);
    cdj = cd(j);
    crj = cr(j);
    
    for i = 2:H
        V_av = (Vi(i)+Vi(i-1))/2;
        A_z = (Vi(i)-Vi(i-1))/h; % changed acceleration h = 1
        Fa = 1/2*Pa*Afj*cdj*V_av^2;
        Fr = crj*mvj*g;
        Fg = mvj*g * sin(0);
        Ft = mvj*A_z + (Fa + Fr + Fg);
        
        if Ft > 0 && V_av~= 0%~=
            V3 = V3 + V_av^3;
            V = V + V_av;
            AV = AV + A_z*V_av;%((Ai(i)-Ai(i-1))/h)*((Vi(i)+Vi(i-1))/2);% A*V mean
            
            T_trac = T_trac + 1*h;
            Xtrac = Xtrac + Vi(i)*(Ti(i)-Ti(i-1));
        end
        if V_av == 0
%             disp(Ai(i))
%             disp(Vi(i))
            T_idle = T_idle + 1*h;
        
        end
        
    end
    % Task 1
    F_traca(j) = 1/Xtot*1/2*Pa*Afj*cdj*V3;
    alpha(j) = V3/Xtot;

    F_tracr(j) = 1/Xtot*mvj*g*crj*V;
    Beta(j) = V/Xtot;

    F_tracm(j) = 1/Xtot*mvj*AV;
    Gamma(j) = AV/Xtot;
    
    T_frac(j) = T_trac/(length(T_z)*h);
    T_idles(j) = T_idle%H-T_trac;%T_idle;
end
disp('task1.1  abg')
disp(alpha)
disp(Beta)
disp(Gamma)

disp(T_frac)
disp(T_idles)

%% task 2
%ändra t_frac för config 3
mv = 1400 ;
Af = 1.9;
cd = 0.3;
cr = 0.01;
    Alpha2 = mean(alpha)
    Beta2 = mean(Beta)
    Gamma2 = mean(Gamma)
 
    T_frac = T_trac/(length(T_z)*h);
    F_traca = Alpha2*1/2*Pa*Af*cd;
    F_tracr = Beta2*mv*g*cr;
    F_tracm = Gamma2*mv;
    F_trac = F_traca + F_tracr + F_tracm;
    
    E_aero = F_traca*100/(qLHV*pf);
    E_roll = F_tracr*100/(qLHV*pf);
    E_mass = F_tracm*100/(qLHV*pf);
    E_cycle = F_trac*100/(qLHV*pf);
    
    
%     T_frac(j) = T_trac/(length(T_z)*h);
%     T_idles(j) = H-T_trac;%T_idle;
%     %Task 2 liter per 100
%     F_trac(j) = F_traca(j) + F_tracr(j) + F_tracm(j);
%     E_aero(j) = F_traca(j)*100/(qLHV*pf);
%     E_roll(j) = F_tracr(j)*100/(qLHV*pf);
%     E_mass(j) = F_tracm(j)*100/(qLHV*pf);
%     E_cycle(j) = F_trac(j)*100/(qLHV*pf);
    
disp('task 2')
disp(E_aero)
disp(E_roll)
disp(E_mass)
disp(E_cycle)

disp('task 3')
%% task 3
% Task 3 fuel consumption estimate
run('Car_data')
load('NEDC_MAN')
% Constants
Vd_c = pi*Bore^2*Stroke/4;
VD_e = Cyl*Vd_c;
cm_idle = (w_e*2*pi/60)*Stroke/pi;
cm_mean = 5.9;
Ratio_aux = 0.08;%power consumed by auxiliary units is 8 %
Ti = T_z;
Vi = V_z;
Ai = D_z;
G = G_z;
H = length(T_z);
V3 = 0;
V = 0;
AV = 0;
h = 1;
T_trac = 0;
T_idle = 0;
Xtrac = 0;
T_prop = 0;
Xtot = trapz(Vi);
for i = 2:H
        V_av = (Vi(i)+Vi(i-1))/2;
        A_z = (Vi(i)-Vi(i-1))/h;
        Fa = 1/2*Pa*Af*cd*V_av^2;
        Fr = cr*mv*g;
        Fg = mv*g * sin(0);
        Ft(i) = mv*A_z + (Fa + Fr + Fg);
        
        if Ft(i) > 0 && V_av~= 0%~=
            V3 = V3 + V_av^3;
            V = V + V_av;
            AV = AV + A_z*V_av;%((Ai(i)-Ai(i-1))/h)*((Vi(i)+Vi(i-1))/2);% A*V mean
            
            T_trac = T_trac + 1*h;
            %Xtrac = Xtrac + Vi(i)*(Ti(i)-Ti(i-1));
            
        end
       if Ft(i) >0 && Ft(i-1) <= 0
           
           T_prop = T_prop + 1;
       end
end
    T_frac = T_trac/(length(T_z)*h);
    F_traca = Alpha2*1/2*Pa*Af*cd;
    F_tracr = Beta2*mv*g*cr;
    F_tracm = Gamma2*mv;
    F_trac = F_traca + F_tracr + F_tracm;
    
  
 
%task 3
P_trac = F_trac*mean(Vi)/T_frac;
P1 = 1/egb*(P_trac+P0gb);
%Paux = P1*Ratio_aux;
Tprop = length(T_z)/T_prop;

%Ec = 1/2*mv*cm_mean^2;
%Ecycle = Ec/Tprop;

Pe = P1 /(1-0.08); %Ecycle = 0
Pme = 16 * Pe/(pi*Cyl*Bore^2*cm_mean);

ne = Pme * Ieff_trac/(Pme + LME_trac*10^5);
Pf = T_frac *Pe/ne;
E_vehicle = Pf/(qLHV*pf);
V_v100 = E_vehicle/mean(Vi)*100*10^3%100/100
% idling stuff
Pmf0 = LME_idle*10^5/Ieff_idle;
Vf_idle = Pmf0*(VD_e/(qLHV*pf))*(cm_idle/(4*Stroke));
% idling time
T_idle = sum(T_idles)/2;%length(T_z)-T_trac;%-T_prop;
Vf_idle100 = T_idle*Vf_idle*100/(Xtot/(10^3))

E_prop = Vf_idle100 + V_v100

%task 4
P_idle = Pmf0*VD_e*(cm_idle/(4*Stroke))
