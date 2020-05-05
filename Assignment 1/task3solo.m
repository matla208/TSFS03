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
% for i = 2:H
%         V_av = (Vi(i)+Vi(i-1))/2;
%         Fa = 1/2*Pa*Af*cd*V_av^2;
%         Fr = cr*mv*g;
%         Fg = mv*g * sin(0);
%         Ft(i) = mv*Ai(i) + (Fa + Fr + Fg);
%         
%         if Ft(i) > 0 && Vi(i)~= 0%~=
%             V3 = V3 + V_av^3;
%             V = V + V_av;
%             AV = AV + Ai(i)*V_av;%((Ai(i)-Ai(i-1))/h)*((Vi(i)+Vi(i-1))/2);% A*V mean
%             
%             T_trac = T_trac + 1*h;
%             %Xtrac = Xtrac + Vi(i)*(Ti(i)-Ti(i-1));
%             
%         end
%        if Ft(i) >0 && Ft(i-1) <= 0
%            
%            T_prop = T_prop + 1;
%        end
% end
%task3
    T_frac = T_trac/(length(T_z)*h);
    F_traca = Alpha2*1/2*Pa*Af*cd;
    F_tracr = Beta2*mv*g*cr;
    F_tracm = Gamma2*mv;
    F_trac = F_traca + F_tracr + F_tracm;
    
  
 
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
