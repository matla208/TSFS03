% Drive Cycles 
%Task 1
clear all
close all
clc
%run('Car_data')
%load('NEDC_MAN');
load('FTP_75')

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
