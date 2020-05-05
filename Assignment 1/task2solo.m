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
    
    
    
disp('task 2')
disp(E_aero)
disp(E_roll)
disp(E_mass)
disp(E_cycle)