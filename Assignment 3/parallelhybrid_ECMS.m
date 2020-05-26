function [U]=parallelhybrid_ECMS(x)
global data
%Given inputs from the driving cycle and controller this script should return
%the torques on the electric motor and combustion engine.


w_ice=x(1);
dw_ice=x(2);
T_req=x(3);
lambda=x(4);
%% Your code here..........
if w_ice <= 0
    T_ice = 0;
    T_em = 0;
else
    %PF = Treq*w_w = Tm+ Tice = Pc*w_e + Pb/w_e
    
    Ib = -data.Imax:0.1:data.Imax; % Imax = 200
    Pb = (data.Uoc^2-(data.Uoc-2*Ib*data.Ri).^2)/(4*data.Ri);
    Pem = Pb.*data.n_em.^(sign(Ib));
    
    T_em = Pem/w_ice ; % w_ice = w_em since parallell
    
    T_ice = T_req - T_em;
    %T_ice(T_ice < 0) = 0;
    mdot_f = w_ice/(data.e*data.Hlhv)*(T_ice + data.pme0*data.Vd/(pi*4)+data.Je*dw_ice);
    
    Pf = mdot_f*data.Hlhv;
    Pf(Pf<0) = 0;pme0
%     P_ice = T_ice*w_ice;
    P_ech = data.Uoc*Ib;
    
    
    
    % Controller using hamiltonian
    H = Pf + lambda*P_ech; % Pb rätt?
    
    H((T_ice+data.Je*dw_ice)>data.Ticemax) = inf;
%     H(abs(Ib)>data.Imax) = inf;
    H(abs(Pem)>data.Pemmax) = inf;
    H(abs(T_em)>data.Temmax) = inf;
    H(dw_ice > data.wdotmax) = inf;
    H(w_ice > data.w_max) = inf; 
    H(T_ice < 0)= inf;
    
    [Y Idx] = min(H);
    
    T_ice = T_ice(Idx);
    T_em  = T_em(Idx);

end

U=[T_ice;T_em];
end  


