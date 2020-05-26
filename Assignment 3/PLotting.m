% LAbb3 main 

run('init_HEV_ECMS')
tic
sim('HEV_ECMS.mdl')
toc

yyaxis left
plot(SOC.time,SOC.data)
ylabel('SOC [%]')
yyaxis right
plot(V_z.time,V_z.data)
ylabel('Velocity [m/s]')
xlabel('Time [s]')
title('NEDC MAN')