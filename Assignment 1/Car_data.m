%Constants

g = 9.81; % m/s2
Pa = 1.18; % kg/m3
qLHV = 44.6*10^6 ;%J/kg
pf = 737.2*10^-3; % kg/m3

%Car Data 
%Body

Mt = 1400;%kg %Total mass
mv = Mt;
Rm = 8; %Fraction of total mass rotating
Af = 1.9 ;%m2 Frontal area
cd = 0.3 ;%-Air drag coeff.
cr = 0.01; % -Roll. res. coeff.
r_w = 0.3 ;%m Wheel radius

%Engine geometry

Cyl = 6 ;%-Cylinders
Stroke = 79.5*10^-3; %mm Stroke
Bore = 82.4*10^-3; %mm Bore
J_e = 0.2; %kg m2 inertia engine

%Performance during traction

Ieff_trac = 0.35 ;%-Indicated engine efficiency
LME_trac = 1.5 ;%bar Loss mean effective pressure

%Performance during idling
Ieff_idle= 0.3 ;%-Indicated engine efficiency
LME_idle= 1.8 ;%bar Loss mean effective pressure
w_e = 750; %rpm Engine speed

% Transmission Data
egb = 0.98;% -Mechanical efficiency
P0gb = 0.3 *10^3;%kW Constant term in losses
G1 = 13.0529 ;%-Ratio gear 1
G2 = 8.1595 ;%-2
G3 = 5.6651 ;%-3
G4 = 4.2555 ;%-4
G5 = 3.2623 ;%-5

% %Parameter intervals
% mv = [1000,2100];% kg Mass mv
% Af =[1,2];% m2 Frontal area
% cd = [0.1,0.4];% -Air drag coeff.
% cr = [0.005,0.02];% -Roll. res. coeff.

% Transmission Data Task 5
%egb = 0.98;% -Mechanical efficiency
%P0gb = 0.3 *10^3;%kW Constant term in losses
G1_a = 2.89 ;%-Ratio gear 1
G2_a = 1.57 ;%-2
G3_a = 1.00 ;%-3
G4_a = 0.69 ;%-4
G5_a = G4_a ;%-5