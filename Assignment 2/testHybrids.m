%% Test case
clear all
close all
clc
format shortg
load('EUDC_MAN_DDP')
%load('City_MAN_DDP')
run('Car_Data_2')

global V_z
global Gear
Gear = zeros(length(G_z),1);
Gear(G_z==0) = 1;
Gear(G_z==1) = Gear_ratio(1);
Gear(G_z==2) = Gear_ratio(2);
Gear(G_z==3) = Gear_ratio(3);
Gear(G_z==4) = Gear_ratio(4);
Gear(G_z==5) = Gear_ratio(5);
xtot = trapz(V_z);
%%
% Uncomment the lines as applicable
% CaseToRun='parallel';
CaseToRun='series';
%CaseToRun='testParallel';
%CaseToRun='testPseries';

switch CaseToRun,
 case 'parallel',
  % Set up the problem:
  % Select discretization for the SOC.
  % Assign cost for the final state.
  T_grid    = T_z;
  SOC_grid  = 0.49:0.0001:0.51;
  %SOC_grid = 0.5; single value simulates conventional as we cant use
  %battery.
  finalCost = zeros(size(SOC_grid));
  finalCost(SOC_grid<0.5) = inf;% makes it charge sustaining 
  % Solve the problem
  tic
  [value,SOC_path]=dynProg1D(@parallelHybrid,T_grid,SOC_grid,finalCost);
  toc
  % bra att kunna clearvars variabelnamn
  % Postprocess and analyze the result
  
  %fuel=value(1,(SOC_grid==0.5))/rho f/x tot*1e8;
  SOC_start = 0.50
  pos = find(SOC_grid==SOC_start)
  %SOC_value(1) = SOC_start
  
 for i = 1:length(T_grid)-1
     SOC_value(i) = SOC_grid(SOC_path(i,pos));
     mdot_f(i) = value(i,pos);
     pos = find(SOC_grid==SOC_value(i));
     
 end
 figure(1)

 plot(SOC_value)
 ylabel('State of Charge [%]')
 xlabel('Time [s]')
 figure(2)
 
 plot(mdot_f)
 ylabel('Cost in fuel mass [Kg]')
 xlabel('Time [s]')
%   mn = length(value)-1;
%   Soc_value = [Soc_start SOC_grid(pos)] ;
 case 'series',
  % Set up the problem:
  % Select discretization for the SOC and engine speed.
  % Assign cost for the final states.
  T_grid    =T_z ;
  Ne_vector = [0 800:200:5000];
 SOC_vector = 0.4:0.0001:0.6;
%   SOC_vector = repmat(0.5,length(Ne_vector),1)'; % ingen batteri användning
  
  SOC_grid  = repmat(SOC_vector,length(Ne_vector),1)';% make them the same to get all possible final points
  
  Ne_grid   = repmat(Ne_vector,length(SOC_vector),1);
  finalCost = zeros(length(SOC_vector),length(Ne_vector));
  finalCost(SOC_grid<0.5) = inf;
  finalCost(Ne_grid ~= 0) = inf;
  % Solve the problem§
  tic
  [value,SOC_path,Ne_path]=dynProg2D(@seriesHybrid,T_grid,SOC_vector,Ne_vector,finalCost);
  toc
  % Postprocess and analyze the result
 
  SOC_start = 0.50;
  Ne_start = 0;
  SOC_pos = find(SOC_vector==SOC_start )%& Ne_vector == Ne_start);
  %SOC_value(1) = SOC_start
  Ne_pos = find(Ne_vector== Ne_start )%& SOC_vector == SOC_start);
  F100 = value(1,find(SOC_vector==SOC_start ),find(Ne_vector== Ne_start ))*100/(pl*10^-3*xtot*10^-3);
 for i = 1:length(T_grid)-1
     SOC_value(i) = SOC_grid(SOC_pos,Ne_pos);
     Ne_value(i) = Ne_grid(SOC_pos,Ne_pos);
     
     mdot_f(i) = value(i,SOC_pos,Ne_pos);
     
     SOC_pos = SOC_path(i,SOC_pos,Ne_pos);
     Ne_pos = Ne_path(i,SOC_pos,Ne_pos);
 end
 
 figure(1)

 plot(SOC_value)
 ylabel('State of Charge [%]')
 xlabel('Time [s]')
 
 figure(3)
 plot(Ne_value)
 ylabel('engine_speed [%]')
 xlabel('Time [s]')
 figure(2)
 
 plot(mdot_f)
 ylabel('Cost in fuel mass [Kg]')
 xlabel('Time [s]')
  
  
 case 'testParallel',
   % Test the arc costs of the parallel hybrid to see that they are reasonable
   % Add the code for the tests that you want to make
  
   
 case 'testSeries',
   % Test the arc costs of the series hybrid and
   % Add the code for the tests that you want to make
   
 otherwise,
  error('Unknown case')
end