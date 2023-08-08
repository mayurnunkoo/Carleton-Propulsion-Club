%% Author
%Name: Mayur Nunkoo
%Role: Founder/President/Propulsion Lead CPC

%% clean
clear all
close
clc
%% Global Parameters Function call
T0 = 3030 %Stagnation Temperature (K)
k = 1.4  %Specific Heat Ratio
P_exit = 101325 %Exit ambient pressure [equal expansion][Pa]
P0 = 1e6 %Stagnation Pressure [Pa]
Rs = 8.314 %Ideal Gas constant [J/K.mol]
F = 1000 %Thrust force [N]
%% Exit Conditions
%Exit Mach Number Calculation
Ma_exit = ((2/(k-1))*((P0/P_exit)^((k-1)/k) - 1))^0.5
%Exit Temperature
Te = T0/(1+((k-1)/2)*Ma_exit^2)
%Exit Speed of sound
c_exit = sqrt(k*Rs*Te)
%Exit velocity
v_exit = Ma_exit*c_exit
%Area Ratio 
A = 1+((k-1)/k)*Ma_exit^2
B = 1+(k-1)/2
Area_Ratio = sqrt((A/B)^((k+1)/(k-1)))*(1/Ma_exit)
%Stagnation density
rho0= P0/Rs*T0
%Exit Density
rho_e = rho0/((1+((k-1)/2)*Ma_exit^2)^(1/(k-1)))
%mass flow rate of propellant
mdot = F/v_exit
%Exit Area
A_e = mdot/(rho_e*v_exit)
%Throat Area
A_t = A_e/Area_Ratio
%% Creating DataFrame
input_Parameters = {'Stagnation Temperature', 'Specific Heat Ratio', 'Exit Ambient Pressure', 'Stagnation Pressure','Ideal Gas Constant','Thrust'}
exit_parameters = {'Exit Mach Number', 'Exit Temperature', 'Exit Speed of sound','Exit velocity','Area Ratio', 'Stagnation Density','Exit Density','Mass Flow Rate of Propellant','Exit Area','Throat Area'} 
Parameters = [input_Parameters,exit_parameters]
Magnitudes = [T0,k,P_exit,P0,Rs,F,Ma_exit,Te,c_exit,v_exit,Area_Ratio,rho0,rho_e,mdot,A_e,A_t]
units = {'K','','Pa','Pa','J/K.mol','N','','K','m/s','m/s','','kg/m^3','kg/m^3','kg/s','m','m'}
df = table(Parameters',Magnitudes',units','VariableNames',{'Performance Parameters','Magnitude','Units'})
disp(df)
%% Export df to CSV
filename = 'CPC_Performance_Sheet.csv'
writetable(df,filename)


%% Final Remarks
%Please make sure to clone the CPC repo on your local machine, and add any
%scripts you deem useful!