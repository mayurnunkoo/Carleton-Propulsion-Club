%% Author
%Name: Mayur Nunkoo
%Role: Founder/President/Propulsion Lead CPC

%% clean
clear all
close
clc
%% Global Parameters
T0 = 3030 %Stagnation Temperature (K)
k = 1.4  %Specific Heat Ratio
P_exit = 101325 %Exit ambient pressure [equal expansion][Pa]
P0 = 1e6 %Stagnation Pressure [Pa]
Rs = 8.314 %Ideal Gas constant [J/K.mol]
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

%% Creating DataFrame
Parameters = {'Stagnation Temperature', 'Specific Heat Ratio', 'Exit Ambient Pressure', 'Stagnation Pressure','Ideal Gas Constant'}
Magnitudes = [T0,k,P_exit,P0,Rs]
units = {'K','','Pa','Pa','J/K.mol'}
df = table(Parameters',Magnitudes',units','VariableNames',{'Performance Parameters','Magnitude','Units'})
disp(df)
%% Export df to CSV
filename = 'CPC_Performance_Sheet.csv'
writetable(df,filename)


%% Final Remarks
%Please make sure to clone the CPC repo on your local machine, and add any
%scripts you deem useful!