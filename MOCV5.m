%% Author 
%Name: Mayur Nunkoo
%Role: Founder/President/Propulsion Exec. @Carleton Propulsion Club
%Date: 2023-05-19
%Revision: 
%1-Adapted for multi-wallpoint
%2-Code Revision
%3-Whole Script re-written into smaller for loops
%% Clean
clear 
close all
clc
%% Global Parameters 
k = 1.4; %Specific Heat Ratio, Ideal Gas Assumption --> Needs Correction for mixture
Ma_e = 2.4%Exit Mach number---> needs correction
n = 7; % Number of characteristic lines
gamma = 1.4; %Specific Heat ratio ---> needs correction
nsave = n;
%% Symmetry Point
%labelling
k = 0;
symm(1) = 1;
zeta = n+1;
stopsymm = n-1;
for i = 1:stopsymm
    symm(i+1) = symm(i) +(zeta-k);
    k = k+1;
end
symm = symm';
%Find the corresponding angles for each LRC on the symmetry axis
[mach,nu,mu] = flowprandtlmeyer(gamma,Ma_e,'mach')
nuMe = nu;
thetamax = nuMe/2;
thetaA1 = 0.375; %Probably should be corrected
del = (thetamax-thetaA1)/(stopsymm);
thetasymm(1,1) = thetaA1;
for i = 1:stopsymm
    thetasymm(i+1,1) = thetasymm(i,1) + del;
    kminus(i,1) = 2*thetasymm(i,1);
    kplus(i,1) = 0;
end
kminus(i+1,1) = 2*thetasymm(i+1,1);
kplus(i+1,1) = 0;
%% Points on the first characteristic line
char1theta(1) = thetaA1;
char1nu(1) = thetaA1;
[mach,nu,mu] = flowprandtlmeyer(gamma,char1nu(1),'nu')
char1mu(1) = mu;
char1theta = thetasymm;
char1nu = char1theta;
%point a definition
mu_a = char1mu(1);
theta_a = thetaA1;
nsave = n;
for i = 1:nsave
  [mach,nu,mu] = flowprandtlmeyer(gamma,char1nu(i,1),'nu')  
   char1mu(i,1) = mu;
end
nsave = nsave-1
%% Points on Remaining Characteristic Lines
for j = 2:n
    for i = 1:nsave
       adder(1) = 0;
       s = 2*thetasymm(j:end,1)
       kminus(i,j) = s(i)
       s2 = -2*thetasymm(j:end,1)
       kplus(i,j) = -2*thetasymm(j,1)
       char1nu(i,j) = 0.5*(kminus(i,j)-kplus(i,j))
       [mach,nu,mu] = flowprandtlmeyer(gamma,char1nu(i,j),'nu')  
       char1mu(i,j) = mu;
    end
    nsave = nsave-1
end
%% Fixing Angles for A-Internal pt
for i =1:n
    char1theta(i,1) = thetasymm(i)
end
%For the remaining points:
char1theta(1,2:n) = 0
nsave = n-2;
for j = 2:n
    for i =1:nsave
        char1theta(i+1,j) = char1theta(i,j) + del;
    end
    nsave = nsave-1
end
%% Finding slope of line connecting a and points on first characteristic line
nsave=n;
for j = 1:n
    for i = 1:nsave
      slopeapt(i,j) = 0.5*(char1theta(i,j) + char1theta(i,j))-0.5*(char1mu(i,j)+char1mu(i,j));
    end
    nsave = nsave-1
end
%% Inter slope between points i and i+1 on first char line
nsave = n-1;
for j = 1:n
    for i = 1:nsave
      slopeinter1(i,j) = 0.5*(char1theta(i+1,j) + char1theta(i,j))+0.5*(char1mu(i+1,j)+char1mu(i,j));
    end
    nsave = nsave-1;
end
%% Coordinates of all the symmetry points
yt = 1;%%%%%%%%MUST CHANGE FOR APPROPRIATE THROAT RADIUS
x0 = 0;
y0 = yt;
for i = 1:n
    xsymm(i) = -y0/tand(slopeapt(i,1))
    ysymm(i) = 0
end
xsymm = xsymm'
ysymm = ysymm'
%% Interior coordinates of the first line
xinter(1,1) = -y0/tand(slopeapt(1,1));
yinter(1,1) = 0;
for i = 2:n
   xinter(i,1) = (yinter(i-1,1)-y0-tand(slopeinter1(i-1,1))*xinter(i-1,1))/(tand(slopeapt(i,1))-tand(slopeinter1(i-1,1)))
   yinter(i,1) = (tand(slopeapt(i,1))*xinter(i,1) + y0);
end
%% Interior points of remaining lines
ns = n-1
for j=2:n
    sav = 0.5*(char1theta(2,j-1))-0.5*(char1mu(2,j-1) + char1mu(1,j));
    slopsymm = tand(sav);
    xinter(1,j) = (xinter(2,j-1)*slopsymm - yinter(2,j-1))/slopsymm;
    yinter(1,j) = 0;
    for i = 2:ns
        a = 0.5*(char1theta(i+1,j-1)+ char1theta(i,j)) - 0.5*(char1mu(i+1,j-1)+char1mu(i,j));
        slope1 = tand(a);
        b = 0.5*(char1theta(i-1,j)+ char1theta(i,j)) + 0.5*(char1mu(i-1,j)+char1mu(i,j));
        slope2 = tand(b);
        x3 = xinter(i+1,j-1);
        y3 = yinter(i+1,j-1);
        x9 = xinter(i-1,j);
        y9 = yinter(i-1,j);
        xinter(i,j) = (slope1*x3 - slope2*x9 + y9 - y3)/(slope1-slope2);
        yinter(i,j) = slope1*(xinter(i,j) - x3)+y3;
        
    end
    ns = ns-1
    
end
%% Wall point- first
m78 = 0.5*(char1theta(n-1,1)+ char1theta(n,1)) + 0.5*(char1mu(n-1,1)+char1mu(n,1));
slopem78 = tand(m78);
ma8 = thetamax;
slopema8 = tand(ma8)
xwall(1) = ((-slopem78*xinter(n,1))+(yinter(n,1))-y0)/(slopema8 - slopem78)
ywall(1) = xwall(1)*slopema8 +y0

%% Remaining Wall points
nx = n-1;
for j = 2:n
   if nx>1
        if j == 2
            x14 = xinter(nx,j);
            y14 = yinter(nx,j);
            x8 = xwall(j-1);
            y8 = ywall(j-1);
            m1415 = 0.5*(char1theta(nx,j) + char1theta(nx-1,j)) + 0.5*(char1mu(nx,j) + char1mu(nx-1,j));
            slopem1415 = tand(m1415);
            m815 = 0.5*(char1theta(nx,j-1) + thetamax);
            slopem815 = tand(m815);
            %nx = nx-1;
            xwall(j) = (-slopem1415*x14 + y14 - y8 +x8*slopem815)/(slopem815-slopem1415);
            ywall(j) = slopem815*xwall(j) - x8*slopem815 + y8;
        else
            x14 = xinter(nx,j);
            y14 = yinter(nx,j);
            x8 = xwall(j-1);
            y8 = ywall(j-1);
            m1415 = 0.5*(char1theta(nx,j) + char1theta(nx-1,j)) + 0.5*(char1mu(nx,j) + char1mu(nx-1,j));
            slopem1415 = tand(m1415);
            m815 = 0.5*(char1theta(nx,j-1) + char1theta(nx,j-2))
            slopem815 = tand(m815);
            %nx = nx-1; 
            xwall(j) = (-slopem1415*x14 + y14 - y8 +x8*slopem815)/(slopem815-slopem1415);
            ywall(j) = slopem815*xwall(j) - x8*slopem815 + y8;
        end
        nx = nx-1
   else
       break
   end
end

%% Last Wall Point
thetaexit = 0;
[mach,nu,mu] = flowprandtlmeyer(gamma,Ma_e,'mach')
mache = mu;
m3435 = thetaexit+mache;
slopem3435 = tand(m3435);
m3335 = 0.5*(char1theta(2,n-1));
slopem3335 = tand(m3335);
x34 = xinter(1,n-1);
y34 = 0;
x33 = xwall(n-1);
y33 = ywall(n-1);
xwall(n) = (x34*slopem3435 - y34 - x33*slopem3335 + y33)/(slopem3435-slopem3335)
ywall(n) = xwall(n)*slopem3435 -x34*slopem3435 +y34;

%% Plot
figure(1)
plot(xwall,ywall,'-*')
xlabel('Distance from throat (m)')
ylabel('Distance from longitudinal axis (m)')
title('Wall points generated through MOC')
axis([0 max(xwall) 0 max(ywall)])


%% Save to file
xwall = xwall'
ywall = ywall'
arr = [xwall,ywall]
%Write all the symmetry point labels on the excel file
columnames = {'X wall Point', 'Y Wall Point'}
filename = 'WALLPOINTS.csv'
fileID = fopen(filename, 'w')
fprintf(fileID, '%s,', columnames{1:end-1}); % Write all column names except the last one
fprintf(fileID, '%s\n', columnames{end}); % Write the last column name and a new line character
dlmwrite(filename,arr,'-append')


