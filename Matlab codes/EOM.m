%% Solving the Equation of Motion with a significance axisymetry
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%% DE Integration %%%%%%%%%%%%%%%%%%
% Initial Conditions
H0 = 100000;
tang0 = -pi/6;
v0 = 3400;
alpha0 = pi/20;
omegax0 = 0.45;
theta0 = -pi/180;
y0 = [v0, tang0, H0, omegax0, theta0, alpha0];

% Calling the Solver
tspan = [0 93];
sol = ode45(@EOMsys,tspan,y0);

%%%%%%%%%%%%%%%% Chartting the Results %%%%%%%%%%%%%%%%%
% Velocity
figure(1); plot(sol.x,sol.y(1,:)); xlabel('Time [sec]'); ylabel('V [m/s]')
% Tang angle
figure(2); plot(sol.x,sol.y(2,:).*57.3); xlabel('Time [sec]'); ylabel('Tang [deg]')
% Height
figure(3); plot(sol.x,sol.y(3,:)); xlabel('Time [sec]'); ylabel('H [m]')
% Omegax
figure(4); plot(sol.x,sol.y(4,:)); xlabel('Time [sec]'); ylabel('Omegax [rad/sec]')
% Theta
figure(5); plot(sol.x,sol.y(5,:).*57.3); xlabel('Time [sec]'); ylabel('Theta [deg]')
% Alpha
figure(6); plot(sol.x,sol.y(6,:).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')


%%%%%%%%%%%%%%% Differential Equations %%%%%%%%%%%%%%%
function out = EOMsys(t,x)
t
v = x(1);      tang=x(2);   h=x(3);
omegax = x(4); theta=x(5);  alpha=x(6);

%% Constants and Unknowns
g = 3.86;  % Mars Gravity Acceleration
rho = 0.019;% Mars atmoshper density
R = 3390000;  % Mars Radius

% Lander Mass-Geometric Characteristics
s = 1.5174;  % Area of the Lander middle-section
L = 2;  % Lander Length?, 1.06
m = 576;  % Lander Mass
Ix = 270;  Iz = 443;
Ixz= 0.0;  Ixy= 0.0;
I = Iz; % Used For normalization

% aerodynamic
% Cxa = ; % coefficient of Drag aerodynamic force D in Speed coordinates
% Cx = ;  % coefficient of Drag aerodynamic force D in Body coordinates
% Cyn = ; % coefficient of Lift aeridynamic force L in body with axicemetry
[Cx, Cy, Mz] = aero_file(alpha*57.3);
Cxa = Cx ./ cos(alpha);
Cyn = Cy;
Mzn = Mz;
% alpha*57.3

% Aerodynamic moments with assymetry
mzn = 0.5; % Coefficient of recoverd aerodynamic moment 
mzf = 0.0; % Coefficients of small aerodynamic coefficients based on the assymetry of KA
myf = 0.5; % Coefficients of small aerodynamic lift based on the assymetry of KA

dz = 0.4; % Z displacment of CG in OXYZ
dy = 0.0; % Y displacment of CG in OXYZ
epsilon = 1; % Parameter related to the ammount of acymetry of mass,inertia and aerodynamic

mxsf = 0; % mxf: momment of shape assymetry
mxcf = 0; % mxf = mxsf*sin(phi) + mxcf*sin(phi);


%% Equations of Moving Center of Gravity
q = 0.5*rho*v^2;

dVdt = -Cxa*q*s/m - g*sin(tang);
dTangdt = -cos(tang)*(g-v^2/(h+R))/v;
dHdt = v*sin(tang);

% [v tang h omegax theta alpha]
%% Supported Equations
dqdt = rho*v*dVdt;
omega= sqrt(mzn*q*s*L*cot(alpha)/I);
omegaa = sqrt(((Ix/I)*omegax)^2/4 + omega^2);
if omegax > 0
   omega12 = ((Ix/I)*omegax)/2 + omegaa
end 
if omegax <= 0
   omega12 = ((Ix/I)*omegax)/2 - omegaa
end

if omegax > 0
   m1a = -(((1+Ix/I)*omegax-3*omega12)/(2*omegaa))*(omega^2/mzn)*(myf-Cx*dz/L)*tan(alpha)-omega12*omega^2*tan(alpha)^2*(mxcf+Cyn*dz/L)/(2*omegaa*mzn) + (Ixz/(2*I*omegaa))*(omegax*omega12*(omegax+omega12*tan(alpha)^2)-omegax^2*(omegax+2*omegaa)); 
   m2a = -(((1+Ix/I)*omegax-3*omega12)/(2*omegaa))*(omega^2/mzn)*(mzf+Cx*dy/L)*tan(alpha)+omega12*omega^2*tan(alpha)^2*(mxsf+Cyn*dy/L)/(2*omegaa*mzn) + (Ixy/(2*I*omegaa))*(omegax*omega12*(omegax+omega12*tan(alpha)^2)-omegax^2*(omegax+2*omegaa)); 
end
if omegax <= 0
   m1a = -((1+Ix/I)*omegax-3*omega12)/(2*omegaa)*(omega^2/mzn)*(myf-Cx*dz/L)*tan(alpha)-omega12*omega^2*tan(alpha)^2*(mxcf+Cyn*dz/L)/(2*omegaa*mzn) - (Ixz/(2*I*omegaa))*(omegax*omega12*(omegax+omega12*tan(alpha)^2)-omegax^2*(omegax-2*omegaa)); 
   m2a = -((1+Ix/I)*omegax-3*omega12)/(2*omegaa)*(omega^2/mzn)*(mzf+Cx*dy/L)*tan(alpha)+omega12*omega^2*tan(alpha)^2*(mxsf+Cyn*dy/L)/(2*omegaa*mzn) - (Ixy/(2*I*omegaa))*(omegax*omega12*(omegax+omega12*tan(alpha)^2)-omegax^2*(omegax-2*omegaa)); 
end
ma = sqrt(m1a^2 + m2a^2);
theta1 = asin(m1a/ma); % theta1 = acos(-m2a/ma);

mx1a = -(omega^2/mzn)*(mxsf+Cyn*dy/L)*tan(alpha) - (Ixy/I)*omega12^2*tan(alpha);
mx2a = -(omega^2/mzn)*(mxcf+Cyn*dz/L)*tan(alpha) - (Ixz/I)*omega12^2*tan(alpha);
mxa = sqrt(mx1a^2 + mx2a^2);
theta2 = asin(-mx1a/mxa); % theta2 = acos(mx2a/mxa);

Fa = omega12^2/cos(alpha)^2 + (omegax*Ix/I - omega12)*(omegax*Ix/I - 2*omega12);

%% Equations of Moving Around Center of Gravity
dOmegaxdt = -epsilon*mxa*sin(theta+theta2)*(I/Ix);
dThetadt = omegax - omega12;
% dOmegadt = epsilon*omega*dqdt/(2*q);
if omegax > 0
   dAlphadt = epsilon*(tan(alpha)*omega^2*dqdt/(4*omegaa^2*q) - ma*cos(theta+theta1)/(2*omegaa))/(4*omegaa^2/Fa);
end 
if omegax <= 0
   dAlphadt = epsilon*(tan(alpha)*omega^2*dqdt/(4*omegaa^2*q) + (ma/omega^2)*omega^2*cos(theta+theta1)/(2*omegaa))/(4*omegaa^2/Fa);
end

out = [dVdt;dTangdt;dHdt;dOmegaxdt;dThetadt;dAlphadt];
% out = [dVdt;dTangdt;dHdt];
end
