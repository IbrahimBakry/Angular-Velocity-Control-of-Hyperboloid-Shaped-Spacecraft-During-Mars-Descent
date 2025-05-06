%% 
clc
clear
close all

%% Solving the DE
v0 = 3400; 
tang0 = -1/57.3; 
h0 = 120000; 
omegax0 = 22/57.3; 
theta0 = 5/57.3; 
alpha0 = 1/57.3; 
y0 = [v0, tang0, h0, omegax0, theta0, alpha0];
tspan = [0 500];

sol = ode45(@DE,tspan,y0);

% Velocity
figure(1); plot(sol.x,sol.y(1,:)); xlabel('Time [sec]'); ylabel('V [m/s]')
% Tang angle
figure(2); plot(sol.x,sol.y(2,:).*57.3); xlabel('Time [sec]'); ylabel('Tang [deg]')
% Height
figure(3); plot(sol.x,sol.y(3,:)); xlabel('Time [sec]'); ylabel('H [m]')
% Omegax, Theta, Alpha
figure(7)
subplot(311); plot(sol.x,sol.y(4,:)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
subplot(312); plot(sol.x,sol.y(5,:).*57.3); xlabel('Time [sec]'); ylabel('Theta [deg]')
subplot(313); plot(sol.x,sol.y(6,:).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')
% Omega_xr
[Cx, Cy, mzn] = aero_file(sol.y(6,:).*57.3);
omega_xr = sqrt(-0.5.*0.019.*(sol.y(1,:)).^2.*mzn.*1.5174*1.06.*cot((sol.y(6,:)))./443)./sqrt(1-(270/443));
figure(8); plot(sol.x,omega_xr,'--'); hold on
           plot(sol.x,sol.y(4,:)); hold off
xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
legend('OmegaxResonance','Omegax')


%% DE Function
function [dxdt]=DE(t,x)
v = x(1); tang = x(2); h = x(3); omegax = x(4); theta = x(5); alpha = x(6);

%% Constants
g = 3.86;     % Mars Gravity Acceleration
rho = 0.019;  % Mars atmoshper density
R0 = 3390000; % Mars Radius
omega0 = 7.088e-5; % Mars rotation speed

s = 1.5174;   % Area of the Lander middle-section 
L = 1.06;     % Lander Length?, 1.06
m = 576;      % Lander Mass
Ix = 270;  
Iz = 443;
I = (Iz+Iz)/2; % Used For normalization

[Cx, Cy, Mz] = aero_file(alpha*57.3);
Cxa = Cx/alpha; % coefficient of Drag aerodynamic force D in Speed coordinates
Cya = Cy/alpha;
Cyn = Cy; % coefficient of Lift aeridynamic force L in body with axicemetry
mzn = Mz; % Cx: coefficient of Drag aerodynamic force D in Body coordinates

mzf = 0.0; % Coefficients of small aerodynamic coefficients based on the asymmetry of KA
myf = 0.6; % Coefficients of small aerodynamic lift based on the asymmetry of KA
dz = 0.4; % Z displacement of CG in OXYZ, dz/L > 0.3
dy = 0.0; % Y displacement of CG in OXYZ, dy/L > 0.3
epsilon = 0.1; % Parameter related to the amount of asymmetry of mass, inertia and aerodynamic

mxsf = 0; % mxf: moment of shape asymmetry0
mxcf = 0; % mxf = mxsf*sin(phi) + mxcf*sin(phi);

%% Coefficients
Ixd = Ix/I;
dyd = dy/L;
dzd = dz/L;

%% DE of CG
q = 0.5*rho*v^2;
% dvdt = g*(-q*Cxa*s/(m*g) - sin(tang) + omega0^2*(R0+h)*sin(tang)/(m*g));
% dtangdt = (g/v)*(q*Cya*s/(m*g) - cos(tang) + v^2*cos(tang)/(g*(R0+h)) + 2*omega0*v/(m*g) + omega0^2*(R0+h)*cos(tang)/(m*g));
% dhdt = v*sin(tang);
dvdt = g*(-q*Cxa*s/(m*g) - sin(tang));
dtangdt = (g/v)*(- cos(tang)+ v^2*cos(tang)/(g*(R0+h)));
dhdt = v*sin(tang);

%% Sublimentary equations
omega = sqrt(-mzn*q*s*L*cot(alpha)/I);
omegaa = sqrt(Ixd^2*omegax^2/4 + omega^2);
if omegax >= 0, omega12 = Ixd*omegax/2 + omegaa; end
if omegax < 0, omega12 = Ixd*omegax/2 - omegaa; end

% mxa1 = -omega^2*(mxsf+Cyn*dyd)*tan(alpha)/mzn;
% mxa2 = -omega^2*(mxcf+Cyn*dzd)*tan(alpha)/mzn;
% mxa = sqrt(mxa1^2 + mxa2^2);
% theta2 = asin(-mxa1/mxa);
% 
% ma1 = - (((1+Ixd)*omegax - 3*omega12)*omega^2*(myf-Cx*dzd)*tan(alpha) - omega12*omega^2*tan(alpha)^2*(mxcf+Cyn*dzd))/(2*omegaa*mzn);
% ma2 = - (((1+Ixd)*omegax - 3*omega12)*omega^2*(mzf+Cx*dyd)*tan(alpha) + omega12*omega^2*tan(alpha)^2*(mxsf+Cyn*dyd))/(2*omegaa*mzn);
% ma = sqrt(ma1^2 + ma2^2);
% theta1 = asin(ma1/ma);
ma = 0.03;
mxa= 0.03;
theta1 = pi;
theta2 = 0;

mzna = mzn/alpha;
Mzna = mzna*q*s*L;
Fa = - Mzna/I + omega12^2/cos(alpha)^2 + (Ixd*omegax-omega12)*(Ixd*omegax-2*omega12);
% Fa = omega12^2/cos(alpha)^2 + (Ixd*omegax-omega12)*(Ixd*omegax-2*omega12);

disp(['Time = ' num2str(t) ', Alpha = ' num2str(alpha*57.3) ', T1 =' num2str(theta1*57.3) ', T2 =' num2str(theta2*57.3)])

%% Calculating the probability of catching or resonance
% D = omegax - omega12; 
% dbdomega
% dDdomegax = 1 - omega12;
% dDdalpha = 0;
% domegadt = -L*mzn*q*s*dalphadt*(1+cot(alpha)^2)/(2*I*sqrt(-mzn*q*s*L*cot(alpha)/I));
% dDdomega = ;
% 
% if omegax > 0; mf1 =-dDdomegax*(mxa2/Ixd) + dDdalpha*(2*omegaa*ma1/Fa);end
% if omegax < 0; mf1 =-dDdomegax*(mxa2/Ixd) - dDdalpha*(2*omegaa*ma1/Fa);end
% if omegax > 0; mf2 = dDdomegax*(mxa1/Ixd) + dDdalpha*(2*omegaa*ma2/Fa);end
% if omegax < 0; mf2 = dDdomegax*(mxa1/Ixd) - dDdalpha*(2*omegaa*ma2/Fa);end
% 
% b = sqrt(mf1^2 + mf2^2);
% C1 = -dbdalpha*(2*omega*tan(alpha)*domegadt/Fa) + dbdomega*domegadt;
% if omegax > 0; C2 = dbdalpha*(2*ma*omegaa/Fa);end
% if omegax > 0; C2 =-dbdalpha*(2*ma*omegaa/Fa);end
% 
% theta3 = acos(mf1/b);
% theta4 = theta1 - theta3;
% Q = abs(24*C1 + 8*C2*cos(theta4))/(3*sqrt(b));
% a = -dDdalpha*(2*omega*tan(alpha)*domegadt/Fa) + dDdomega*domegadt;
% 
% if Q < 4*pi*abs(a); Pr = Q/(2*pi*abs(a)+Q/2);end
% if Q >=4*pi*abs(a); Pr = 1; end

%% DE about CG
dqdt = rho*v*dvdt;
domegaxdt = -epsilon*mxa*sin(theta+theta2)/Ixd;
dthetadt = omegax - omega12;

if omegax > 0, dalphadt = epsilon*((tan(alpha)/(2*omegaa^2))*(omega/(2*q))*dqdt + ma*cos(theta+theta1)/(2*omegaa))/(Fa/(4*omegaa^2)); end
if omegax < 0, dalphadt = epsilon*((tan(alpha)/(2*omegaa^2))*(omega/(2*q))*dqdt - ma*cos(theta+theta1)/(2*omegaa))/(Fa/(4*omegaa^2)); end

dxdt = [dvdt; dtangdt; dhdt; domegaxdt; dthetadt; dalphadt];
end

