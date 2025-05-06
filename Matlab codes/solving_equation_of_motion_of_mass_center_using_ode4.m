% Initial Conditions
H0 = 100000;
tang0 = -pi/6;
v0 = 3400;
y0 = [v0, tang0, H0];

% Calling the Solver
tspan = [0 300];
sol = ode113(@EOMsys,tspan,y0);

%%%%%%%%%%%%%%%% Chartting the Results %%%%%%%%%%%%%%%%%
% Velocity
figure(1); plot(sol.x,sol.y(1,:)); xlabel('Time [sec]'); ylabel('V [m/s]')
% Tang angle
figure(2); plot(sol.x,sol.y(2,:).*57.3); xlabel('Time [sec]'); ylabel('Tang [deg]')
% Height
figure(3); plot(sol.x,sol.y(3,:)); xlabel('Time [sec]'); ylabel('H [m]')
%%%%%%%%%%%%%%% Differential Equations %%%%%%%%%%%%%%%
function out = EOMsys(t,x)
t
v = x(1);  tang=x(2);   h=x(3);
[rho,T,P] = marsatmoshper(h);

%% Constants and Unknowns
g = 3.86;  % Mars Gravity Acceleration
% rho = 0.019;% Mars atmoshper density
R = 3390000;  % Mars Radius

% Lander Mass-Geometric Characteristics
s = 1.5174;  % Area of the Lander middle-section
L = 2;  % Lander Length?, 1.06
m = 576;  % Lander Mass
Ix = 270;  Iz = 443;
I = Iz; % Used For normalization
Cxa = 0.1;

%% Equations of Moving Center of Gravity
q = 0.5*rho*v^2;

dVdt = -Cxa*q*s/m - g*sin(tang);
dTangdt = -cos(tang)*(g-v^2/(h+R))/v;
dHdt = v*sin(tang);
out = [dVdt;dTangdt;dHdt];
end