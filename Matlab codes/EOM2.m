%% Solving the Equation of Motion of the Lander with a significance axisymetry
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%% DE Integration %%%%%%%%%%%%%%%%%%
% Initial Conditions
H0 = 100000;     % Height
tang0 = -1/57.3;   % Trajectory Angle
v0 = 3400;       % Speed
alpha0 = 5/57.3;  % Angle of Attack
omegax0 = 25/57.3;  % Angular velocity of the Lander relative to OX axis;
theta0 = 30/57.3;% Fast phase -pi/180
% omega0 = 1/57.3;    % precession frequency with angular velocity ?x
% q0 = 0.1*109820;
% y0 = [v0, tang0, H0, omegax0, theta0, alpha0, omega0];
% y0 = [v0, tang0, H0, omegax0, theta0, alpha0, q0];
y0 = [v0, tang0, H0, omegax0, theta0, alpha0];

% Calling the Solver
tspan = [0 230]; %93
sol = ode45(@EOMsys,tspan,y0);

% Calculating omega_x_r
v = sol.y(1,:); alpha = sol.y(6,:); %q = sol.y(7,:);
rho = 0.019; s = 1.5174; L = 1.06;
Ix = 270;  Iz = 443; I = (Iz+Iz)/2;

q = 0.5.*rho.*v.^2;
omega_xr = sqrt(0.2.*q.*s.*L.*cot(alpha)./I)./sqrt(1-(Ix/I));

%%%%%%%%%%%%%%%% Chartting the Results %%%%%%%%%%%%%%%%%
% Velocity
figure(1); plot(sol.x,sol.y(1,:)); xlabel('Time [sec]'); ylabel('V [m/s]')
% Tang angle
figure(2); plot(sol.x,sol.y(2,:).*57.3); xlabel('Time [sec]'); ylabel('Tang [deg]')
% Height
figure(3); plot(sol.x,sol.y(3,:)); xlabel('Time [sec]'); ylabel('H [m]')
% % Omegax
% figure(4); plot(sol.x,sol.y(4,:).*57.3); xlabel('Time [sec]'); ylabel('Omegax [deg/sec]')
% % Theta
% figure(5); plot(sol.x,sol.y(5,:).*57.3); xlabel('Time [sec]'); ylabel('Theta [deg]')
% % Alpha
% figure(6); plot(sol.x,sol.y(6,:).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')
% Omegax, Theta, Alpha
figure(7)
subplot(311); plot(sol.x,sol.y(4,:)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
subplot(312); plot(sol.x,sol.y(5,:).*57.3); xlabel('Time [sec]'); ylabel('Theta [deg]')
subplot(313); plot(sol.x,sol.y(6,:).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')
% Omega_xr
figure(8); plot(sol.x,omega_xr,'--'); hold on
           plot(sol.x,sol.y(4,:)); hold off
xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
legend('OmegaxResonance','Omegax')


%%%%%%%%%%%%%%% Differential Equations Function %%%%%%%%%%%%%%%
function out = EOMsys(t,x)

v = x(1);      tang=x(2);   h=x(3);
omegax = x(4); theta=x(5);  alpha=x(6);
% q = x(7)
% omega = x(7);

%% Constants and Unknowns
g = 3.86;  % Mars Gravity Acceleration
rho = 0.019;% Mars atmoshper density
R = 3390000;  % Mars Radius

% Lander Mass-Geometric Characteristics
s = 1.5174;  % Area of the Lander middle-section
L = 2;  % Lander Length?, 1.06
m = 576;  % Lander Mass
Ix = 270;  
Iz = 443;
Ixz= 0.0;  
Ixy= 0.0;
I = (Iz+Iz)/2; % Used For normalization

% aerodynamic File
[Cx, Cy, Mz] = aero_file(alpha*57.3);
Cxa = Cx/alpha; % coefficient of Drag aerodynamic force D in Speed coordinates
Cya = Cy/alpha;
Cyn = Cy; % coefficient of Lift aeridynamic force L in body with axicemetry
mzn = Mz; % Cx: coefficient of Drag aerodynamic force D in Body coordinates

% Aerodynamic moments with assymetry
% mzn = -0.5*rho*v^2*s*L*(2*sin(alpha));%-0.3*sin(2*alpha));
% Mzn = Mz; % Coefficient of recovered aerodynamic moment 
mzf = 0.5; % Coefficients of small aerodynamic coefficients based on the asymmetry of KA
myf = 0.5; % Coefficients of small aerodynamic lift based on the asymmetry of KA
dz = 0.4; % Z displacement of CG in OXYZ, dz/L > 0.3
dy = -0.8; % Y displacement of CG in OXYZ, dy/L > 0.3
epsilon = 0.1; % Parameter related to the amount of asymmetry of mass, inertia and aerodynamic

mxsf = 0; % mxf: moment of shape asymmetry0
mxcf = 0; % mxf = mxsf*sin(phi) + mxcf*sin(phi);

%% Supporting Equations
q = 0.5*rho*v^2;
dqdt = -rho*v*(Cxa*q*s/m + g*sin(tang)); % dvdt = - (Cxa*q*s/m + g*sin(tang))
omega= sqrt(-mzn*q*s*L*cot(alpha)/I);
omegaa = sqrt((Ix*omegax/I)^2/4 + omega^2);

if omegax > 0, omega12 = ((Ix/I)*omegax)/2 + omegaa; end 
if omegax <= 0,omega12 = ((Ix/I)*omegax)/2 - omegaa; end

m1a = -(((1+Ix/I)*omegax-3*omega12)/(2*omegaa))*(omega^2/mzn)*(myf-Cx*dz/L)*tan(alpha)-omega12*omega^2*tan(alpha)^2*(mxcf+Cyn*dz/L)/(2*omegaa*mzn); 
m2a = -(((1+Ix/I)*omegax-3*omega12)/(2*omegaa))*(omega^2/mzn)*(mzf+Cx*dy/L)*tan(alpha)+omega12*omega^2*tan(alpha)^2*(mxsf+Cyn*dy/L)/(2*omegaa*mzn); 
ma = sqrt(m1a^2 + m2a^2);
theta1 = asin(m1a/ma); 
% theta1 = acos(-m2a/ma); % theta1 = atand(-m1a/m2a); % ma = 0.002; theta1 = pi;

mx1a = -(omega^2/mzn)*(mxsf+Cyn*dy/L)*tan(alpha) - (Ixy/I)*omega12^2*tan(alpha);
mx2a = -(omega^2/mzn)*(mxcf+Cyn*dz/L)*tan(alpha) - (Ixz/I)*omega12^2*tan(alpha);
mxa = sqrt(mx1a^2 + mx2a^2);
theta2 = asin(-mx1a/mxa); 
% theta2 = acos(mx2a/mxa); % theta2 = atand(-mx1a/mx2a); % mxa = 0.003; theta2 = 0;

Fa = - mzn/(alpha*I) + omega12^2/cos(alpha)^2 + (omegax*Ix/I - omega12)*(omegax*Ix/I - 2*omega12);

disp(['Time = ' num2str(t) ', Alpha = ' num2str(alpha*57.3) ', T1 =' num2str(theta1*57.3) ', T2 =' num2str(theta2*57.3)])

%% Equations of Motion of CG of Lander
dVdt = -Cxa*q*s/m - g*sin(tang);
dTangdt = -cos(tang)*g/v; % -cos(tang)*(g-v^2/(h+R))/v
dHdt = v*sin(tang);

%% Equations of Movtion Around CG of Lander
dOmegaxdt = -epsilon*mxa*sin(theta+theta2)/(Ix/I); % First eq. of OmegaX
dThetadt = omegax - omega12; % Second Eq. of Theta 

if omegax > 0, dAlphadt = -epsilon*(tan(alpha)*omega^2*dqdt/(4*omegaa^2*q) - epsilon*ma*cos(theta+theta1)/(2*omegaa))/(Fa/(4*omegaa^2)); end % Third Eq. of Alpha
if omegax <= 0,dAlphadt = -epsilon*(tan(alpha)*omega^2*dqdt/(4*omegaa^2*q) + epsilon*ma*cos(theta+theta1)/(2*omegaa))/(Fa/(4*omegaa^2)); end
% dOmegadt = epsilon*omega*dqdt/(2*q); % Fourth Eq of Omega

% out = [dVdt;dTangdt;dHdt;dOmegaxdt;dThetadt;dAlphadt;dOmegadt];
% out = [dVdt;dTangdt;dHdt;dOmegaxdt;dThetadt;dAlphadt;dqdt];
out = [dVdt;dTangdt;dHdt;dOmegaxdt;dThetadt;dAlphadt];
end
