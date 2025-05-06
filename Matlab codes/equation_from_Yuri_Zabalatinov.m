%%%
clc; close all; clear

v0 = 7000; tang0 = -10/57.3; h0 = 110000;
a10 = 1/57.3; omegax0 = 10; theta0 = -10/57.3;
% omega0 = 10; 

y0 = [omegax0, theta0, a10, v0, tang0];
% y0 = [omegax0, theta0, a10];
tf = 400;
tspan = [0 tf];
% vpoly = polyfit([0:5:tf],v0.*exp(-0.02.*[0:5:tf]),2);

[x,y] = ode23(@(t,x)DE(t,x),tspan,y0);

% Omegax, Theta, Alpha
figure(1)
subplot(211); plot(x,y(:,1)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
subplot(212); plot(x,y(:,3).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')
figure(3)
subplot(211); plot(x,y(:,4)); xlabel('Time [sec]'); ylabel('V [m/sec]')
subplot(212); plot(x,y(:,5).*57.3); xlabel('Time [sec]'); ylabel('tang [deg]')

s = 1; L = 1.5; Ixd = 0.5; I = 20; mz1 = -0.05; rho = 0.019;

a1 = y(:,3); v = y(:,4);
% omega_xr = sqrt(-0.5*mz1.*rho.*v.^2.*s.*L)/sqrt((1-Ixd)*I);
omega_xr = -0.0027.*v.*sqrt(a1)+11.2;
% omega_xr = omega/sqrt(1-Ixd);

figure(2); plot(x,omega_xr,'--',x,y(:,1)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
legend('OmegaxResonance','Omegax'); 

function [dxdt]=DE(t,x)
omegax = x(1);
theta = x(2);
a1 = x(3);
v = x(4);
tang = x(5);

L = 1.5; Ixd = 0.5; I = 20; mz1 = -0.05; rho = 0.019; s = 1; g = 9.81;
ma = 0.002; mxa = 0.003; theta1 = pi/2; theta2 = 0; m = 1000; cx1 = 1.1;
thetads0 = pi/10; r = 3390000;

% v = polyval(vpoly,t);

cxa = cx1./cos(a1);
q = 0.5.*rho.*v.^2;
omega = sqrt(-mz1.*q.*s.*L.*a1./I);
omegaa = sqrt(Ixd.^2.*omegax.^2./4 + omega.^2);
omega1 = Ixd.*omegax./2 + omegaa;

% t = 0:5:400; y = 3500.*exp(-0.02.*t); plot(t,y)
dvdt = -q.*cxa.*s./m  - g.*sin(tang);
dtangdt = -(g-v^2/r).*cos(tang)./v;

% da1dt = -omega*a1*(omega*(rho*v*dvdt)/(2*q))/(2*omegaa^2) - ma*cos(theta+theta1)/(2*omegaa);
% da1dt = -omega*a1*(omega*(rho*v)/(2*q))/(2*omegaa^2) - ma*cos(theta+theta1)/(2*omegaa);
% da1dt = 2*(1-Ixd)*a1*dwdt/(omega*(2-Ixd)^2) - mad*omega*sqrt(1-Ixd)*cos(thetads0 + theta1-theta2);
da1dt = (1-Ixd)*a1.*((rho.*v.*dvdt)/(2.*q))./(omega*(2-Ixd)^2) - 5*(ma/omega^2)*omega*sqrt(1-Ixd)*cos(thetads0 + theta1-theta2);
domegaxdt = -a1*mxa*sin(theta+theta2)/Ixd;
% dthetadt = (1-Ixd)*omegax;
% dthetadt = omegax-omega1;
dthetadt = omegax-omega1 + ma*sin(theta+theta1)/(2*a1*omegaa);
% domegadt = q/(2*omega);

dxdt = [domegaxdt; dthetadt; da1dt; dvdt; dtangdt];
% dxdt = [domegaxdt; dthetadt; da1dt];
t
a1
theta 
omegax
end
