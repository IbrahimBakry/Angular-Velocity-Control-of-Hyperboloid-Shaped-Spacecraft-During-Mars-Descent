
clc; close all; clear

v0 = 7000; a10 = 10/57.3; omegax0 = 10; theta0 = -10/57.3; %a10 = 1.57079632679;
y0 = [omegax0, theta0, a10];
tf = 500;
tspan = [0 tf];
vpoly = polyfit([0:5:tf],v0.*exp(-0.04.*[0:5:tf]),2);
% vpoly = polyfit([0:0.005:0.955],7000.*[0:0.01:1 1:-0.01:0.1].^2,2);
ma = 0.05; mxa = 0.1; theta1 = pi/2; theta2 = 0;

[x,y] = ode45(@(t,x)DE(t,x,vpoly,ma,mxa,theta1,theta2),tspan,y0);

s = 1; L = 1.5; Ixd = 0.5; I = 20; mz1 = -0.05; rho = 0.019; %mz1=-0.05
a1 = y(:,3); v = polyval(vpoly,x);
omega_xr = sqrt(-0.5*mz1.*rho.*v.^2.*s.*L.*a1)/sqrt((1-Ixd)*I); %
% omega_xr = 0.0398.*v0.*exp(-0.005.*[0:5:tf]);
% omega_xr = 4*sqrt(-0.5*mz1.*rho.*v.^2.*s.*L)/sqrt((4-4*Ixd-3*Ixd^2)*I);

% Omegax, Theta, Alpha
figure(1); plot(x,y(:,1)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
figure(2);
% figure(3); plot(x,y(:,2)); xlabel('Time [sec]'); ylabel('Theta')
subplot(211); plot(x,omega_xr,'--',x,y(:,1)); xlabel('Time [sec]'); %,'--',x,y(:,1) ylabel('Omegax [1/sec]'); legend('OmegaxResonance','Omegax'); 
subplot(212); plot(x,y(:,3).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')

function [dxdt]=DE(t,x,vpoly,ma,mxa,theta1,theta2)
omegax = x(1); theta = x(2); a1 = x(3);

L = 1.5; Ixd = 0.5; I = 20; mz1 = -0.05; rho = 0.019; s = 1;

v = polyval(vpoly,t);

q = 0.5.*rho.*v.^2;
omega = sqrt(-mz1.*q.*s.*L.*a1./I);
omegaa = sqrt(Ixd.^2.*omegax.^2./4 + omega.^2);
omega1 = Ixd.*omegax./2 + omegaa;

% da1dt = -omega*a1*(omega*(rho*v*dvdt)/(2*q))/(2*omegaa^2) - ma*cos(theta+theta1)/(2*omegaa);
da1dt = -1.*omega*a1*(omega*(rho*v)/(2*q))/(2*omegaa^2) - ma*cos(theta+theta1)/(2*omegaa);
domegaxdt = -10.*a1*mxa*sin(theta+theta2)/Ixd;
dthetadt = (1-Ixd)*omegax;
% dthetadt = omegax-omega1 + ma*sin(theta+theta1)/(2*a1*omegaa);

dxdt = [domegaxdt; dthetadt; da1dt];
end
