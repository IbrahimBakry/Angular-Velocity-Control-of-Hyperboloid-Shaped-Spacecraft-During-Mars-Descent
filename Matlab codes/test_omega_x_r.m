clc; clear; close all
v0 = 7000;
tf = 500;
t = 0:1:tf;
s = 1; L = 1.5; Ixd = 0.5; I = 20; mz1 = -0.5; rho = 0.019;

vpoly = polyfit(t,v0.*exp(0.01.*t),2);
v = polyval(vpoly,t);

omega_xr = sqrt(-0.5*mz1.*rho.*v.^2.*s.*L.*0.1)/sqrt((1-Ixd)*I); %
figure; plot([0:1:tf],omega_xr); xlabel('Time [sec]'); 
