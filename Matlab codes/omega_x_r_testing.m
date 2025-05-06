% NOT WORKING
clc;clear; close all

t = -2:0.01:2;
alpha = -t.^2;
mz = -0.02;
rho = 0.019; s = 1.5174; L = 2; 
Ix = 270; Iz = 443; I = (Iz+Iz)/2; Ixd = Ix/I;
v = 240*exp(-t);
q = 0.5.*rho.*v.^2;
omegaxr = sqrt(-mz.*q.*s.*L.*cot(alpha)./(I.*(1-Ixd)));

plot(t,omegaxr)