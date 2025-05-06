clc; clear; close all

alpha = [0 5 10 15 20];
Cy = [0 0.15 0.28 0.42 0.53] ;
Cx = [0.24 0.22 0.185 0.13 0.04];
Mz = [-0.001 -0.02 -0.04 -0.05 -0.06];

A = 0:90;
Mzfit = -0.06*sind(4.*A);% + Mz2*cosd(A);
Cyfit = 0.83*sind(2.*A);
Cxfit = 0.24*cosd(4.*A);

plot(alpha,Mz,A,Mzfit)

