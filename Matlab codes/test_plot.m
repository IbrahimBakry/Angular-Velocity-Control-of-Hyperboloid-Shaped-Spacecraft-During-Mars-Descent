clc;
clear;
close all;

% test plot
A = 0:180;
% Mzfit = -1.2*sind(A) - 0.8*sind(2.*A);
% Cyfit = 0.83*sind(2.*A);
% Cxfit = 1.58*cosd(4.*A);

alpha = [0 5 10 15 20 90];
Cy = [0 0.15 0.28 0.42 0.53 0];
Cx = [1.58 1.56 1.5 1.47 1.4 1.58];
Mz = [0 -0.02 -0.04 -0.05 -0.06 0];

Cypoly = polyfit(alpha,Cy,2);
Cyfit = polyval(Cypoly,A*57.3);
figure; plot(alpha,Cy,'o',[0:90],polyval(Cypoly,[0:90]));

Cxpoly = polyfit(alpha,Cx,2);
Cxfit = polyval(Cxpoly,A*57.3);
figure; plot(alpha,Cx,'o',[0:90],polyval(Cxpoly,[0:90]));

Cmpoly = polyfit(alpha,Mz,2);
Mzfit = polyval(Cmpoly,A*57.3);
figure; plot(alpha,Mz,'o',[0:90],polyval(Cmpoly,[0:90]));


% plot(A,Mzfit,alpha,Mz,'o'); title('Cx(alpha)'); ylabel('Cx'); xlabel('alpha')

% Mzfit = -0.06*sind(4.*A);
% Cyfit = 0.83*sind(2.*A);
% Cxfit = 0.24*cosd(4.*A);

% subplot(311);plot(A,Mzfit,alpha,Mz,'o'); title('Mz(alpha)'); ylabel('Mz')
% subplot(312);plot(A,Cyfit,alpha,Cy,'o'); title('Cy(alpha)'); ylabel('Cy')
% subplot(313);plot(A,Cxfit,alpha,Cx,'o'); title('Cx(alpha)'); ylabel('Cx'); xlabel('alpha')
