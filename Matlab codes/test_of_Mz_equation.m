clc; clear; close all

alpha = [0 5 10 15 20];
Mz = [0 -0.02 -0.04 -0.05 -0.06];
Mzs1 = 0.723417;
Mzs2 = -0.4782656;
Mzc1 = -0.269692;
Mzc2 = 0.25250197;
Mzseq = Mzs1.*sind([0:0.1:90]) + Mzs2.*sind(2.*[0:0.1:90]);
% Mzceq = Mzc1.*cosd([0:0.1:90]) + Mzc2.*cosd(2.*[0:0.1:90]);
% syms Mz1 Mz2
% [M1,M2]=solve([-0.02==Mz1*sind(5)+Mz2*sind(10),-0.06==Mz1*sind(20)+Mz2*sind(40)],[Mz1,Mz2])

Cmpoly = polyfit(alpha,Mz,2);
Mfitted = polyval(Cmpoly,[0:0.1:90]);
figure(3); plot(alpha,Mz,'o',[0:0.1:90],Mfitted,[0:0.1:90],Mzseq);
legend('Given','X^2','sin')
