function [Cxfit, Cyfit, Mzfit]=aero_file(A)
%% Aerodynamic File

% Mzfit = -0.6*sin(A);% - 0.4*sind(2.*A); % lubimov book
% Cyfit = 2*sin(A); 
% Cxfit = -1*cos(A); 

% Mzfit = -0.2*sin(A);% + Mz2*sind(2*A); %-0.2
% Cyfit = 1.6*sin(A); % 1.6
% Cxfit = 1.58*cos(A); %0.24

Mzfit = -0.06*sin(4.*A);% + Mz2*sind(2*A);
Cyfit = 0.83*sin(2.*A);
Cxfit = 1.58*cos(4.*A);

% alpha = [0 5 10 15 20 90];
% Cy = [0 0.15 0.28 0.42 0.53 0];
% Cx = [1.58 1.56 1.5 1.47 1.4 1.58];
% Mz = [0 -0.02 -0.04 -0.05 -0.06 0];

% % Mz1 = 0.723417; Mz2 = -0.4782656;
% % Mzseq = Mz1*sind(A) + Mz2*sind(2*A);
% % syms Mz1 Mz2
% % [M1,M2]=solve([-0.02==Mz1*sind(5)+Mz2*sind(10),-0.06==Mz1*sind(20)+Mz2*sind(40)],[Mz1,Mz2])

% % Cypoly = polyfit(alpha,Cy,2);
% % Cypoly = [-0.002104477, 0.27910447]; % First order
% Cypoly = [-0.00037755, 0.0340727, -0.0081435]; % Second order
% % Cypoly = [-0.00000226, -0.000130697, 0.03008568, 0.00040468];
% Cyfit = polyval(Cypoly,A*57.3);
% % figure; plot(alpha,Cy,'o',[0:90],polyval(Cypoly,[0:90]));
% 
% % Cxpoly = polyfit(alpha,Cx,2);
% Cxpoly = [0.000528, 1.5026716];  % First order
% Cxpoly = [0.00012378, -0.0113329, 1.59685];  % Second order
% % Cxpoly = [0.00000407677, -0.00032098, -0.004149, 1.5814489];
% Cxfit = polyval(Cxpoly,A*57.3);  % First order
% % figure; plot(alpha,Cx,'o',[0:90],polyval(Cxpoly,[0:90]));
% 
% % Cmpoly = polyfit(alpha,Mz,2);
% % Cmpoly = [0.0002716, -0.0346716];
% Cmpoly = [0.000043395, -0.00388, -0.00165];  % Second order
% % Cmpoly = [-0.0000005328, 0.0001015269, -0.004825, 0.00035774];
% Mzfit = polyval(Cmpoly,A*57.3);
% % figure; plot(alpha,Mz,'o',[0:90],polyval(Cmpoly,[0:90]));
% 



end