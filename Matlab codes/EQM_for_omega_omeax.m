%% The working version
clc; clear; close all

% Solving the DE
omegax0 = 0.51;  theta0 = 1/57.3;  alpha0 = 0.01/57.3;
y0 = [omegax0];
tf = 335;
tspan = [0 tf];

[x1,y1]=Integrating_the_CG_equations_of_motion(tf);
vpoly = polyfit(x1,y1(:,1),2);

[x,y] = ode45(@(t,x)DE(t,x,vpoly),tspan,y0);

% Omegax, Theta, Alpha
figure(2)
% subplot(311); 
plot(x,y(:,1)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
% subplot(312); plot(x,y(:,2).*57.3); xlabel('Time [sec]'); ylabel('Theta [deg]')
% subplot(313); plot(x,y(:,3).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')

% [Cx, Cy, mzn] = aero_file(y(:,6));
% omega_xr = sqrt(-0.5.*0.019.*(y(:,1)).^2.*mzn.*1.5174*1.06.*cot((y(:,6)))./443)./sqrt(1-(270/443));
% % rho = 0.019; s = 1.5174; L = 2; Ix = 270; Iz = 443; I = (Iz+Iz)/2; Ixd = Ix/I;
% % omega_xr = sqrt(-mzn.*(0.5.*rho.*y(:,1).^2).*s.*L./I)/sqrt(1-Ixd);
% figure(3); plot(x,omega_xr,'--',x,y(:,4)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
% legend('OmegaxResonance','Omegax'); 

% axis([tspan 0 4]); % For Catching
% axis([0 tspan(2) 0 2])

function [dxdt]=DE(t,x,vpoly) %% DE Function
%% Calling variabls
omegax = x(1);

%% Constants
g = 3.86; rho = 0.019; epsilon = 1;
s = 1.5174; L = 2; m = 576; Ix = 270; Iz = 443; I = (Iz+Iz)/2; Ixd = Ix/I;

% Aysmetry Parameters
mxw = 0.01;

v = polyval(vpoly,t*57.3);

%% Sublimentary equations
q = 0.5.*rho.*v.^2;

% Diff. Eq. about CG
dvdt = -q.*Cxa.*s./m - g.*sin(tang);
dtangdt = -g.*cos(tang)./v;
domegaxdt = epsilon.*(mxw*q*s^2*L).*omegax./(Ixd*v);

dxdt = [domegaxdt];
disp(['Time = ' num2str(t)])
end
