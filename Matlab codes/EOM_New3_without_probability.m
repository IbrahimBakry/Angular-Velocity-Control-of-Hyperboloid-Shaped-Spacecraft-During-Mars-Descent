%% The working version
clc; clear; close all

%% Solving the DE
v0 = 3400; tang0 = -1/57.3; h0 = 100000; 
omegax0 = 5;  theta0 = 1/57.3;  alpha0 = 0.01/57.3;  %% catching

y0 = [v0, tang0, h0, omegax0, theta0, alpha0];
tspan = [0 400];

mad = 0.4; mxad = 0.2; theta1 = pi; % catching behaviour
% mad = 0.008; mxad = 0.008; theta1 = pi/2;

[x,y] = ode113(@(t,x)DE(t,x,mad,mxad,theta1),tspan,y0);

% Velocity, Tang angle, Height
figure(1); 
subplot(311); plot(x,y(:,1)); xlabel('Time [sec]'); ylabel('V [m/s]')
subplot(312); plot(x,y(:,2).*57.3); xlabel('Time [sec]'); ylabel('Tang [deg]')
subplot(313); plot(x,y(:,3)); xlabel('Time [sec]'); ylabel('H [m]')
% Omegax, Theta, Alpha
figure(2)
subplot(311); plot(x,y(:,4)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
subplot(312); plot(x,y(:,5).*57.3); xlabel('Time [sec]'); ylabel('Theta [deg]')
subplot(313); plot(x,y(:,6).*57.3); xlabel('Time [sec]'); ylabel('Alpha [deg]')

[Cx, Cy, mzn] = aero_file(y(:,6));
% K = -0.1; omega_xr = (2*K+1).*sqrt(-0.5.*0.019.*(y(:,1)).^2.*mzn.*1.5174*1.06.*cot((y(:,6)))./443)./sqrt(1-(270/443)-(K+1)*K*(270/443)^2);
omega_xr = sqrt(-0.5.*0.019.*(y(:,1)).^2.*mzn.*1.5174*1.06.*cot((y(:,6)))./443)./sqrt(1-(270/443));

figure(3); plot(x,omega_xr,'--',x,y(:,4)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
% figure(3); plot(x,omega_xr); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
legend('OmegaxResonance','Omegax'); 
% axis([0 tspan(2) 0 2])

function [dxdt]=DE(t,x,mad,mxad,theta1) %% DE Function
%% Calling variabls
v = x(1); tang = x(2); h = x(3); 
omegax = x(4); theta = x(5); alpha = x(6);

%% Constants
g = 3.86; rho = 0.019; epsilon = 1; r = 3396.2*10^3;
s = 1.5174; L = 2; m = 576; 
Ix = 270; Iz = 443; Iy = 443; Iyz = 100; deltaI = 0; 
I = (Iz+Iz)/2; 
Ixd = Ix/I;
deltaId = deltaI/I;
Iyzd = Iyz/I;

theta2 = 0;  % 
theta3 = 0;
mdelta = sqrt(Iyzd.^2+deltaId.^2);
% theta3 = asin(deltaId/mdelta)/.2;

%% Aerodynamics
[Cxa, Cyn, mzn] = aero_file(alpha);
% Cxa = Cx.*cot(alpha); 

%% Sublimentary equations
q = 0.5.*rho.*v.^2;
omega = sqrt(-mzn.*q.*s.*L.*cot(alpha)./I);
omegaa = sqrt(Ixd.^2.*omegax.^2./4 + omega.^2);
if omegax > 0, omega12 = Ixd.*omegax./2 + omegaa; 
else,          omega12 = Ixd.*omegax./2 - omegaa; 
end
Fa = - (mzn.*cot(alpha)).*q.*s.*L./I + omega12.^2./cos(alpha).^2 + (Ixd.*omegax-omega12).*(Ixd.*omegax-2.*omega12);
% Fa = + omega12.^2./cos(alpha).^2 + (Ixd.*omegax-omega12).*(Ixd.*omegax-2.*omega12);

%% Diff. Eq. about CG
dvdt = -q.*Cxa.*s./m - g.*sin(tang);
% dtangdt = -g.*cos(tang)./v;
dtangdt = -cos(tang)*(g - v^2/(h-r))/v;
dhdt = v.*sin(tang);

domegaxdt = epsilon.*(-mxad.*omegax.^2.*sin(theta+theta2) + mdelta.*omega12.^2.*tan(alpha).^2.*cos(2*theta+2*theta3))./Ixd;
dthetadt = omegax - omega12;

PSI = epsilon.*pi.*rho.*v.*dvdt./(q.*omega);
if omegax > 0
%     dalphadt = (-PSI.*omega.^2.*tan(alpha)./(4.*omegaa.^2.*pi) + epsilon.*(0.5.*mad.*omegax.^2./omegaa).*cos(theta+theta1) - epsilon.*(0.25.*omega12.*tan(alpha)./omegaa.^2).*((10+Ixd).*omegax.*omega12-2.*(2+Ixd).*omegax.^2+(tan(alpha).^2-4).*omega12.^2).*mdelta.*cos(2.*theta+2.*theta3))./(Fa./(4.*omegaa.^2));   
    dalphadt = -epsilon.*(epsilon.*(tan(alpha)./(2.*omegaa.^2)).*(omega./(2.*q)).*rho.*v.*dvdt + mad.*omegax.^2.*cos(theta+theta1)./(2.*omegaa))./(Fa./(4.*omegaa.^2));
else        
%     dalphadt = (-PSI.*omega.^2.*tan(alpha)./(4.*omegaa.^2.*pi) - epsilon.*(0.5.*mad.*omegax.^2./omegaa).*cos(theta+theta1) - epsilon.*(0.25.*omega12.*tan(alpha)./omegaa.^2).*((10+Ixd).*omegax.*omega12-2.*(2+Ixd).*omegax.^2+(tan(alpha).^2-4).*omega12.^2).*mdelta.*cos(2.*theta+2.*theta3))./(Fa./(4.*omegaa.^2));   
    dalphadt = -epsilon.*(epsilon.*(tan(alpha)./(2.*omegaa.^2)).*(omega./(2.*q)).*rho.*v.*dvdt - mad.*omegax.^2.*cos(theta+theta1)./(2.*omegaa))./(Fa./(4.*omegaa.^2));
end
dxdt = [dvdt; dtangdt; dhdt; domegaxdt; dthetadt; dalphadt];
disp(['Time = ' num2str(t) ', Alpha = ' num2str(alpha*57.3)])
end