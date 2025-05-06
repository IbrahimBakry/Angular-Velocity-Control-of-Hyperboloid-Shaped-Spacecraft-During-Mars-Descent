%% The working version
clc; clear; close all

Pr_out = 0.0; T_out = 0; Thetad_out = 0; Theta_out = 0; Prp_out = 0;

%% Solving the DE
v0 = 3400;  tang0 = pi/6; h0 = 100000; 
omegax0 = 0.45;  theta0 = pi/180;  alpha0 = pi/20;  %% catching
% omegax0 = 0.6;  theta0 = 1/57.3;  alpha0 = 0.01/57.3;  %% catching
y0 = [v0, tang0, h0, omegax0, theta0, alpha0];
tspan = [0 325];

mad = 0.002; mxad = 0.003; theta1 = pi/2; % 1.1/0.2
% mad = 0.008; mxad = 0.008; theta1 = pi/2;

[x,y] = ode113(@(t,x)DE(t,x,mad,mxad,theta1),tspan,y0);
% i = 0;
% for mad = 0.001:0.001:0.009
%     for mxad = 0.001:0.001:0.009
%         [x,y] = ode113(@(t,x)DE(t,x,mad,mxad,theta1),tspan,y0);
%         
%         [Cx, Cy, mzn] = aero_file(y(:,6));
%         omega_xr = sqrt(-0.5.*0.019.*(y(:,1)).^2.*mzn.*1.5174*1.06.*cot((y(:,6)))./443)./sqrt(1-(270/443));
%         figure; plot(x,omega_xr,'--',x,y(:,4)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]'); title(['mad = ' num2str(mad) 'mxad = ' num2str(mxad)])
%         legend('OmegaxResonance','Omegax');  axis([tspan 0 4]); % For Catching
%         
%         i = i + 1;
%         disp(['fig_num = ' num2str(i) ', mad = ' num2str(mad) ', mxad = ' num2str(mxad)])
%     end
% end

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

% [Cx, Cy, mzn] = aero_file(y(:,6));
omega_xr = sqrt(-0.5.*0.019.*(y(:,1)).^2.*(-0.001).*1.5174*1.06.*cot((y(:,6)))./443)./sqrt(1-(270/443));
figure(3); plot(x,omega_xr,'--',x,y(:,4)); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
% figure(3); plot(x,omega_xr); xlabel('Time [sec]'); ylabel('Omegax [1/sec]')
legend('OmegaxResonance','Omegax'); 
% axis([tspan 0 4]); % For Catching
% Pr
% figure(4); plot(T_out,Pr_out); xlabel('Time [sec]'); ylabel('Pr [deg]')
axis([0 tspan(2) 0 2])

function [dxdt]=DE(t,x,mad,mxad,theta1) %% DE Function
%% Calling variabls
v = x(1); tang = x(2); h = x(3); omegax = x(4); theta = x(5); alpha = x(6);

%% Constants
g = 3.86; rho = 0.019; epsilon = 1;
s = 1.5174; L = 1.25; m = 576;  r = 3390000;
Ix = 270; Iz = 443; Iy = 443; Iyz = 50; deltaI = 0.5; 
I = (Iz+Iz)/2; 
Ixd = Ix/I;
deltaId = deltaI/I;
Iyzd = Iyz/I;

%% Without Pr 
% mad=0.0896; mxad=0.052; theta1=pi; theta2=0; % catching behaviour 
% mad=0.003; mxad=0.002; theta1=pi/2; theta2=0; % escabing behaviour 
%% With Pr
% mad=0.089; mxad=0.05; theta1=pi; theta2=0; % catching behaviour 
% mad=0.003; mxad=0.002; theta1=pi/2; theta2=0; % escabing behaviour 
%% For Testing
% mad = 0.01; % Moves omegax on the vertical 
% mxad= 0.06;  % Shapes the resonance frequency
% theta1 = pi; % 
theta2 = 0;  % 
mdelta = sqrt(Iyzd.^2+deltaId.^2);
theta3 = asin(deltaId/mdelta)/.2;

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
dtangdt = -(g-v^2/r).*cos(tang)./v;
dhdt = v.*sin(tang);

domegaxdt = epsilon.*(-mxad.*omegax.^2.*sin(theta+theta2));% + mdelta.*omega12.^2.*tan(alpha).^2.*cos(2*theta+2*theta3))./Ixd;
dthetadt = omegax - omega12;

PSI = epsilon.*pi.*rho.*v.*dvdt./(q.*omega);
if omegax > 0
%     dalphadt = (-PSI.*omega.^2.*tan(alpha)./(4.*omegaa.^2.*pi) + epsilon.*(0.5.*mad.*omegax.^2./omegaa).*cos(theta+theta1) - epsilon.*(0.25.*omega12.*tan(alpha)./omegaa.^2).*((10+Ixd).*omegax.*omega12-2.*(2+Ixd).*omegax.^2+(tan(alpha).^2-4).*omega12.^2).*mdelta.*cos(2.*theta+2.*theta3))./(Fa./(4.*omegaa.^2));   
     dalphadt = epsilon.*((tan(alpha)./(2.*omegaa.^2)).*(omega.^2./(2.*q)).*rho.*v.*dvdt + mad.*omega.^2.*cos(theta+theta1)./(2.*omegaa))./(Fa./(4.*omegaa.^2));
else          
%      dalphadt = (-PSI.*omega.^2.*tan(alpha)./(4.*omegaa.^2.*pi) - epsilon.*(0.5.*mad.*omegax.^2./omegaa).*cos(theta+theta1) - epsilon.*(0.25.*omega12.*tan(alpha)./omegaa.^2).*((10+Ixd).*omegax.*omega12-2.*(2+Ixd).*omegax.^2+(tan(alpha).^2-4).*omega12.^2).*mdelta.*cos(2.*theta+2.*theta3))./(Fa./(4.*omegaa.^2));   
     dalphadt = epsilon.*((tan(alpha)./(2.*omegaa.^2)).*(omega.^2./(2.*q)).*rho.*v.*dvdt - mad.*omega.^2.*cos(theta+theta1)./(2.*omegaa))./(Fa./(4.*omegaa.^2));
end
dxdt = [dvdt; dtangdt; dhdt; domegaxdt; dthetadt; dalphadt];
% disp(['Time = ' num2str(t) ', Alpha = ' num2str(alpha*57.3)])

% %% Calculating the probability of catching or resonance
% mxa1 = -mxad.*omegax.^2.*sin(theta2);
% mxa2 = sqrt((mxad.*omegax.^2).^2 - mxa1.^2);
% ma1 = mad.*omegax.^2.*sin(theta1);
% ma2 = sqrt((mad.*omegax.^2).^2 - ma1.^2);
% ma = mad.*omegax.^2;
% 
% if omegax > 0, dbdalpha = -(2.*((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma1.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) - (L.*mxa2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*Ixd.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) - 2.*((mxa1.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma2.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) + (L.*mxa1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*Ixd.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))./(2.*(((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).^2 + ((mxa1.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).^2).^(1/2));
%                dbdomega = (2.*((mxa1.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa1.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) + (L.*ma2.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) + 2.*((mxa2.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa2.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) - (L.*ma1.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))/(2.*(((mxa1.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2 + ((mxa2.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2).^(1/2));
%                dDdalpha = -(L.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2));
%                dDdomega = -omega./sqrt((Ixd.^2.*omegax.^2)/4 + omega.^2);
%                     mf1 = -(1 - omega12).*(mxa2./Ixd) + dDdalpha.*(2.*omegaa.*ma1./Fa);
%                     mf2 = (1 - omega12).*(mxa1./Ixd) + dDdalpha.*(2.*omegaa.*ma2./Fa);
%          else, dbdalpha = -(2.*((mxa1.*(((Ixd.^2.*omegax^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma2.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) - (L.*mxa1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*Ixd.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) - 2.*((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma1.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) + (L.*mxa2.*mzn.*q.*s.*(cot(alpha).^2 + 1))/(2.*I.*Ixd.*((Ixd.^2*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))./(2.*(((mxa1.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))/(Fa.*I)).^2 + ((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).^2).^(1/2));
%                dbdomega = (2.*((mxa2.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa2.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) + (L.*ma1.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) + 2.*((mxa1.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa1.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) - (L.*ma2.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))./(2.*(((mxa2.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))/(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2 + ((mxa1.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2).^(1/2));
%                dDdalpha = (L.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*((Ixd.^2.*omegax.^2)/4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2));
%                dDdomega = omega./sqrt((Ixd.^2*omegax.^2)./4 + omega.^2);
%                     mf1 =-(1 - omega12).*(mxa2./Ixd) - dDdalpha.*(2.*omegaa.*ma1./Fa);
%                     mf2 = (1 - omega12).*(mxa1./Ixd) - dDdalpha.*(2.*omegaa.*ma2./Fa);
% end
% domegadt = -L.*mzn.*q.*s.*dalphadt.*(1+cot(alpha).^2)./(2.*I.*sqrt(-mzn.*q.*s.*L.*cot(alpha)./I));
%      
% b = sqrt(mf1.^2 + mf2.^2);
% C1 = -dbdalpha.*(2.*omega.*tan(alpha).*domegadt./Fa) + dbdomega.*domegadt;
% if omegax > 0; C2 = dbdalpha.*(2.*ma.*omegaa./Fa);
%          else, C2 =-dbdalpha.*(2.*ma.*omegaa./Fa);end
% 
% theta3 = acos(mf1./b);
% theta4 = theta1 - theta3;
% Q = abs(24.*C1 + 8.*C2.*cos(theta4))./(3.*sqrt(b));
% a = -dDdalpha.*(2.*omega.*tan(alpha).*domegadt./Fa) + dDdomega.*domegadt;
% 
% % Prc - Catching Propability
% if Q < 4*pi*abs(a); Pr = Q./(2.*pi.*abs(a)+Q./2); else, Pr = 1; end
% % Prp - Twice Passing Propability
% Prp = (1 - Pr)*(1 - Pr); 
% 
% assignin('base','Pr_step',Pr); evalin('base','Pr_out(end+1) = Pr_step;');
% assignin('base','Prp_step',Prp); evalin('base','Prp_out(end+1) = Prp_step;');
% assignin('base','T_step',t); evalin('base','T_out(end+1) = T_step;');
% assignin('base','Thetad_step',dthetadt); evalin('base','Thetad_out(end+1) = Thetad_step;');
% assignin('base','Theta_step',theta); evalin('base','Theta_out(end+1) = Theta_step;');
% 
% dxdt = [dvdt; dtangdt; dhdt; domegaxdt; dthetadt; dalphadt];
% disp(['Time = ' num2str(t) ', Alpha = ' num2str(alpha*57.3) ', Pr = ' num2str(Pr)])
end