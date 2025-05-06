function [dxdt]=DE_test(t,x) %% DE Function
%% Calling variabls
v = x(:,1); tang = x(:,2); h = x(:,3); omegax = x(:,4); theta = x(:,5); alpha = x(:,6);

%% Constants
g = 3.86; rho = 0.019; epsilon = 0.9; 
s = 1.5174; L = 1.06; m = 576; Ix = 270; Iz = 443; I = (Iz+Iz)/2; Ixd = Ix/I;

%% Without Pr 
% mad=0.0896; mxad=0.052; theta1=pi; theta2=0; % catching behaviour 
% mad=0.003; mxad=0.002; theta1=pi/2; theta2=0; % escabing behaviour 
%% With Pr
mad=0.089; mxad=0.05; theta1=pi; theta2=0; % catching behaviour 
% mad=0.003; mxad=0.002; theta1=pi/2; theta2=0; % escabing behaviour 
%% For Testing
% mad = 0.003; % Moves omegax on the vertical 
% mxad= 0.002; % Shapes the resonance frequency
% theta1 = pi/2; % pi/2
% theta2 = 0; % 0

%% Aerodynamics
[Cx, Cyn, Mz] = aero_file(alpha*57.3);
Cxa = Cx./alpha; 
mzn = Mz; 

%% Sublimentary equations
q = 0.5.*rho.*v.^2;
omega = sqrt(-mzn.*q.*s.*L.*cot(alpha)./I);
omegaa = sqrt(Ixd.^2.*omegax.^2./4 + omega.^2);
if omegax > 0, omega12 = Ixd.*omegax./2 + omegaa; else, omega12 = Ixd.*omegax./2 - omegaa; end
Fa = - (mzn./alpha).*q.*s.*L./I + omega12.^2./cos(alpha).^2 + (Ixd.*omegax-omega12).*(Ixd.*omegax-2.*omega12);

%% DE about CG
dvdt = -q.*Cxa.*s./m - g.*sin(tang);
dtangdt = -g.*cos(tang)./v;
dhdt = v.*sin(tang);

domegaxdt = -epsilon.*mxad.*omegax.^2.*sin(theta+theta2)./Ixd;
dthetadt = omegax - omega12;
if omegax > 0, dalphadt = epsilon.*((tan(alpha)./(2.*omegaa.^2)).*(omega./(2.*q)).*rho.*v.*dvdt + mad.*omegax.^2.*cos(theta+theta1)./(2.*omegaa))./(Fa./(4.*omegaa.^2)); end
if omegax < 0, dalphadt = epsilon.*((tan(alpha)./(2.*omegaa.^2)).*(omega./(2.*q)).*rho.*v.*dvdt - mad.*omegax.^2.*cos(theta+theta1)./(2.*omegaa))./(Fa./(4.*omegaa.^2)); end

% dxdt = [dvdt; dtangdt; dhdt; domegaxdt; dthetadt; dalphadt];
% disp(['Time = ' num2str(t) ', Alpha = ' num2str(alpha*57.3)])

%% Calculating the probability of catching or resonance
mxa1 = -mxad.*omegax.^2.*sin(theta2);
mxa2 = sqrt((mxad.*omegax.^2).^2 - mxa1.^2);
ma1 = mad.*omegax.^2.*sin(theta1);
ma2 = sqrt((mad.*omegax.^2).^2 - ma1.^2);
ma = mad.*omegax.^2;

% D = omegax - omega12; 
if omegax > 0, dbdalpha = -(2.*((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma1.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) - (L.*mxa2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*Ixd.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) - 2.*((mxa1.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma2.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) + (L.*mxa1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*Ixd.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))./(2.*(((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).^2 + ((mxa1.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).^2).^(1/2));
               dbdomega = (2.*((mxa1.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa1.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) + (L.*ma2.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) + 2.*((mxa2.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa2.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) - (L.*ma1.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))/(2.*(((mxa1.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2 + ((mxa2.*((Ixd.*omegax)./2 + ((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2).^(1/2));
               dDdalpha = -(L.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2));
               dDdomega = -omega./sqrt((Ixd.^2.*omegax.^2)/4 + omega.^2);
                    mf1 = -(1 - omega12).*(mxa2./Ixd) + dDdalpha.*(2.*omegaa.*ma1./Fa);
                    mf2 = (1 - omega12).*(mxa1./Ixd) + dDdalpha.*(2.*omegaa.*ma2./Fa);
         else, dbdalpha = -(2.*((mxa1.*(((Ixd.^2.*omegax^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma2.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) - (L.*mxa1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*Ixd.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) - 2.*((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).*((2.*L.*ma1.*mzn.*q.*s.*cot(alpha).*(cot(alpha).^2 + 1))./(Fa.*I) + (L.*mxa2.*mzn.*q.*s.*(cot(alpha).^2 + 1))/(2.*I.*Ixd.*((Ixd.^2*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))./(2.*(((mxa1.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd + (L.*ma2.*mzn.*q.*s.*(cot(alpha).^2 + 1))/(Fa.*I)).^2 + ((mxa2.*(((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2) + (Ixd.*omegax)./2 - 1))./Ixd - (L.*ma1.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I)).^2).^(1/2));
               dbdomega = (2.*((mxa2.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa2.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) + (L.*ma1.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))) + 2.*((mxa1.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).*((mxa1.*omega)./(Ixd.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2)) - (L.*ma2.*mzn.*omega.*q.*s.*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))))./(2.*(((mxa2.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd + (L.*ma1.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2 + ((mxa1.*(((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2) - (Ixd.*omegax)./2 + 1))./Ixd - (L.*ma2.*mzn.*q.*s.*((Ixd.^2.*omegax.^2)./4 + omega.^2).^(1/2).*(cot(alpha).^2 + 1))./(Fa.*I.*((Ixd.^2.*omegax.^2)./4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2))).^2).^(1/2));
               dDdalpha = (L.*mzn.*q.*s.*(cot(alpha).^2 + 1))./(2.*I.*((Ixd.^2.*omegax.^2)/4 - (L.*mzn.*q.*s.*cot(alpha))./I).^(1/2));
               dDdomega = omega./sqrt((Ixd.^2*omegax.^2)./4 + omega.^2);
                    mf1 =-(1 - omega12).*(mxa2./Ixd) - dDdalpha.*(2.*omegaa.*ma1./Fa);
                    mf2 = (1 - omega12).*(mxa1./Ixd) - dDdalpha.*(2.*omegaa.*ma2./Fa);
end
domegadt = -L.*mzn.*q.*s.*dalphadt.*(1+cot(alpha).^2)./(2.*I.*sqrt(-mzn.*q.*s.*L.*cot(alpha)./I));

b = sqrt(mf1.^2 + mf2.^2);
C1 = -dbdalpha.*(2.*omega.*tan(alpha).*domegadt./Fa) + dbdomega.*domegadt;
if omegax > 0; C2 = dbdalpha.*(2.*ma.*omegaa./Fa);
         else, C2 =-dbdalpha.*(2.*ma.*omegaa./Fa);end

theta3 = acos(mf1./b);
theta4 = theta1 - theta3;
Q = abs(24.*C1 + 8.*C2.*cos(theta4))./(3.*sqrt(b));
a = -dDdalpha.*(2.*omega.*tan(alpha).*domegadt./Fa) + dDdomega.*domegadt;

if Q < 4*pi*abs(a); Pr = Q./(2.*pi.*abs(a)+Q./2); else, Pr = 1; end

% Pr_save(i) = Pr
% assignin('base','Pr_save',Pr_save); % get the data to the workspace.

size(dvdt), size(dtangdt), size(dhdt),size(domegaxdt),size(dthetadt),size(dalphadt)     
dxdt = [Pr];
end