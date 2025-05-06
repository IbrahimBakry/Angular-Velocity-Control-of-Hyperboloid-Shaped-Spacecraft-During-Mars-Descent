% function [x,y]=Integrating_the_CG_equations_of_motion(tf)
%% Integrating the CG equations of motion
% Solving the DE
v0 = 3500; tang0 = -10/57.3; h0 = 120000; 
[x,y] = ode23(@(t,x)DE(t,x),[0 350],[v0, tang0, h0]);

% Velocity, Tang angle, Height
figure(1); 
subplot(311); plot(x,y(:,1)); xlabel('Time [sec]'); ylabel('V [m/s]')
subplot(312); plot(x,y(:,2).*57.3); xlabel('Time [sec]'); ylabel('Tang [deg]')
subplot(313); plot(x,y(:,3)); xlabel('Time [sec]'); ylabel('H [m]')

function [dxdt]=DE(t,x) %% DE Function
%% Calling variabls
v = x(1); tang = x(2); h = x(3);

%% Constants
g = 3.86; rho = 0.019; s = 1.5174; m = 576;

%% Aerodynamics
Cxa = 0.95967; 

%% Sublimentary equations
q = 0.5.*rho.*v.^2;

%% Diff. Eq. of CG
dvdt = -q.*Cxa.*s./m - g.*sin(tang);
dtangdt = -g.*cos(tang)./v;
dhdt = v.*sin(tang);

%% function output
dxdt = [dvdt; dtangdt; dhdt];
end
% end