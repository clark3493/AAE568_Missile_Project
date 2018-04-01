%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script Name:            Simple_2D_TPBVP.m
%   Script Author:          Steve Clark
%   Modified By:            Reese Robinson
%
%   SCRIPT DESCRIPTION
%       This script solves a simplified TP-BVP for the missile project.
%       - Only 2D dynamics are considered
%       - Assumes that the missile is well controlled by arbitrary control 
%         surface inputs such that the missile's pitch Euler angle can be 
%         considered a control and control surface deflections are not 
%         considered.
%       - No inequality constraints (i.e. no limits on vehicle or control
%         surface angle of attack)
%       - Constant thrust
%       - Aerodynamic forces are not considered
%
%       ** SOLUTION IS THE SAME AS ORIGINAL SIMPLE_2D_TPBVP EXCEPT IMPACT
%       ANGLE IS NOT A CONSTRAINT
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all;

%% ----- USER INPUTS ------------------------------------------------------

% initial conditions
global x0 z0 xdot0 zdot0;
x0 = 0.;        % x coord, absolute coordinates, meters
z0 = -4000.;    % z coord, absolute coordinates, meters  (z = -altitude)
xdot0 = 500.;   % x velocity, absolute coordinates, m/s
zdot0 = 0.;     % z velocity, absolute coordinates, m/s

% final conditions
global X Z eta;
X = 10000;       % target x absolute coordinates, meters
Z = 0.;         % target z absolute coordinates, meters

% constant parameters
global m T g;
m = 100.;       % missile mass, kg
T = 2000.;       % propulsive thrust, N
g = 9.81;       % gravity, m/s^2

% numerical solution parameters
global usign;
npts = 5000;     % number of time steps between t0 and tf (inclusive)
reltol = 10.;     % relative tolerance for bvp4c solution
abstol = 1.;    % absolute tolerance for bvp4c solution
usign = -1;      % 1 = num+,den-; -1 = num-,den+
mu = [1. 1. 1.];   % initial guesses for the value of mu1, mu2, mu3

%% ----- END USER INPUTS --------------------------------------------------

%% INITIALIZATION
tau0 = 0.;      % parameterized time variable
tauf = 1.;
tau_init = linspace(tau0, tauf, npts); 

param_guess = [ (X - x0) / xdot0, mu(1), mu(2) ]; % [ tf, mu1, mu2, mu3 ]

var_init = [ x0;
             z0;
             xdot0;
             zdot0;
             0.5;
             0.5;
             0.5;
             0.5 ];

solinit = bvpinit(tau_init,var_init,param_guess);
options = bvpset('Stats','on','RelTol',reltol,'AbsTol',abstol);

%% SOLUTION

sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
tau = sol.x;
y = sol.y;
tf = sol.parameters(1);
mu1 = sol.parameters(2);
mu2 = sol.parameters(3);
etaf = atan2( y(4,end), y(3,end) );
disp(['tf = ',num2str(tf),' s'])
disp(['mu1 = ',num2str(mu1),'; mu2 = ',num2str(mu2),'; mu3 = '])
disp(['The terminal impact angle is ',num2str(etaf*180/pi),' deg'])
disp(['The velocity at impact is ',num2str(sqrt(y(4,end)^2+y(3,end)^2)),' m/s'])

t = tau*tf;
%u = atan( y(8,:) ./ y(7,:) );
u = atan2( usign*y(8,:), -usign*y(7,:) );

%% PLOTTING
fig = figure;

subplot(311)
plot(y(1,:),-y(2,:))
xlabel('x (m)','FontSize',16')
ylabel('Altitude (m)','FontSize',16)
grid on;

subplot(312)
plot(t,u*180/pi)
xlabel('time (s)','FontSize',16)
ylabel('$\theta$ (deg)','FontSize',16,'Interpreter','latex')
grid on;

subplot(313)
plot(t,atan2( y(4,:), y(3,:) ) * 180/pi)
xlabel('time (s)','FontSize',16)
ylabel('Flight direction (deg)','FontSize',16)
grid on;

%% TP-BVP DEFINITION

% y(1:4) = x(1:4), y(5:8) = lambda(1:4), params = [tf mu1 mu2]

% define the ODE
function dydt = BVP_ode(tau,y,params)
    global m T g usign;
    
    tf = params(1);
    %u = atan( y(8) / y(7) );
    u = atan2( usign*y(8), -usign*y(7) );
    
    dydt = [ y(3)*tf;
             y(4)*tf;
             (T/m*cos(u))*tf;
             (g - T/m*sin(u))*tf;
             0;
             0;
             (-y(5))*tf;
             (-y(6))*tf ];
end

% define the BC's
function res = BVP_bc(ya,yb,params)
    global x0 z0 xdot0 zdot0 X Z;
    
    tf_constraint = 1 +...
		yb(5) * yb(3) +...
		yb(6) * yb(4);
    
    res = [ ya(1) - x0;         % state initial conditions
            ya(2) - z0;
            ya(3) - xdot0;
            ya(4) - zdot0;
            yb(1) - X;          % target position
            yb(2) - Z;
            yb(5) - params(2);  % co-state final conditions
            yb(6) - params(3);
            yb(7) - 0;
            yb(8) - 0;
            tf_constraint - 0 ];% from 4th Euler-Lagrange equation    
end