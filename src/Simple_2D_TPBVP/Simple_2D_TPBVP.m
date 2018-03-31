%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script Name:            Simple_2D_TPBVP.m
%   Script Author:          Steve Clark
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
eta = 45.;      % desired impact angle relative to horizontal, deg

% cost function weights
global A B;
A = 1.;         % impact angle
B = 0.;         % time to target

% constant parameters
global m T g;
m = 100.;       % missile mass, kg
T = 1000.;       % propulsive thrust, N
g = 9.81;       % gravity, m/s^2

% numerical solution parameters
global usign;
npts = 5000;     % number of time steps between t0 and tf (inclusive)
reltol = .1;     % relative tolerance for bvp4c solution
usign = -1;      % 1 = num+,den-; -1 = num-,den+
mu = [.074 .187];   % initial guesses for the value of mu1 and mu2

%% ----- END USER INPUTS --------------------------------------------------

%% INITIALIZATION
tau0 = 0.;      % parameterized time variable
tauf = 1.;
tau_init = linspace(tau0, tauf, npts); 

param_guess = [ (X - x0) / xdot0, mu(1), mu(2) ]; % [ tf, mu1, mu2 ]

var_init = [ x0;
             z0;
             xdot0;
             zdot0;
             0.5;
             0.5;
             0.5;
             0.5 ];

solinit = bvpinit(tau_init,var_init,param_guess);
options = bvpset('Stats','on','RelTol',reltol);

%% SOLUTION

sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
tau = sol.x;
y = sol.y;
tf = sol.parameters(1);
disp(['tf = ',num2str(tf),' s'])
disp(['mu1 = ',num2str(sol.parameters(2)),'; mu2 = ',num2str(sol.parameters(3))])

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
             (m/T*cos(u))*tf;
             (g - m/T*sin(u))*tf;
             0;
             0;
             (-y(5))*tf;
             (-y(6))*tf ];
end

% define the BC's
function res = BVP_bc(ya,yb,params)
    global x0 z0 xdot0 zdot0 X Z A B eta m T g usign;
    
    %uf = atan( yb(8) / yb(7) );
    uf = atan2( usign*yb(8), -usign*yb(7) );
    
    % intermediate calculations
    lambda3_tf = 2*A * yb(4) / ( yb(4)^2 + yb(3)^2 ) * ...
        ( eta*pi/180 - atan( yb(4) / yb(3) ) );
        %( eta*pi/180 - atan2( yb(4), yb(3) ) );
    lambda4_tf = 2*A / ( yb(4)^2 / yb(3) + yb(3) ) * ...
        ( eta*pi/180 + atan( yb(4) / yb(3) ) );
        %( eta*pi/180 + atan2( yb(4), yb(3) ) );
    tf_constraint = A * (...
		(1/yb(3)*(g-m/T*sin(uf))-yb(4)/yb(3)^2*m/T*cos(uf)) *...
			(2*atan(yb(4)/yb(3))+eta*pi/180)/((yb(4)/yb(3))^2+1) ) +...
		B +...
		yb(5) * yb(3) +...
		yb(6) * yb(4) +...
		yb(7) * m/T*cos(uf) +...
		yb(8) * ( g - m/T*sin(uf) );
    
    res = [ ya(1) - x0;         % state initial conditions
            ya(2) - z0;
            ya(3) - xdot0;
            ya(4) - zdot0;
            yb(1) - X;          % target position
            yb(2) - Z;
            yb(5) - params(2);  % co-state final conditions
            yb(6) - params(3);
            yb(7) - lambda3_tf;
            yb(8) - lambda4_tf;
            tf_constraint - 0 ];% from 4th Euler-Lagrange equation    
end