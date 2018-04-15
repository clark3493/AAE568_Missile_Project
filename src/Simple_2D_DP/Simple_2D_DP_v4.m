%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script Name:            Simple_2D_DP_v2.m
%   Script Author:          Steve Clark
%
%   SCRIPT DESCRIPTION
%       This script solves a simplified dynamic programming problem for the
%		missile project.
%       - Only 2D dynamics are considered
%       - Assumes that the missile is well controlled by arbitrary control 
%         surface inputs such that the missile's pitch Euler angle can be 
%         considered a control and control surface deflections are not 
%         considered.
%       - No inequality constraints (i.e. no limits on vehicle or control
%         surface angle of attack
%       - Constant thrust
%       - Aerodynamic forces are not considered
%
%       In contrast to the first version, this version of the DP
%       implementation has the time between each time step as a second
%       control input.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clearvars;

%% ----- USER INPUTS ------------------------------------------------------

% initial conditions
global x0 z0;
x0 = 0.;        % x coord, absolute coordinates, meters
z0 = -2001.;    % z coord, absolute coordinates, meters  (z = -altitude)
xdot0 = 100.;   % x velocity, absolute coordinates, m/s
zdot0 = 0.;     % z velocity, absolute coordinates, m/s

% final conditions
global x_target z_target;
x_target = 2000;       % target x absolute coordinates, meters
z_target = 1.;         % target z absolute coordinates, meters
eta = 45.;             % desired impact angle relative to horizontal, deg

% cost function weights
A = 1.;         % impact angle
B = 1.;         % time to target

% constant parameters
global m T g;
m = 100.;       % kg
T = 1000.;      % N
g = 9.81;       % m/s^2

% discretization parameters
Nx = 21;        % number of discretized points for each state variable
Nz = 21;
Nxdot = 21;
Nzdot = 21;

Nu_theta = 19;    % number of discrete missile attitude control input options
Nu_dt = 11;       % number of discrete time step size control input options

Ntau = 101;     % number of normalized time steps from initial to final state

rmin_v = 0.75;
rmax_v = 1.50;

rmin_tf = 0.75;
rmax_tf = 2.00;

% numerical soluton parameters
md_tol = .001;      % tolerance for ration between tf solution for x and z
                    % minimum distance problems before throwing error
                    
MyInf = 10000;      % finite high cost value to allow interpolation

%% ----- END USER INPUTS --------------------------------------------------

%% INITIALIZATION

% calculate solution to minimum distance problem
theta_md = fzero(@theta_md_fun,0);
V0 = sqrt(xdot0^2 + zdot0^2);
eta_md = atand((z_target-z0)/(x_target-x0));
tf_mdx = max(roots([0.5*T/m*cosd(theta_md), V0*cosd(eta_md), x0-x_target]));
tf_mdz = max(roots([0.5*(g+T/m*sind(theta_md)), V0*sind(eta_md), z0-z_target]));

if abs(tf_mdx-tf_mdz)-1 > md_tol
    error('Ratio of tf solutions for x and z min distance problems is out of tolerance')
else
    tf_md = mean([tf_mdx tf_mdz]);
end

xf_md = x0 + tf_md*V0*cosd(eta_md) + 0.5*tf_md^2*T/m*cosd(theta_md);
zf_md = z0 + tf_md*V0*sind(eta_md) + 0.5*tf_md^2*(g+T/m*sind(theta_md));
xdotf_md = V0*cosd(eta_md) + tf_md*T/m*cosd(theta_md);
zdotf_md = V0*sind(eta_md) + tf_md*(g+T/m*sind(theta_md));

if abs(xf_md-x_target) - 1 > md_tol
    error('Ratio of xf/xtarget for min distance problem is out of tolerance')
elseif abs(zf_md-z_target) - 1 > md_tol
    error('Ratio of zf/ztarget for min distance problem is out of tolerance')
elseif abs(atand(zdotf_md/xdotf_md)-eta_md) - 1 > md_tol
    error('Ratio of etaf_md/eta_md for min distance problem is out of tolerance')
end

% determine min and max values for xdot,zdot,tf
xdot_min = rmin_v * xdot0;
xdot_max = rmax_v * xdotf_md;

zdot_min = rmin_v * zdot0;
zdot_max = rmax_v * zdotf_md;

tf_min = rmin_tf * tf_md;
tf_max = rmax_tf * tf_md;

dt_max = tf_max / (Ntau-1);
dt_min = tf_min / (Ntau-1);

% discretize the state, control and time variables;
dx = (x_target - x0) / (Nx-1);
dz = (z_target - z0) / (Nz-1);
dxdot = (xdot_max - xdot_min) / (Nxdot-1);
dzdot = (zdot_max - zdot_min) / (Nzdot-1);
du_theta = 180 / (Nu_theta-1);
du_dt = (dt_max - dt_min) / (Nu_dt-1);
dtau = 1 / (Ntau-1);

x_vec = x0:dx:x_target;
z_vec = z0:dz:z_target;
xdot_vec = xdot_min:dxdot:xdot_max;
zdot_vec = zdot_min:dzdot:zdot_max;
u_theta_vec = -90:du_theta:90;
u_dt_vec = dt_min:du_dt:dt_max;
tau_vec = 0:dtau:1;

[X,Z,XDOT,ZDOT,U_THETA,U_DT] = ndgrid(x_vec,z_vec,xdot_vec,zdot_vec,...
        u_theta_vec,u_dt_vec);
[x_state,z_state,xdot_state,zdot_state] = ndgrid(x_vec,z_vec,xdot_vec,zdot_vec);
valid_final_state_inds = x_state == x_target & z_state == z_target;

%% SOLUTION

% initialize cost and control matrices
V = MyInf*ones(Nx,Nz,Nxdot,Nzdot,Ntau);
u_theta = NaN*ones(Nx,Nz,Nxdot,Nzdot,Ntau-1);
u_dt = NaN*ones(Nx,Nz,Nxdot,Nzdot,Ntau-1);

% set final costs for valid final states (all other final state costs=infinity)
impact_angle = MyInf*ones(Nx,Nz,Nxdot,Nzdot);
impact_angle(valid_final_state_inds) = atand(zdot_state(valid_final_state_inds) ./...
        xdot_state(valid_final_state_inds));
V(:,:,:,:,end) = A .* ( impact_angle - eta ) .^ 2;
clear impact_angle valid_final_state_inds;

% calculate next state for all combinations of state and control
[x_next,z_next,xdot_next,zdot_next] = f(X,Z,XDOT,ZDOT,U_THETA,U_DT);

% find optimal control and cost for every state at every time step
tstart = tic;
for k = Ntau-1:-1:1
    tstart_local = tic;
    fprintf('Calculating time step %d/%d\n',Ntau-k,Ntau-1);     
    Vfuture = MyInf*ones(Nx,Nz,Nxdot,Nzdot,Nu_theta,Nu_dt);
    for i = 1:Nu_theta
        for j = 1:Nu_dt
            Vfuture(:,:,:,:,i,j) = interpn(x_state,z_state,xdot_state,...
                zdot_state,V(:,:,:,:,k+1),x_next(:,:,:,:,i,j),z_next(:,:,:,:,i,j),...
                xdot_next(:,:,:,:,i,j),zdot_next(:,:,:,:,i,j),'linear',MyInf);
        end
    end
    cost = B*U_DT + Vfuture;
    
    [Vmin_dt,Imin_dt_theta] = min(cost,[],6,'omitnan');
    [V(:,:,:,:,k),Imin_theta] = min(Vmin_dt,[],5,'omitnan');
    
    Imin_dt = NaN*ones(Nx,Nz,Nxdot,Nzdot);
    for i = 1:Nx
        for j = 1:Nz
            for l = 1:Nxdot
                for n = 1:Nzdot
                    Imin_dt(i,j,l,n) = Imin_dt_theta( Imin_dt_theta(i,j,l,n,:) == ...
                        Imin_theta(i,j,l,n) );
                end
            end
        end
    end
    
    u_theta(:,:,:,:,k) = u_theta_vec(Imin_theta);
    u_dt(:,:,:,:,k) = u_dt_vec(Imin_dt);
    fprintf('  Calculation completed in %f seconds (%f total minutes elapsed)',...
                toc(tstart_local),toc(tstart)/60);
end




%% FUNCTION DEFINITIONS

% dynamics
function [x_,z_,xdot_,zdot_] = f(x,z,xdot,zdot,u_theta,u_dt)
    global m T g;
    
    x_ = x + u_dt .* xdot;
    z_ = z + u_dt .* zdot;
    xdot_ = xdot + u_dt .* (  T/m .* cosd(u_theta) );
    zdot_ = zdot + u_dt .* ( -T/m .* sind(u_theta) + g );
end

% theta for minimum distance solution
function out = theta_md_fun(theta)
    global m T g x0 z0 x_target z_target;
    out = (m*g + T*sind(theta))/(T*cosd(theta)) - (z_target - z0)/(x_target - x0);
end






