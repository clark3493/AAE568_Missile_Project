%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script Name:            Simple_2D_DP.m
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;

%% ----- USER INPUTS ------------------------------------------------------

% user switches for portions of the solution
step1 = 0;      % full state + control space evaluation
if step1 == 1
    clearvars;
    step1 = 1;
end
step2 = 1;      % solution for specific initial conditions
step3 = 1;      % post processing

% initial conditions
global x0 z0 xdot0 zdot0;
x0 = 0.;        % x coord, absolute coordinates, meters
z0 = -2000.;    % z coord, absolute coordinates, meters  (z = -altitude)
xdot0 = 100.;   % x velocity, absolute coordinates, m/s
zdot0 = 0.;     % z velocity, absolute coordinates, m/s

% final conditions
global x_target z_target eta;
x_target = 2000;       % target x absolute coordinates, meters
z_target = 0.;         % target z absolute coordinates, meters
eta = 45.;             % desired impact angle relative to horizontal, deg

% cost function weights
global A B;
A = 1.;         % impact angle
B = 1.;         % time to target

% constant parameters
global m T g;
m = 100.;       % kg
T = 1000.;      % N
g = 9.81;       % m/s^2

% discretization and numerical solution parameters
Nx = 21;       % number of discretized points for each variable
Nz = 21;
Nxdot = 21;
Nzdot = 21;

Nu = 19;
Ntau = 101;
Ntf = 11;

rmin = 1.0;    % ratio to use when determining the minimum value of x3,x4,tf
rmax = 1.5;     % ratio to use when determining the maximum value of x3,x4,tf

MyInf = 10000;  % high cost value to allow interpolation

global x_vec z_vec xdot_max xdot_min zdot_max zdot_min;

%% ----- END USER INPUTS --------------------------------------------------

%% FULL STATE + CONTROL SPACE COST AND OPTIMAL CONTROL EVALUATION
if step1 == 1
    %% INITIALIZATION

    % calculate solution to minimum distance problem
    eta_md = atand((z_target-z0)/(x_target-x0));
    tf_md = max(roots([.5*T/m*cosd(eta_md), xdot0, x0-x_target]));
    xdot_mdf = xdot0 + tf_md*T/m*cosd(eta_md);
    zdot_mdf = zdot0 + tf_md*(g + T/m*sind(eta_md));

    % determine min and max value for x3,x4,tf
    xdot_min = rmin * xdot0;
    xdot_max = rmax * xdot_mdf;

    zdot_min = rmin * zdot0;
    zdot_max = rmax * zdot_mdf;

    tf_min = rmin * tf_md;
    tf_max = rmax * tf_md;

    % discretize the state, control and time variables
    dx = (x_target - x0) / (Nx-1);
    dz = (z_target - z0) / (Nz-1);
    dxdot = (xdot_max - xdot_min) / (Nxdot-1);
    dzdot = (zdot_max - zdot_min) / (Nzdot-1);
    du = 180 / (Nu-1);
    dtau = 1 / (Ntau-1);
    dtf = (tf_max - tf_min) / (Ntf-1);

    x_vec = x0:dx:x_target;
    z_vec = z0:dz:z_target;
    xdot_vec = xdot_min:dxdot:xdot_max;
    zdot_vec = zdot_min:dzdot:zdot_max;
    u_vec = -90:du:90;
    tau = 0:dtau:1;
    tf_vec = tf_min:dtf:tf_max;

    [X,Z,XDOT,ZDOT,U] = ndgrid(x_vec,z_vec,xdot_vec,zdot_vec,u_vec);
    [x_state_grid,z_state_grid,xdot_state_grid,zdot_state_grid] = ...
        ndgrid(x_vec,z_vec,xdot_vec,zdot_vec);
    f_inds = x_state_grid == x_target & z_state_grid == z_target;

    %% SOLUTION

    % initialize cost and control matrices
    V = MyInf*ones(Nx,Nz,Nxdot,Nzdot,Ntau,Ntf);
    uind = NaN*ones(Nx,Nz,Nxdot,Nzdot,Ntau-1,Ntf);
    u = NaN*ones(Nx,Nz,Nxdot,Nzdot,Ntau-1,Ntf);

    % set final costs for valid final states (all other final state costs=infinity)
    temp = MyInf*ones(Nx,Nz,Nxdot,Nzdot);
    temp(f_inds) = atand(zdot_state_grid(f_inds)./xdot_state_grid(f_inds));
    for i = 1:Ntf
        V(:,:,:,:,end,i) = ...
            A .* ( temp - eta ) .^ 2;
    end
    clear temp f_inds;

    tstart = tic;
    % iterate over each value of tf
    for itf = 1:length(tf_vec)
        tf = tf_vec(itf);
        dt = tf/(Ntau-1);
        disp(['Performing solution for tf = ',num2str(tf),' sec'])
        tstart_local = tic;

        % calculate immediate cost for each time step
        L = B * dtau * tf;      % same for every combination of state and control

        % calculate next state for all combinations of state and control
        [x_next,z_next,xdot_next,zdot_next] = f(X,Z,XDOT,ZDOT,U,dt);
        
        % find optimal control and cost for every state at every time step
        for k = Ntau-1:-1:1
            Vfuture = MyInf*ones(Nx,Nz,Nxdot,Nzdot,Nu);
            for i = 1:Nu
                Vfuture(:,:,:,:,i) = interpn(x_state_grid,z_state_grid,xdot_state_grid,...
                    zdot_state_grid,V(:,:,:,:,k+1,itf),x_next(:,:,:,:,i),z_next(:,:,:,:,i),...
                    xdot_next(:,:,:,:,i),zdot_next(:,:,:,:,i),'linear',10000);
            end
            cost = L + Vfuture;
            
            [V(:,:,:,:,k,itf),uind(:,:,:,:,k,itf)] = min(cost,[],5);
            u(:,:,:,:,k,itf) = u_vec(uind(:,:,:,:,k,itf));
        end

        disp(['Finished local solution in ',num2str(round(toc(tstart_local),1)),' sec'])
        disp([num2str(round(toc(tstart)/60,2)),' minutes elapsed - ',num2str(round(itf/Ntf*100)),'% complete with full solution'])
        disp(' ')
    end

    tend = toc(tstart);
    disp(['TOTAL STATE + CONTROL SPACE SOLUTION CALCULATED IN ',num2str(round(toc(tstart)/60,1)),' MINUTES'])
end

%% SOLUTION FOR SPECIFIED INITIAL CONDITIONS
if step2 == 1
    
    % find the value of tf for which the minimum cost occurs
    Vtotal = NaN*ones(Ntf,1);
    for i = 1:Ntf
        Vtotal(i) = interpn(x_state_grid,z_state_grid,xdot_state_grid,zdot_state_grid,...
                V(:,:,:,:,1,i),x0,z0,xdot0,zdot0);
    end
    
    [Vmin,Imin] = min(Vtotal);
    
    u_actual = NaN*ones(1,Ntau-1);
    x = NaN*ones(4,Ntau);
    x(:,1) = [x0;z0;xdot0;zdot0];
    
    for k = 1:Ntau-1
        u_actual(k) = interpn(x_state_grid,z_state_grid,xdot_state_grid,zdot_state_grid,...
                u(:,:,:,:,k,Imin),x(1,k),x(2,k),x(3,k),x(4,k));
        [x(1,k+1),x(2,k+1),x(3,k+1),x(4,k+1)] = f(x(1,k),x(2,k),x(3,k),x(4,k),u_actual(k),dtau*tf_vec(Imin));
    end
    
    V_plot = NaN*ones(1,Ntau);
    for k = 1:Ntau
        V_plot(k) = interpn(x_state_grid,z_state_grid,xdot_state_grid,zdot_state_grid,...
                V(:,:,:,:,k,Imin),x(1,k),x(2,k),x(3,k),x(4,k));
    end
end

if step3 == 1
    fig = figure;
    
    subplot(221)
    plot(x(1,:),-x(2,:)); hold on;
    xlabel('x')
    ylabel('altitude')
    
    subplot(222)
    plot(tf_vec(Imin)*tau,x(3,:)); hold on;
    plot(tf_vec(Imin)*tau,x(4,:));
    legend('xdot','zdot')
    xlabel('Time (s)')
    
    subplot(223)
    plot(tf_vec(Imin)*tau(1:end-1),u_actual)
    xlabel('Time (s)')
    ylabel('u')
    
    subplot(224)
    plot(tf_vec(Imin)*tau,V_plot)
    xlabel('Time (s)')
    ylabel('cost')
end



%% FUNCTION DEFINITIONS

% dynamics
function [x_,z_,xdot_,zdot_] = f(x,z,xdot,zdot,u,dt)
    global m T g;
    global x_vec z_vec xdot_max xdot_min zdot_max zdot_min;
    
    x_ = x + dt .* xdot;
    z_ = z + dt .* zdot;
    xdot_ = xdot + dt .* ( T/m .* cosd(u) );
    zdot_ = zdot + dt .* ( g - T/m .* sind(u) );
    
    %x_(x_ > max(x_vec)) = max(x_vec);
    %x_(x_ < min(x_vec)) = min(x_vec);
    %z_(z_ > max(z_vec)) = max(z_vec);
    %z_(z_ < min(z_vec)) = min(z_vec);
    %xdot_(xdot_>xdot_max) = xdot_max;
    %xdot_(xdot_<xdot_min) = xdot_min;
    %zdot_(zdot_>zdot_max) = zdot_max;
    %zdot_(zdot_<zdot_min) = zdot_min;
end