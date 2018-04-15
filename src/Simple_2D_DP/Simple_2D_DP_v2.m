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
%       This version fixes the calculations of the minimum distance problem
%       and modifies some parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;

%% ----- USER INPUTS ------------------------------------------------------

% user switches for portions of the solution
step1 = 1;      % full state + control space evaluation
if step1 == 1
    clearvars;
    step1 = 1;
end
step2 = 1;      % solution for specific initial conditions
step3 = 1;      % post processing

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
A = 0.;         % impact angle
B = 1.;         % time to target
C = 1000.;      % final position

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

Ntau = 101;     % number of normalized time steps from initial to final state

Ntf = 11;

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
if step1 == 1
    % calculate solution to minimum distance problem
    theta_md = fzero(@theta_md_fun,0);
    V0 = sqrt(xdot0^2 + zdot0^2);
    eta_md = atand((z_target-z0)/(x_target-x0));
    tf_mdx = max(roots([0.5*T/m*cosd(theta_md), V0*cosd(eta_md), x0-x_target]));
    tf_mdz = max(roots([0.5*(g+T/m*sind(theta_md)), V0*sind(eta_md), z0-z_target]));

    if abs(tf_mdx-tf_mdz) > md_tol
        error('Ratio of tf solutions for x and z min distance problems is out of tolerance')
    else
        tf_md = mean([tf_mdx tf_mdz]);
    end

    xf_md = x0 + tf_md*V0*cosd(eta_md) + 0.5*tf_md^2*T/m*cosd(theta_md);
    zf_md = z0 + tf_md*V0*sind(eta_md) + 0.5*tf_md^2*(g+T/m*sind(theta_md));
    xdotf_md = V0*cosd(eta_md) + tf_md*T/m*cosd(theta_md);
    zdotf_md = V0*sind(eta_md) + tf_md*(g+T/m*sind(theta_md));

    if abs(xf_md-x_target) > md_tol
        error('Ratio of xf/xtarget for min distance problem is out of tolerance')
    elseif abs(zf_md-z_target) > md_tol
        error('Ratio of zf/ztarget for min distance problem is out of tolerance')
    elseif abs(atand(zdotf_md/xdotf_md)-eta_md) > md_tol
        error('Ratio of etaf_md/eta_md for min distance problem is out of tolerance')
    end

    % determine min and max values for xdot,zdot,tf
    x_min = x0 -100;
    x_max = x_target + 100;
    
    z_min = z0 - 100;
    z_max = z_target + 100;
    
    xdot_min = rmin_v * xdot0;
    xdot_max = rmax_v * xdotf_md;

    zdot_min = rmin_v * zdot0 - 100;
    zdot_max = rmax_v * zdotf_md;

    tf_min = rmin_tf * tf_md;
    tf_max = rmax_tf * tf_md;

    % discretize the state, control and time variables;
    dx = (x_max - x_min) / (Nx-1);
    dz = (z_max - z_min) / (Nz-1);
    dxdot = (xdot_max - xdot_min) / (Nxdot-1);
    dzdot = (zdot_max - zdot_min) / (Nzdot-1);
    du_theta = 180 / (Nu_theta-1);
    dtau = 1 / (Ntau-1);
    dtf = (tf_max - tf_min) / (Ntf-1);

    x_vec = x_min:dx:x_max;
    z_vec = z_min:dz:z_max;
    xdot_vec = xdot_min:dxdot:xdot_max;
    zdot_vec = zdot_min:dzdot:zdot_max;
    u_theta_vec = -90:du_theta:90;
    tau_vec = 0:dtau:1;
    tf_vec = tf_min:dtf:tf_max;

    [X,Z,XDOT,ZDOT,U_THETA] = ndgrid(x_vec,z_vec,xdot_vec,zdot_vec,u_theta_vec);
    [x_state,z_state,xdot_state,zdot_state] = ndgrid(x_vec,z_vec,xdot_vec,zdot_vec);

    %% SOLUTION

    % initialize cost and control matrices
    V = MyInf*ones(Nx,Nz,Nxdot,Nzdot,Ntau,Ntf);
    uind = NaN*ones(Nx,Nz,Nxdot,Nzdot,Ntau-1,Ntf);
    u_theta = NaN*ones(Nx,Nz,Nxdot,Nzdot,Ntau-1,Ntf);

    % set final costs for valid final states (all other final state costs=infinity)
    impact_angle = atand(zdot_state ./ xdot_state);
    for i = 1:Ntf
        V(:,:,:,:,end,i) = A .* ( impact_angle - eta ) .^ 2 + ...
            C * sqrt( (x_state-x_target).^2 + (z_state-z_target).^2 );
    end
    clear impact_angle valid_final_state_inds;

    tstart = tic;
    for itf = 1:Ntf
        tf = tf_vec(itf);
        dt = tf/(Ntau-1);
        disp(['Performing solution for tf = ',num2str(tf),' sec'])
        tstart_local = tic;
        
        % calculate immediate cost
        L = B * dtau * tf;      % same for every combination of state and control
        
        % calculate next state for all combinations of state and control
        [x_next,z_next,xdot_next,zdot_next] = f(X,Z,XDOT,ZDOT,U_THETA,dt);
        
        % find optimal control and cost for every state at every time step
        for k = Ntau-1:-1:1
            Vfuture = MyInf*ones(Nx,Nz,Nxdot,Nzdot,Nu_theta);
            for i = 1:Nu_theta
                Vfuture(:,:,:,:,i) = interpn(x_state,z_state,xdot_state,zdot_state,...
                    V(:,:,:,:,k+1,itf),x_next(:,:,:,:,i),z_next(:,:,:,:,i),...
                    xdot_next(:,:,:,:,i),zdot_next(:,:,:,:,i),'linear',MyInf);
            end
            cost = L + Vfuture;
            
            [V(:,:,:,:,k,itf),uind(:,:,:,:,k,itf)] = min(cost,[],5);
            u_theta(:,:,:,:,k,itf) = u_theta_vec(uind(:,:,:,:,k,itf));
        end
        
        disp(['Finished local solution in ',num2str(round(toc(tstart_local),1)),' sec'])
        disp([num2str(round(toc(tstart)/60,2)),' minutes elapsed - ',num2str(round(itf/Ntf*100)),'% complete with full solution'])
        disp(' ')
    end

    tend = toc(tstart);
    disp(['TOTAL STATE + CONTROL SPACE SOLUTION CALCULATED IN ',num2str(round(tend/60,1)),' MINUTES'])
end

%% SOLUTION FOR SPECIFIED INITIAL CONDITIONS
if step2 == 1
    % find the value of tf for which the minimum cost occurs
    Vtotal = NaN*ones(Ntf,1);
    for i = 1:Ntf
        Vtotal(i) = interpn(x_state,z_state,xdot_state,zdot_state,V(:,:,:,:,1,i),...
            x0,z0,xdot0,zdot0,'linear',MyInf);
    end
    
    [Vmin,Imin] = min(Vtotal);
    
    u_actual = NaN*ones(1,Ntau-1);
    alpha = NaN*ones(1,Ntau);
    x = NaN*ones(4,Ntau);
    x(:,1) = [x0;z0;xdot0;zdot0];
    
    for k = 1:Ntau-1
        u_actual(k) = interpn(x_state,z_state,xdot_state,zdot_state,...
            u_theta(:,:,:,:,k,Imin),x(1,k),x(2,k),x(3,k),x(4,k));
        alpha(k) = u_actual(k) + atand( x(3,k) / x(4,k) );
        [x(1,k+1),x(2,k+1),x(3,k+1),x(4,k+1)] = f(x(1,k),x(2,k),x(3,k),x(4,k),u_actual(k),dtau*tf_vec(Imin));
    end
    
    V_plot = NaN*ones(1,Ntau);
    for k = 1:Ntau
        V_plot(k) = interpn(x_state,z_state,xdot_state,zdot_state,V(:,:,:,:,k,Imin),...
            x(1,k),x(2,k),x(3,k),x(4,k));
    end
end

if step3 == 1
    fig = figure;
    
    subplot(321)
    plot(x(1,:),-x(2,:)); hold on;
    xlabel('x')
    ylabel('altitude')
    
    subplot(322)
    plot(tf_vec(Imin)*tau_vec,x(3,:)); hold on;
    plot(tf_vec(Imin)*tau_vec,x(4,:));
    legend('xdot','zdot')
    xlabel('Time (s)')
    
    subplot(323)
    plot(tf_vec(Imin)*tau_vec(1:end-1),u_actual)
    xlabel('Time (s)')
    ylabel('u')
    
    subplot(324)
    plot(tf_vec(Imin)*tau_vec,V_plot)
    xlabel('Time (s)')
    ylabel('cost')
    
    subplot(325)
    plot(tf_vec(Imin)*tau_vec,alpha)
    xlabel('Time (s)')
    ylabel('$\alpha','Interpreter','latex')
end

%% FUNCTION DEFINITIONS

% dynamics
function [x_,z_,xdot_,zdot_] = f(x,z,xdot,zdot,u_theta,dt)
    global m T g;
    
    x_ = x + dt .* xdot;
    z_ = z + dt .* zdot;
    xdot_ = xdot + dt .* (  T/m .* cosd(u_theta) );
    zdot_ = zdot + dt .* ( -T/m .* sind(u_theta) + g );
end

% theta for minimum distance solution
function out = theta_md_fun(theta)
    global m T g x0 z0 x_target z_target;
    out = (m*g + T*sind(theta))/(T*cosd(theta)) - (z_target - z0)/(x_target - x0);
end






