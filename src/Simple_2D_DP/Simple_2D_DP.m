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
clearvars; close all;

%% ----- USER INPUTS ------------------------------------------------------

% initial conditions
global x0 z0 xdot0 zdot0;
x0 = 0.;        % x coord, absolute coordinates, meters
z0 = -2000.;    % z coord, absolute coordinates, meters  (z = -altitude)
xdot0 = 100.;   % x velocity, absolute coordinates, m/s
zdot0 = 0.;     % z velocity, absolute coordinates, m/s

% final conditions
global X Z eta;
X = 2000;       % target x absolute coordinates, meters
Z = 0.;         % target z absolute coordinates, meters
eta = 45.;      % desired impact angle relative to horizontal, deg

% cost function weights
global A B;
A = 1.;         % impact angle
B = 1.;         % time to target

% constant parameters
global m T g;
m = 100.;       % kg
T = 1000.;       % N
g = 9.81;       % m/s^2

% numerical solution parameters


%% ----- END USER INPUTS --------------------------------------------------

%% INITIALIZATION


%% SOLUTION


%% PLOTTING


%% FUNCTION DEFINITIONS


