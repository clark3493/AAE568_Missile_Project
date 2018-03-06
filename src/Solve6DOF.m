function x = Solve6DOF(F,x0,P,dt,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function Name:            Solve6DOF
%   Function Author:          Steve Clark
%   Class:                    AAE 568
%   Term:                     Sp2018
%
%   SCRIPT DESCRIPTION:
%       Calculates the 6 degree-of-freedom vehicle state variables at the
%       subsequent timestep given the summed forces in each of the 6 
%       directions and initial states
%   
%   INPUTS:
%       F = 6 element vector containing the summed forces and moments along
%           the 3 body axes
%       x0 = 6 element vector containing initial state for each of the 6 
%           DOF's in body axes
%       P = vector of vehicle properties [m, Ix, Iy, Iz]
%       dt = timestep to next frame
%       Solution method (optional) = define the solution method to use
%           Default = 'Euler'
%           Additional options = 
%
%   OUTPUTS:
%       x = 6 element vector containing the linear and angular velocities
%           of the vehicle in body axes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input error checking
if length(F) ~= 6 || length(x0) ~= 6 || length(P) ~= 4
    error('Input force and initial state vector must have 6 elements');

% get specified solution method if any
solution_options = {};
if nargin > 0
    sol = varargin{1};
    if ~any(strcmp(solution_options,sol))
        error('Specified solution method is not a valid option');
    end
else
    sol = 'Euler';
end

% extract inputs for readability
u = x0(1);
v = x0(2);
w = x0(3);
p = x0(4);
q = x0(5);
r = x0(6);

X = F(1);
Y = F(2);
Z = F(3);
L = F(4);
M = F(5);
N = F(6);

m = P(1);
Ix = P(2);
Iy = P(3);
Iz = P(4);

%% EULER STEP SOLUTION
if strcmp(sol,'Euler')
    unew = u + dt * (X/m - w*q + v*r);
    vnew = v + dt * (Y/m - u*r + w*p);
    wnew = w + dt * (Z/m - v*p + u*q);
    
    pnew = p + dt * (q*r*(Iy-Iz)/Ix + L/Ix);
    qnew = q + dt * (p*r*(Iz-Ix)/Iy + M/Iy);
    rnew = r + dt * (p*q*(Ix-Iy)/Iz + N/Iz);
end

x = [unew,vnew,wnew,pnew,qnew,rnew];

