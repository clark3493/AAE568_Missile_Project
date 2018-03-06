classdef Missile < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Class Name:               Missile
%   Class Author:             Steve Clark
%   Class:                    AAE 568
%   Term:                     Sp2018
%
%   SCRIPT DESCRIPTION:
%       Defines a missile object for use in simulations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        m         % mass 
        d         % diameter
        Ix        % Moments of inertia
        Iy
        Iz
        l         % length
        cg        % cg, distance from tip
        
        CS        % vector of control surface objects
        
        T_max      % maximum thrust
        alpha_max  % maximum angle of attack
        
        CLalpha    % slope of lift curve wrt AoA
        CL         % lift coefficient
        CD         % drag coefficient
        
        x          % absolute position in inertial space
        y
        z
        
        u          % linear velocity in body axes
        v
        w
        
        p          % angular velocity in body axes
        q
        r
    end
    methods
        
        function obj = AddCS(obj,angle,area,x,varargin)
            cs = ControlSurface(angle,area,x,varargin{:});
            obj.CS = [obj.CS, cs];
        end
        
    end
end