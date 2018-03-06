classdef ControlSurface < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Class Name:               ControlSurface
%   Class Author:             Steve Clark
%   Class:                    AAE 568
%   Term:                     Sp2018
%
%   SCRIPT DESCRIPTION:
%       Defines a control surface object to be added to a missile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        angle        % angle of the undeflected control surface wrt x-x axis
        area         % surface area
        x            % x-coordinate of center of lift, as measured from tip
        
        def          % instantaneous deflection from neutral position
        
        movable      % boolean to define mmovable or stationary
        
        CLalpha      % slope of list curve wrt AoA
        CL           % lift coefficient
        CD           % drag coefficient
       
        alpha_max    % maximum AoA
        def_rate_max % maximum time rate of deflection
    end
    methods
        function obj = ControlSurface(angle,area,x,varargin)
            obj.angle = angle;
            obj.area = area;
            obj.x = x;
            
            if nargin > 3
                obj.movable = varargin{1};
            end
            if nargin > 4
                obj.CLalpha = varargin{2};
            end
            if nargin > 5
                obj.alpha_max = varargin{3};
            end
            if nargin > 6
                obj.def_rate_max = varargin{4};
            end
        end
    end
end