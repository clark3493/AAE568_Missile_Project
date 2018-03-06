function [ue,ve,we,ub,vb,wb] = compute_velocities1(V,phi,theta,psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function Name:            compute_velocities1
%   Function Author:          Steve Clark
%   Class:                    AAE 568
%   Term:                     Sp2018
%
%   SCRIPT DESCRIPTION:
%       This functions converts the given velocity magnitude and aircraft
%       heading into earth fixed axis velocity components. Quaternions are
%       then computed from the heading components and compiled into a
%       direction cosine matrix. The direction cosine matrix is then used
%       to compute the body axis velocity components. Singularities in the
%       Euler angle conversion matrix are avoided.
%   
%   NOTES:
%       Input Euler angles should be in radians.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % compute the axis velocity components
    ue = V * cos(theta) * cos(psi);
    ve = V * cos(phi)   * sin(psi);
    we =-V * sin(theta);
    
    % compute the direction cosine matrix
    Cm = dircosmat(phi,theta,psi);
    
    % compute body axis velocity components
    earth_vel = [ue; ve; we];
    body_vel = Cm * earth_vel;
    
    ub = body_vel(1);
    vb = body_vel(2);
    wb = body_vel(3);