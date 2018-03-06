function Cm = dircosmat(phi,theta,psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function Name:            dircostmat
%   Function Author:          Steve Clark
%   Class:                    AAE 568
%   Term:                     Sp2018
%
%   SCRIPT DESCRIPTION:
%       Quaternions are computed from the heading components and 
%       compiled into a direction cosine matrix.
%   
%   NOTES:
%       Input Euler angles should be in radians.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute quaternion parameters
    A =  sin(psi/2)*sin(theta/2)*cos(phi/2) - ...
         cos(psi/2)*cos(theta/2)*sin(phi/2);
    B = -cos(psi/2)*sin(theta/2)*cos(phi/2) - ...
         sin(psi/2)*cos(theta/2)*sin(phi/2);
    C = -sin(psi/2)*cos(theta/2)*cos(phi/2) + ...
         cos(psi/2)*sin(theta/2)*sin(phi/2);
    D = -cos(psi/2)*cos(theta/2)*cos(phi/2) - ...
         sin(psi/2)*sin(theta/2)*sin(phi/2);
     
    % compute direction cosine matrix
    Cm = zeros(3);
    Cm(1,1) = A^2 - B^2 - C^2 + D^2;
    Cm(1,2) = 2 * (A*B - C*D);
    Cm(1,3) = 2 * (A*C + B*D);
    Cm(2,1) = 2 * (A*B + C*D);
    Cm(2,2) = -A^2 + B^2 - C^2 + D^2;
    Cm(2,3) = 2 * (B*C - A*D);
    Cm(3,1) = 2 * (A*C - B*D);
    Cm(3,2) = 2 * (B*C + A*D);
    Cm(3,3) = -A^2 - B^2 + C^2 + D^2;