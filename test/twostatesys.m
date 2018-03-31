function [X,C,I,out] = twostatesys(inp,par)

% *************************************************************************
% Example Problem From 
% Implementation of Dynamic Programming for n-Dimensional Optimal Control
% Problems with Final State Constraints
% Philipp Elbert, Soren Ebbesen, Lino Guzzella 
% IEEE Transactions on Control Systems Technology 
% DOI: 10.1109/TCST.2012.219035
% *************************************************************************

% State Update
% Euler discretization
%X{1} = inp.X{1} + (-0.5*inp.X{1} + inp.U{1}.*inp.U{2})*inp.Ts;
%X{2} = inp.X{2} + (-0.5*inp.X{2} + inp.U{1}.*(1-inp.U{2}))*inp.Ts;

% Exact Discretization:
et = exp(-0.5*inp.Ts);
X{1} = et*inp.X{1} - 2*(et-1)*inp.U{1}.*inp.U{2};
X{2} = et*inp.X{2} - 2*(et-1)*inp.U{1}.*(1-inp.U{2});

% Cost
C{1} = (inp.U{1} + 0.1*abs(inp.U{2}-0.5))*inp.Ts;

% Feasibility
I = zeros(size(X{1}));
% constraints taken care of by dpm.m

% Output
out.u1 = inp.U{1};
out.u2 = inp.U{2};
