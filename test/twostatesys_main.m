% *************************************************************************
% Example Problem From 
% Implementation of Dynamic Programming for n-Dimensional Optimal Control
% Problems with Final State Constraints
% Philipp Elbert, Soren Ebbesen, Lino Guzzella 
% IEEE Transactions on Control Systems Technology 
% DOI: 10.1109/TCST.2012.219035
% *************************************************************************
addpath('../src')
%--------------------------------------------------------------------------
% CLEAR ALL:
%--------------------------------------------------------------------------
clear all
close all
clc

%--------------------------------------------------------------------------
% SIMULATION MODEL:
%--------------------------------------------------------------------------
model = 'twostatesys';
par = [];

%--------------------------------------------------------------------------
% PROBLEM & OPTIONS STRUCTURE:
%--------------------------------------------------------------------------
prb.Ts = 0.01;
prb.N  = 2/prb.Ts;
prb.N0 = 1;

options    = dpm();
Nx = 51;

%--------------------------------------------------------------------------
% PREPARE DP:
%--------------------------------------------------------------------------
grd.X0{1} = 0;
grd.Xn{1}.lo = 0;
grd.Xn{1}.hi = 1;
grd.XN{1}.lo = 0.5;
grd.XN{1}.hi = 1;
grd.Nx{1} = Nx;

grd.X0{2} = 0;
grd.Xn{2}.lo = 0;
grd.Xn{2}.hi = 1;
grd.XN{2}.lo = 0.5;
grd.XN{2}.hi = 1;
grd.Nx{2} = Nx;

grd.Un{1}.hi = 1;
grd.Un{1}.lo = 0;
grd.Nu{1} = 21;

grd.Un{2}.hi = 1;
grd.Un{2}.lo = 0;
grd.Nu{2} = 21;

options.MyInf = sum(grd.Un{1}.hi*ones(prb.N,1)*prb.Ts);


%--------------------------------------------------------------------------
% BASIC DP:
%--------------------------------------------------------------------------
tic
options.BoundaryMethod = 'none';
[outb dynb] = dpm(model,par,grd,prb,options);
tmeb = toc;

%--------------------------------------------------------------------------
% LEVEL SET DP:
%--------------------------------------------------------------------------
options.BoundaryMethod = 'LevelSet';
x1 = linspace(grd.Xn{1}.lo,grd.Xn{1}.hi,grd.Nx{1});
x2 = linspace(grd.Xn{2}.lo,grd.Xn{2}.hi,grd.Nx{2});
[X1,X2] = ndgrid(x1,x2);
options.gN{1} = max(grd.XN{1}.lo-X1,0) + max(grd.XN{2}.lo-X2,0);
tic
[out dyn] = dpm(model,par,grd,prb,options);
tme = toc;

%SAVE RESULTS
save results

%%
%--------------------------------------------------------------------------
% PLOT RESULTS:
%--------------------------------------------------------------------------
%Calculate Optimal Solution:
t = [0:prb.Ts:2];
ton = 2+2*log(0.5);
xopt = max(0,1-exp(-0.5*(t-ton)));
u1 = t>=ton;
u2 = 0.5*(t>=0);

%Plot State Trajectory
fgr = figure(1);
clf
ax(1) = subplot(211);
k=2;
plot(t,xopt,'Color',[0.8 0.8 0.8],'LineWidth',2)
hold on
plot(t,outb.X{1})
plot(t,out.X{1},'--','Color',[0 0.5 0]);
plot([.35 0.45],[0.3 0.3])
text(.5,.3,'basic DP')
plot([1.2 1.3],[0.2 0.2],'--','Color',[0 0.5 0])
text(1.35,.2,'level-set DP')
plot([1.2 1.3],[0.1 0.1],'Color',[0.8 0.8 0.8],'LineWidth',2)
text(1.35,.1,'analytic solution')
ylabel('x_1, x_2 [-]')

%Plot Control Inputs:
ax(2) = subplot(212);
plot(t,u2,'Color',[0.8 0.8 0.8],'LineWidth',2);
hold on
plot(t(1:end-1),outb.u2,'b')
plot(t(1:end-1),out.u2,'--','Color',[0 0.5 0])
ylabel('$u_2$ [-]')
xlabel('time [s]')
plot(t,u1,'Color',[0.8 0.8 0.8],'LineWidth',2);
plot(t(1:end-1),outb.u1,'b')
plot(t(1:end-1),out.u1,'--','Color',[0 0.5 0]);
ylabel('u_1, u_2 [-]')
