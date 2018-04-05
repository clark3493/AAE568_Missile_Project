clear all
close all

par = [];

prb.Ts = 0.01;
prb.N  = .78/prb.Ts;

options = dpm();
Nx1 = 51;
Nx2 = 51;
Nu  = 51;

grd.X0{1}    =  .05;
grd.Xn{1}.lo = -.25;
grd.Xn{1}.hi =  .25;
grd.XN{1}.lo = -.25;
grd.XN{1}.hi =  .25;
grd.Nx{1}    =  Nx1;

grd.X0{2}    =  .00;
grd.Xn{2}.lo = -.25;
grd.Xn{2}.hi =  .25;
grd.XN{2}.lo = -.25;
grd.XN{2}.hi =  .25;
grd.Nx{2}    =  Nx2;

grd.Un{1}.lo =  .00;
grd.Un{1}.hi = 1.10;
grd.Nu{1}    =   Nu;

tic
options.BoundaryMethod = 'none';
[out dyn] = dpm(@dxdt,par,grd,prb,options);
tb = toc;

t = 0:prb.Ts:0.78;

fig = figure;

subplot(211)
plot(t,out.X{1}); hold on;
plot(t,out.X{2});
legend('x1','x2')

subplot(212)
plot(t(1:end-1),out.U{1})
legend('u')




function [X,C,I,signals] = dxdt(inp,par)
    % state update
    x1 = inp.X{1};
    x2 = inp.X{2};

    
    f1 = -2*(x1+.25)+(x2+.5).*exp(25*x1./(x1+2))-(x1+.25).*inp.U{1};
    f2 = 0.5 - x2 - (x2+0.5).*exp(25*x1./(x1+2));
    
    X{1} = inp.Ts*f1 + x1;
    X{2} = inp.Ts*f2 + x2;
    
    % cost
    C{1} = x1.^2 + x2.^2 + 0.1*inp.U{1}.^2;
    
    % infeasibility
    I = 0;
    signals.U{1} = inp.U{1};
end