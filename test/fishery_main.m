% THIS SCRIPT IS AN EXAMPLE THAT CAME WITH THE DPM SOURCE CODE

addpath('../src')

% create grid
clear grd
grd.Nx{1} = 201;
grd.Xn{1}.lo = 0;
grd.Xn{1}.hi = 1000;
grd.Nu{1} = 21;
grd.Un{1}.lo = 0;
grd.Un{1}.hi = 10;

% set initial state
grd.X0{1} = 250;

% set final state constraints
grd.XN{1}.hi = 1000;
grd.XN{1}.lo = 750;

% define problem
clear prb
prb.Ts = 1/5;
prb.N = 200*1/prb.Ts + 1;

% set options
options = dpm();
options.BoundaryMethod = 'Line'; % also possible: 'none' or 'LevelSet';
if strcmp(options.BoundaryMethod,'Line')
    options.FixedGrid = 1;
    options.Iter = 10;
    options.Tol = 1e-8;
end
[res dyn] = dpm(@fishery,[],grd,prb,options);