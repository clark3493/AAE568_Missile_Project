clearvars;

X1 = -0.1:0.01:0.3;
X2 = -0.1:0.01:0.3;
U  =  0. :0.1:2.0;

nx1 = length(X1);
nx2 = length(X2);
nu  = length(U);

tf = 0.78;
dt = .01;

% costs
gxu = Inf*ones(nx1,nx2,nu);

for i = 1:nx1
    for j = 1:nx2
        for k = 1:nu
            gxu(i,j,k) = 1.*X1(i)^2 + 1.*X2(j)^2 + 0.1*U(k)^2;
        end
    end
end

V = Inf*ones(nx1,nx2,tf/dt+1);
V(:,:,end) = zeros(nx1,nx2);
u = NaN*ones(nx1,nx2,tf/dt);

for t = tf/dt:-1:1
    [V(:,:,t),u(:,:,t)] = min(gxu + repmat(V(:,:,t+1),1,1,nu),[],3);
end

x1 = NaN*ones(1,tf/dt+1);
x2 = NaN*ones(1,tf/dt+1);
u_opt = NaN*ones(1,tf/dt);

x1(1) = 0.05;
x2(1) = 0.;

[x1_grid, x2_grid] = meshgrid(X1,X2);

for t = 1:tf/dt
    u_opt(t) = interp2(x1_grid,x2_grid,u(:,:,t)',x1(t),x2(t));
    [x1(t+1),x2(t+1)] = f(x1(t),x2(t),u_opt(t),dt);
end

t = 0:dt:tf;

fig = figure;

subplot(211)
plot(t,x1,'r'); hold on;
plot(t,x2,'g');
legend('x1','x2');

subplot(212)
plot(t(1:end-1),u_opt,'k');


function [x1,x2] = f(x1_,x2_, u, dt)
    x1 = dt * ( -2*(x1_+.25) + (x2_ + .5)*exp(25*x1_/(x1_+2)) - (x1_+.25)*u);
    x2 = dt * ( 0.5 - x2_ - (x2_+0.5)*exp(25*x1_/(x1_+2)));
end