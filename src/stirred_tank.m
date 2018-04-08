clearvars; close all;

%% discretize the state, control and time
x1_min = -0.2;
x1_max =  0.2;
dx1 = .01;

x2_min = -0.2;
x2_max =  0.2;
dx2 = .01;

u_min = 0;
u_max = 2;
du = .01;

t0 = 0;
tf = .78;
dt = .01;

x1_vec = x1_min:dx1:x1_max;
x2_vec = x2_min:dx2:x2_max;
u_vec = u_min:du:u_max;
t = t0:dt:tf;

Nx1 = length(x1_vec);
Nx2 = length(x2_vec);
Nu = length(u_vec);
Nt = length(t);

%% calculate the immediate cost for every combination of state and control
[X1,X2,U] = ndgrid(x1_vec,x2_vec,u_vec);
L = X1.^2 + X2.^2 + .1*U.^2;

%% initialize the cost storage and control matrices
V = Inf*ones(Nx1,Nx2,Nt);
V(:,:,end) = 0;     % no terminal state cost
u = NaN*ones(Nx1,Nx2,Nt-1);
uind = NaN*ones(Nx1,Nx2,Nt-1);

% calculate the next state for every combination of current state and control
[x1_next,x2_next] = f(X1,X2,U,dt);

% find optimal control and cost for every state at every time step
[x1_state_grid, x2_state_grid] = ndgrid(x1_vec,x2_vec);
for k = Nt-1:-1:1
    Vfuture = Inf*ones(Nx1,Nx2,Nu);
    for i = 1:Nu
        Vfuture(:,:,i) = interpn(x1_state_grid,x2_state_grid,V(:,:,k+1),...
            x1_next(:,:,i),x2_next(:,:,i));
    end
    cost = L+Vfuture;
    [V(:,:,k),uind(:,:,k)] = min(cost,[],3);
    u(:,:,k) = u_vec(uind(:,:,k));
end

% calculate solution for given initial conditions
u_actual = NaN*ones(1,Nt-1);
x = NaN*ones(2,Nt);
x(:,1) = [.05; 0];

for k = 1:Nt-1
    u_actual(k) = interpn(x1_state_grid,x2_state_grid,u(:,:,k),x(1,k),x(2,k));
    [x(1,k+1), x(2,k+1)] = f(x(1,k),x(2,k),u_actual(k),dt);
end

% get cost at every time step for plotting
V_plot = NaN*ones(1,Nt);
for k = 1:Nt
    V_plot(k) = interpn(x1_state_grid,x2_state_grid,V(:,:,k),x(1,k),x(2,k));
end

fig = figure;

subplot(311)
plot(t,x(1,:)); hold on;
plot(t,x(2,:));
legend('x1','x2')
ylabel('state','FontSize',16)

subplot(312)
plot(t(1:end-1),u_actual)
ylabel('u','FontSize',16)

subplot(313)
plot(t,V_plot)
xlabel('Time (s)','FontSize',16)
ylabel('Cost','FontSize',16)

% dynamics
function [x1_,x2_] = f(x1,x2,u,dt)
    x1_ = x1 + dt .* (-2.*(x1+0.25) + (x2+0.5).*exp(25.*x1./(x1+2)) - (x1+.25).*u);
    x2_ = x2 + dt .* (0.5 - x2 - (x2+0.5).*exp(25.*x1./(x1+2)));
end