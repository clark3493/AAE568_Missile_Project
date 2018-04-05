% This is the example starting on Slide 4 of lecture 27

n = 4;      % number of states
T = 3;      % final time increment

% costs
% u       1    2    3    4    states
gxu = [ Inf,   2,   6, Inf;   % 1
        Inf, Inf,   3,   8;   % 2
        Inf,   1, Inf,   4;   % 3
        Inf, Inf, Inf,   0 ]; % 4
    
% final costs
% state   1    2    3    4
gT  = [ Inf, Inf, Inf,   0 ]';

% allocate space for value function and control
V = [ Inf*ones(n,T), gT ];
u = NaN*zeros(n,T);

for t = T:-1:1
    cost = gxu + repmat(V(:,t+1)',n,1)
    [V(:,t),u(:,t)] = min(gxu + repmat(V(:,t+1)',n,1), [], 2);
    a = repmat(V(:,t+1)',n,1)
    [V(:,t),u(:,t)]
end

x = [1; NaN*ones(T-1,1)];
for t = 1:T
    x(t+1) = u(x(t),t);
end

x