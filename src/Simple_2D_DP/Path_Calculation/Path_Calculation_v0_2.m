%% Optimal Path Calculation
%  AAE 568 Semester Project
%  DP Method
%
%  Uses csv sheets produced from C code
%  R L Robinson
%
%  v0_2 updates the state grids from 
%% Initialize params
%  Trial will be done for tf4
close all, clear all

tic

tau_size = 100;
% flight_Path = [X Z Xdot Zdot InputU time Vel_Angle]
flight_Path = zeros(tau_size, 6);

X0 = 0; Z0 = -2001; Xdot0 = 100; Zdot0 = 0; initial_U = -77.5439463761163;
tf = 17.678714;             %************MUST BE ASSIGNED DYNAMICALLY*****
dt = tf/(tau_size);
T = 1000; m = 100; g = 9.81;
X_Target = 2000;

% Get discretized parameters
X_min = -100; X_max = 3500;
Z_min = -2101; Z_max = 1499;
Xdot_min = 75; Xdot_max = 318.172996;
Zdot_min = -100; Zdot_max = 318.172996;

% NOTE, unless rational, these will not match exactly during interpolation
X_dis = linspace(X_min,X_max, 35);
Z_dis = linspace(Z_min, Z_max,35);
Xdot_dis = linspace(Xdot_min, Xdot_max, 21); 
Zdot_dis = linspace(Zdot_min, Zdot_max, 21);


flight_Path(1, 1) = X0; flight_Path(1, 2) = Z0;
flight_Path(1, 3) = Xdot0; flight_Path(1, 4) = Zdot0;
flight_Path(1, 5) = initial_U; flight_Path(1, 6) = 0;
flight_Path(1, 7) = atan2d(flight_Path(1,4), flight_Path(1, 3));
%% Calculate next state and do interpolation for next u
% Current State
for i = 1:99
    Xc = flight_Path(i, 1);         % c denotes current
    Zc = flight_Path(i, 2);
    Xdotc = flight_Path(i, 3);
    Zdotc = flight_Path(i, 4);
    u_deg = flight_Path(i, 5);

    % Next State
    X2 = Xc + dt * Xdotc;
    Z2 = Zc + dt * Zdotc;
    Xdot2 = Xdotc + dt * (  T/m * cosd(u_deg) );
    Zdot2 = Zdotc + dt * ( -T/m * sind(u_deg) + g );
    %% Interpolation
    % Get the next u value
    % Open up the next time step and load in the state matrix
    main_path = 'C:\Users\Reese Robinson\Desktop\Purdue\AAE 568 Applied Optimal Control and Estimation\Project\C DP Program\Local_Comp\A1_B0_45\data\';
    tf_value = 'tf4';
    tau_value = strcat('tau', num2str(i));
    full_path = strcat(main_path,tf_value,'\cost_',tf_value, '_',tau_value,'.csv');

    % state_Matrix [X Z Xdot Zdot Cost U]
    state_Matrix = csvread(full_path, 1);

    % Get upper and lower bounds, for interpolation
    clear X0 Z0 Xdot0 Zdot0 N    % Repurposing
    for j = 1:35
        if X_dis(j) < X2
            X0 = X_dis(j);
        end

        if Z_dis(j) < Z2
            Z0 = Z_dis(j);
        end
        
        if j < 22 % dot values discretized to 21 values
            if Xdot_dis(j) < Xdot2
                Xdot0 = Xdot_dis(j);
            end

            if Zdot_dis(j) < Zdot2
                Zdot0 = Zdot_dis(j);
            end
        end
    end

    for j = 35:-1.:1
        if X_dis(j) > X2
            X1 = X_dis(j);
        end

        if Z_dis(j) > Z2
            Z1 = Z_dis(j);
        end
        
        if j < 22 % dot values discretized to 21 values
            if Xdot_dis(j) > Xdot2
                Xdot1 = Xdot_dis(j);
            end

            if Zdot_dis(j) > Zdot2
                Zdot1 = Zdot_dis(j);
            end
        end
    end

    %Calculate the N parameter, to be multiplied agains the control input
    denom = (X1-X0)*(Z1-Z0)*(Xdot1-Xdot0)*(Zdot1-Zdot0);

    X1X2 = X1 - X2;
    X2X0 = X2 - X0;
    Z1Z2 = Z1 - Z2;
    Z2Z0 = Z2 - Z0;
    Xdot12 = Xdot1 - Xdot2;
    Xdot20 = Xdot2 - Xdot0;
    Zdot12 = Zdot1 - Zdot2;
    Zdot20 = Zdot2 - Zdot0;

    N = ones(16, 1);
    N(16) = X1X2*Z1Z2*Xdot12*Zdot12/denom;
    N(15) = X1X2*Z1Z2*Xdot12*Zdot20/denom;
    N(14) = X1X2*Z1Z2*Xdot20*Zdot12/denom;
    N(13) = X1X2*Z1Z2*Xdot20*Zdot20/denom;
    N(12) = X1X2*Z2Z0*Xdot12*Zdot12/denom;
    N(11) = X1X2*Z2Z0*Xdot12*Zdot20/denom;
    N(10) = X1X2*Z2Z0*Xdot20*Zdot12/denom;
    N(9) = X1X2*Z2Z0*Xdot20*Zdot20/denom;
    N(8) = X2X0*Z1Z2*Xdot12*Zdot12/denom;
    N(7) = X2X0*Z1Z2*Xdot12*Zdot20/denom;
    N(6) = X2X0*Z1Z2*Xdot20*Zdot12/denom;
    N(5) = X2X0*Z1Z2*Xdot20*Zdot20/denom;
    N(4) = X2X0*Z2Z0*Xdot12*Zdot12/denom;
    N(3) = X2X0*Z2Z0*Xdot12*Zdot20/denom;
    N(2) = X2X0*Z2Z0*Xdot20*Zdot12/denom;
    N(1) = X2X0*Z2Z0*Xdot20*Zdot20/denom;

    % Get u values from state_matrix
    U = zeros(16, 1);
    state_Length = size(state_Matrix, 1);
    for k = 1:state_Length
        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(1) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(2) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(3) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(4) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(5) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(6) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(7) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X0) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(8) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(9) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(10) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(11) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z0) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(12) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(13) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot0) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(14) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot0) < 1
           U(15) = state_Matrix(k, 6);
        end

        if abs(state_Matrix(k, 1) - X1) < 1 && ...
                abs(state_Matrix(k, 2) - Z1) < 1 && ...
                abs(state_Matrix(k, 3) - Xdot1) < 1 && ...
                abs(state_Matrix(k, 4) - Zdot1) < 1
           U(16) = state_Matrix(k, 6);
        end
    end

    % Interpolate
    nextU = 0;
    for k = 1:16
       nextU = nextU + N(k) * U(k); 
    end

    %% Record params in flight path
    flight_Path(i + 1, 1) = X2;
    flight_Path(i + 1, 2) = Z2;
    flight_Path(i + 1, 3) = Xdot2;
    flight_Path(i + 1, 4) = Zdot2;
    flight_Path(i + 1, 5) = nextU;
    flight_Path(i + 1, 6) = i*dt;
    flight_Path(i + 1, 7) = atan2d(flight_Path(i+1,4), flight_Path(i+1, 3));
    
    % if the missile has reached Z = 0 it has impacted. Break Loop
    if flight_Path(i, 2) > 0
       break
    end
end

%% Plot/Conclude everything
impact_Angle = atan2d(flight_Path(i, 4), flight_Path(i, 3));
flight_time = flight_Path(i, 6);
fprintf('Missile Impact Angle: %f degrees\n', impact_Angle)
fprintf('Flight time: %f sec \n', flight_time)
fprintf('Missile impacted %f meters away from target\n',...
    X_Target - flight_Path(i, 1))

% redefine flight_Path in order to plot it - get rid of the empty rows
flight_Path = flight_Path(1:i, :);

fig = figure;

subplot(311)
plot(flight_Path(:,1), -flight_Path(:,2))
xlabel('x (m)','FontSize',16')
ylabel('Altitude (m)','FontSize',16)
grid on;

subplot(312)
plot(flight_Path(:, 6), flight_Path(:, 5))
xlabel('time (s)','FontSize',16)
ylabel('$\theta$ (deg)','FontSize',16,'Interpreter','latex')
grid on;

subplot(313)
plot(flight_Path(:, 6), flight_Path(:, 7))
xlabel('time (s)','FontSize',16)
ylabel('Flight direction (deg)','FontSize',16)
grid on;

toc