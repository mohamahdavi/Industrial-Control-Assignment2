clc; clear; close all;

%% System parameters
Ts = 0.01;           % sampling time for discrete simulation (Euler)
t = 0:Ts:10;         % simulation time vector (seconds)
N = length(t);

% Plant: P(s) = 1/(s+1)
% Discretised using forward Euler: y(k+1) = y(k) + Ts*(-y(k) + u(k))

Kp = 5;
Ki = 15;

% Actuator saturation limits
u_max = 2;
u_min = -2;

% Reference: unit step at t = 1 s
r = (t >= 1) * 1;   % 0 for t<1, 1 for t>=1

%% ===== (a) No anti-windup (pure integrator) =====
y_a = zeros(size(t));        % plant output
xI_a = zeros(size(t));       % integrator state
u_c_a = zeros(size(t));      % controller output before saturation
u_a = zeros(size(t));        % actual plant input (saturated)
e_a = zeros(size(t));        % tracking error

y_a(1) = 0; xI_a(1) = 0;
for k = 1:N-1
    e_a(k) = r(k) - y_a(k);
    u_c_a(k) = Kp*e_a(k) + Ki*xI_a(k);        % PI control law
    u_a(k) = min(max(u_c_a(k), u_min), u_max); % saturation
    % Integrator state update (standard, no anti-windup)
    xI_a(k+1) = xI_a(k) + Ts * e_a(k);
    % Plant dynamics (Euler)
    y_a(k+1) = y_a(k) + Ts * (-y_a(k) + u_a(k));
end
% Last sample evaluation
e_a(N) = r(N) - y_a(N);
u_c_a(N) = Kp*e_a(N) + Ki*xI_a(N);
u_a(N) = min(max(u_c_a(N), u_min), u_max);

%% Plot results for (a)
figure;
subplot(3,1,1);
plot(t, y_a, 'b', 'LineWidth',1.5);
ylabel('y(t)'); title('(a) Step response – No Anti-Windup'); grid on;

subplot(3,1,2);
plot(t, u_c_a, 'r', t, u_a, 'b--', 'LineWidth',1.5);
ylabel('Control signals'); legend('u_c','u','Location','best'); grid on;

subplot(3,1,3);
plot(t, e_a, 'k', 'LineWidth',1.5);
ylabel('e(t)'); xlabel('Time (s)'); grid on;
sgtitle('Windup behaviour without anti-windup');
saveas(gcf, 'a_windup.png');

%% ===== (b) Back-calculation anti-windup =====
kaw_values = [5, 20];   % two different anti-windup gains
colors = {'b','r'};
y_b_all = cell(length(kaw_values),1);  % store outputs for later comparison
u_b_all = cell(length(kaw_values),1);

figure;
for idx = 1:length(kaw_values)
    kaw = kaw_values(idx);
    y_b = zeros(size(t));
    xI_b = zeros(size(t));
    u_c_b = zeros(size(t));
    u_b = zeros(size(t));
    e_b = zeros(size(t));
    y_b(1)=0; xI_b(1)=0;

    for k = 1:N-1
        e_b(k) = r(k) - y_b(k);
        u_c_b(k) = Kp*e_b(k) + Ki*xI_b(k);
        u_b(k) = min(max(u_c_b(k), u_min), u_max);
        % Back-calculation: integrator update with (u - u_c) feedback
        xI_b(k+1) = xI_b(k) + Ts * ( e_b(k) + kaw*(u_b(k) - u_c_b(k)) );
        y_b(k+1) = y_b(k) + Ts * (-y_b(k) + u_b(k));
    end
    y_b_all{idx} = y_b;
    u_b_all{idx} = u_b;

    % Plot output and control for this kaw
    subplot(2,1,1);
    hold on;
    plot(t, y_b, colors{idx}, 'LineWidth',1.5);
    grid on;
    subplot(2,1,2);
    hold on;
    plot(t, u_b, colors{idx}, 'LineWidth',1.5);
    grid on;
end
subplot(2,1,1);
xlabel('Time (s)'); ylabel('y(t)');
title('(b) Back-calculation anti-windup: Output');
legend(strcat('k_{aw}=', string(kaw_values)), 'Location','best');

subplot(2,1,2);
xlabel('Time (s)'); ylabel('u(t)');
title('Control signal (saturated)');
legend(strcat('k_{aw}=', string(kaw_values)), 'Location','best');
saveas(gcf, 'b_backcalc.png');

% Store y_b for kaw=20 for later comparison with Simulink (will be reused)
y_b_kaw20 = y_b_all{2};

%% ===== (c) Conditional integration =====
y_c = zeros(size(t));
xI_c = zeros(size(t));
u_c_c = zeros(size(t));
u_c = zeros(size(t));
e_c = zeros(size(t));
y_c(1)=0; xI_c(1)=0;

for k = 1:N-1
    e_c(k) = r(k) - y_c(k);
    u_c_c(k) = Kp*e_c(k) + Ki*xI_c(k);
    u_c_sat = min(max(u_c_c(k), u_min), u_max);
    u_c(k) = u_c_sat;

    % Conditional integration logic:
    % If saturated at upper limit and error positive, OR
    % saturated at lower limit and error negative, stop integrating.
    if (u_c_sat == u_max && e_c(k) > 0) || (u_c_sat == u_min && e_c(k) < 0)
        dxI = 0;
    else
        dxI = e_c(k);
    end
    xI_c(k+1) = xI_c(k) + Ts * dxI;
    y_c(k+1) = y_c(k) + Ts * (-y_c(k) + u_c(k));
end

% Plot comparison with no anti-windup
figure;
subplot(2,1,1);
plot(t, y_a, 'b', t, y_c, 'r--', 'LineWidth',1.5); grid on;
ylabel('y(t)'); legend('No anti-windup','Conditional Integration','Location','best');
title('(c) Conditional integration vs. no anti-windup');
subplot(2,1,2);
plot(t, u_a, 'b', t, u_c, 'r--', 'LineWidth',1.5); grid on;
ylabel('u(t)'); xlabel('Time (s)'); legend('No anti-windup','Conditional Integration','Location','best');
title('Control signal');
saveas(gcf, 'c_conditional.png');

%% ===== (d) Simulink PID block with built-in anti-windup =====
modelName = 'PID_antiwindup_simulink';
if bdIsLoaded(modelName)
    close_system(modelName,0);
end
if exist([modelName '.slx'], 'file')
    delete([modelName '.slx']);
end

open_system(new_system(modelName));

% Blocks with correct library paths and parameter names
add_block('simulink/Sources/Step', [modelName '/Step'], ...
    'Time','1', 'Before','0', 'After','1');
add_block('simulink/Continuous/PID Controller', [modelName '/PID'], ...
    'P','5', 'I','15', 'D','0', ...
    'LimitOutput','on', ...
    'UpperSaturationLimit','2', 'LowerSaturationLimit','-2', ...
    'AntiWindupMode','back-calculation', 'Kb','1');   % corrected parameter
add_block('simulink/Continuous/Transfer Fcn', [modelName '/Plant'], ...
    'Numerator','[1]', 'Denominator','[1 1]');
add_block('simulink/Sinks/Scope', [modelName '/Scope']);
add_block('simulink/Sinks/To Workspace', [modelName '/y_out'], ...
    'VariableName','y_sim', 'SaveFormat','Array');
add_block('simulink/Sinks/To Workspace', [modelName '/u_out'], ...
    'VariableName','u_sim', 'SaveFormat','Array');
add_block('simulink/Sinks/To Workspace', [modelName '/t_out'], ...
    'VariableName','t_sim', 'SaveFormat','Array');
add_block('simulink/Sources/Clock', [modelName '/Clock']);

% Connections
add_line(modelName, 'Step/1', 'PID/1');
add_line(modelName, 'PID/1', 'Plant/1');
add_line(modelName, 'Plant/1', 'Scope/1');
add_line(modelName, 'Plant/1', 'y_out/1');
add_line(modelName, 'PID/1', 'u_out/1');
add_line(modelName, 'Clock/1', 't_out/1');

% Solver settings
set_param(modelName, 'StopTime', '10', ...
    'SolverType', 'Fixed-step', 'Solver', 'ode4', 'FixedStep', '0.01');
save_system(modelName);
out = sim(modelName);
t_sim = out.tout;
y_sim = out.y_sim;
u_sim = out.u_sim;

% Plot comparison
figure;
plot(t_sim, y_sim, 'b-', t, y_b_kaw20, 'r--', 'LineWidth',1.5);
grid on;
xlabel('Time (s)'); ylabel('y(t)');
legend('Simulink built-in anti-windup (back-calc, Kb=1)', ...
    'Manual back-calculation (k_{aw}=20)', 'Location','best');
title('(d) Simulink built-in anti-windup vs manual method');
saveas(gcf, 'd_simulink.png');

close_system(modelName,0);