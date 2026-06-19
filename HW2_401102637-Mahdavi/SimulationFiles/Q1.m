clc; clear; close all;

%% Define the original transfer functions
s = tf('s');

% Plant G1 (stable)
G1 = 80/(s+8)^2 * (1/(s+2) + 0.4/(s+0.04));

% Plant G2 (integrating)
G2 = (s+5)^2 / (s * (s+2)^2 * (s+25));

%% 1. Fit reduced models from open‑loop step response
disp('--- Fitting reduced models from step response ---');

% --- G1 : FOPDT (stable process) ---
[K1, tau1, theta1] = fitFOPDT(G1, 30);
fprintf('G1 FOPDT: K = %.4f, tau = %.4f, theta = %.4f\n', K1, tau1, theta1);

% --- G2 : IPDT (integrating process) ---
[K2, theta2] = fitIPDT(G2, 50);
fprintf('G2 IPDT:  K = %.6f, theta = %.4f\n', K2, theta2);

%% 2. CHR PID tuning for set‑point tracking (0% and 20% overshoot)
disp('--- CHR PID gains ---');

% -------- G1 (FOPDT) --------
% 0% overshoot
Kc1_0  = 0.6 * tau1 / (K1 * theta1);
Ti1_0  = tau1;
Td1_0  = 0.5 * theta1;
% 20% overshoot
Kc1_20 = 0.95 * tau1 / (K1 * theta1);
Ti1_20 = 1.4 * tau1;
Td1_20 = 0.47 * theta1;

fprintf('G1 (0%% OS):  Kc = %.4f, Ti = %.4f, Td = %.4f\n', Kc1_0, Ti1_0, Td1_0);
fprintf('G1 (20%% OS): Kc = %.4f, Ti = %.4f, Td = %.4f\n', Kc1_20, Ti1_20, Td1_20);

% -------- G2 (IPDT) --------
% 0% overshoot
Kc2_0  = 0.6 / (K2 * theta2);
Ti2_0  = theta2;
Td2_0  = 0.5 * theta2;
% 20% overshoot
Kc2_20 = 0.95 / (K2 * theta2);
Ti2_20 = 1.4 * theta2;
Td2_20 = 0.47 * theta2;

fprintf('G2 (0%% OS):  Kc = %.4f, Ti = %.4f, Td = %.4f\n', Kc2_0, Ti2_0, Td2_0);
fprintf('G2 (20%% OS): Kc = %.4f, Ti = %.4f, Td = %.4f\n', Kc2_20, Ti2_20, Td2_20);

%% 3. Build PID controllers (ideal parallel form)
% PID(s) = Kc * (1 + 1/(Ti*s) + Td*s)

C1_0  = Kc1_0  * (1 + 1/(Ti1_0*s)  + Td1_0*s);
C1_20 = Kc1_20 * (1 + 1/(Ti1_20*s) + Td1_20*s);

C2_0  = Kc2_0  * (1 + 1/(Ti2_0*s)  + Td2_0*s);
C2_20 = Kc2_20 * (1 + 1/(Ti2_20*s) + Td2_20*s);

%% 4. Form closed‑loop systems (unity feedback)
CL1_0  = feedback(C1_0  * G1, 1);
CL1_20 = feedback(C1_20 * G1, 1);

CL2_0  = feedback(C2_0  * G2, 1);
CL2_20 = feedback(C2_20 * G2, 1);

%% 5. Simulate step responses and plot
t1 = 0:0.01:30;    % time vector for G1
t2 = 0:0.01:50;    % time vector for G2

% ---- Plot for G1 ----
figure('Name','G1 – CHR PID tracking');
[y1_0,  ~] = step(CL1_0,  t1);
[y1_20, ~] = step(CL1_20, t1);
plot(t1, y1_0, 'b-', 'LineWidth', 1.5); hold on;
plot(t1, y1_20, 'r--', 'LineWidth', 1.5);
yline(1, 'k:', 'LineWidth', 1.2);          % reference
xlabel('Time (s)'); ylabel('Output');
title('G1 – PID Set‑point Tracking (CHR method)');
legend('0% OS', '20% OS', 'Reference');
grid on;

% ---- Plot for G2 ----
figure('Name','G2 – CHR PID tracking');
[y2_0,  ~] = step(CL2_0,  t2);
[y2_20, ~] = step(CL2_20, t2);
plot(t2, y2_0, 'b-', 'LineWidth', 1.5); hold on;
plot(t2, y2_20, 'r--', 'LineWidth', 1.5);
yline(1, 'k:', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Output');
title('G2 – PID Set‑point Tracking (CHR method)');
legend('0% OS', '20% OS', 'Reference');
grid on;

disp('--- Done. Check the two figure windows. ---');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K, tau, theta] = fitFOPDT(G, tfinal)
    % Fit FOPDT model (K*exp(-theta*s)/(tau*s+1)) to stable process G
    % using the tangent method on the step response.
    [y, t] = step(G, tfinal);
    K = y(end);                          % steady‑state gain
    dy = gradient(y, t);                 % numerical derivative
    [~, idx] = max(dy);                  % inflection point index
    t_inf = t(idx);
    y_inf = y(idx);
    slope = dy(idx);                     % maximum slope

    % Intersection of tangent with initial level (y = 0)
    theta = t_inf - y_inf / slope;
    if theta < 0, theta = 0; end

    % Intersection with steady‑state line (y = K)
    t_K = t_inf + (K - y_inf) / slope;
    tau = t_K - theta;
    if tau <= 0, tau = 1e-3; end         % avoid division by zero
end

function [K, theta] = fitIPDT(G, tfinal)
    % Fit IPDT model (K*exp(-theta*s)/s) to integrating process G.
    % K is the velocity gain (final ramp slope), theta is the dead time.
    [y, t] = step(G, tfinal);
    dy = gradient(y, t);
    [~, idx] = max(dy);                  % inflection point (max slope)
    t_inf = t(idx);
    y_inf = y(idx);
    slope = dy(idx);                     % the ramp slope

    % Intersection of tangent with time axis (y = 0)
    theta = t_inf - y_inf / slope;
    if theta < 0, theta = 0; end

    K = slope;                           % velocity gain
end