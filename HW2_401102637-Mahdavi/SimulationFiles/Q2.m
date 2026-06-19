clc; clear; close all;

%% Original system
s = tf('s');
G = 1 / (s+1)^8;

%% 1. Step response of the original system and plot
t_final = 50;               % simulation time (until steady state)
[y, t] = step(G, t_final);

figure('Name', 'Step Response of Original System');
plot(t, y, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)'); ylabel('Output');
title('Step Response of G(s) = 1/(s+1)^8');
saveas(gcf, 'original_step.png');   % save figure

%% 2. Reduced-order models (using different methods)

% --- (a) First-order model: 63.2% method ---
idx_a = find(y >= 0.632, 1);        % first time output reaches 63.2%
T_a = t(idx_a);
G_a = 1 / (T_a*s + 1);              % gain = 1
fprintf('(a) First-order: T = %.4f s\n', T_a);

% --- (b) FOPDT model: tangent method (inflection point) ---
dy = gradient(y, t);
[~, idx_max] = max(dy);
t_inf = t(idx_max);
y_inf = y(idx_max);
slope = dy(idx_max);

L_b = t_inf - y_inf / slope;        % intersection with y=0
if L_b < 0, L_b = 0; end
t_K = t_inf + (1 - y_inf) / slope; % intersection with final value (1)
T_b = t_K - L_b;

G_b = tf(1, [T_b 1], 'InputDelay', L_b); % gain = 1
fprintf('(b) FOPDT: T = %.4f s, L = %.4f s\n', T_b, L_b);

% --- (c) SOPDT model: numerical optimization (minimizing step response error) ---
% Model: G_m(s) = 1 / ((1+T1*s)*(1+T2*s)) * exp(-L*s)
cost_fun = @(x) norm( step( tf(1, conv([x(1) 1], [x(2) 1]), 'InputDelay', x(3)), t ) - y );
x0 = [3, 4, 1];                     % initial guess
x_opt = fminsearch(cost_fun, x0);
T1_c = x_opt(1);
T2_c = x_opt(2);
L_c  = x_opt(3);

G_c = tf(1, conv([T1_c 1], [T2_c 1]), 'InputDelay', L_c);
fprintf('(c) SOPDT: T1 = %.4f s, T2 = %.4f s, L = %.4f s\n', T1_c, T2_c, L_c);

%% 3. Comparison of all responses
figure('Name', 'Comparison of Reduced Models');
step(G, 'b-', G_a, 'g--', G_b, 'r-.', G_c, 'm:', t);
legend('Original 1/(s+1)^8', ...
       sprintf('1st-order (T=%.2f)', T_a), ...
       sprintf('FOPDT (T=%.2f, L=%.2f)', T_b, L_b), ...
       sprintf('SOPDT (T1=%.2f, T2=%.2f, L=%.2f)', T1_c, T2_c, L_c), ...
       'Location', 'southeast');
grid on;
xlabel('Time (s)'); ylabel('Output');
title('Comparison of Reduced Models with Original System');
saveas(gcf, 'model_comparison.png'); % save figure

disp('--- Done. Check the figures. ---');