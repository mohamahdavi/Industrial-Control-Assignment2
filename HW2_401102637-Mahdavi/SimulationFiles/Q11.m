clc; clear; close all;

%% Plant model
s = tf('s');
numG = 2;
denG = conv([4 1], [3 1]);   % (4s+1)(3s+1) = 12s^2 + 7s + 1
G = tf(numG, denG, 'InputDelay', 0.3);

%% Simulation parameters
Ts = 0.01;
Tfinal = 30;
t = 0:Ts:Tfinal;
r = ones(size(t));            % unit step reference

% Control signal bounds
u_max = 5; u_min = -5;

%% Objective function (self-contained, with all needed variables defined inside)
function J = objectiveFunction(K)
    % K = [Kp, Ki, Kd, Tf]
    Kp = K(1); Ki = K(2); Kd = K(3); Tf = K(4);
    if Tf <= 0
        J = 1e6; return;
    end

    % Simulation parameters
    Ts = 0.01;
    Tfinal = 30;
    t = (0:Ts:Tfinal)';
    r = ones(size(t));

    % Control bounds
    u_max = 5;

    % Plant
    s = tf('s');
    G = tf(2, conv([4 1], [3 1]), 'InputDelay', 0.3);

    % PIDF controller
    numC = [Kp*Tf + Kd,  Kp + Ki*Tf,  Ki];
    denC = [Tf, 1, 0];
    C = tf(numC, denC);

    % Closed-loop system
    CL = feedback(C * G, 1);

    % Step response
    [y, t_sim] = step(CL, t);
    if length(y) < length(t_sim)
        J = 1e6; return;
    end

    % Error
    e = r - y;

    % Control signal (lsim from error)
    [u, ~] = lsim(C, e, t_sim);

    % Performance metrics
    ISE = sum(e.^2) * Ts;
    ITSE = sum((t_sim(:) .* e(:)).^2) * Ts;
    control_effort = sum(u.^2) * Ts;

    % Constraint violation on u
    u_viol = sum(max(0, abs(u) - u_max).^2) * Ts;

    % Rise time (10% → 90%)
    idx_10 = find(y >= 0.1, 1);
    idx_90 = find(y >= 0.9, 1);
    if isempty(idx_10) || isempty(idx_90)
        rise_time = Tfinal;
    else
        rise_time = t_sim(idx_90) - t_sim(idx_10);
    end

    % Overshoot
    yss = dcgain(CL);
    if isnan(yss), yss = 1; end
    overshoot = max(0, max(y) - yss) / yss;

    % Overshoot penalty (>5%)
    overshoot_penalty = 0;
    if overshoot > 0.05
        overshoot_penalty = 100 * (overshoot - 0.05)^2;
    end

    % Steady-state error penalty (>1%)
    ess = abs(1 - yss);
    if ess > 0.01
        ess_penalty = 1000 * (ess - 0.01)^2;
    else
        ess_penalty = 0;
    end

    % Weights
    w1 = 1;        % ITSE
    w2 = 0.1;      % control effort
    w3 = 100;      % constraint violation
    w4 = 10;       % rise time
    w5 = 1000;     % overshoot penalty
    w6 = 1000;     % steady-state error penalty

    J = w1*ITSE + w2*control_effort + w3*u_viol + ...
        w4*rise_time + w5*overshoot_penalty + w6*ess_penalty;
end

%% Bounds for [Kp, Ki, Kd, Tf]
lb = [0, 0, 0, 0.001];
ub = [10, 10, 10, 0.5];

%% --- Genetic Algorithm ---
options_ga = optimoptions('ga', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 50, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotbestf);
fprintf('Running Genetic Algorithm...\n');
[K_ga, J_ga] = ga(@objectiveFunction, 4, [], [], [], [], lb, ub, [], options_ga);
fprintf('GA optimal: Kp=%.4f, Ki=%.4f, Kd=%.4f, Tf=%.4f, Cost=%.4f\n', ...
    K_ga(1), K_ga(2), K_ga(3), K_ga(4), J_ga);

%% --- Particle Swarm ---
options_ps = optimoptions('particleswarm', ...
    'SwarmSize', 50, ...
    'MaxIterations', 50, ...
    'Display', 'iter');
fprintf('Running Particle Swarm...\n');
[K_ps, J_ps] = particleswarm(@objectiveFunction, 4, lb, ub, options_ps);
fprintf('PS optimal: Kp=%.4f, Ki=%.4f, Kd=%.4f, Tf=%.4f, Cost=%.4f\n', ...
    K_ps(1), K_ps(2), K_ps(3), K_ps(4), J_ps);

%% --- fmincon (interior‑point) ---
x0 = [1, 0.5, 0.1, 0.05];   % initial guess
options_fmincon = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'iter');
fprintf('Running fmincon...\n');
[K_fc, J_fc] = fmincon(@objectiveFunction, x0, [], [], [], [], lb, ub, [], options_fmincon);
fprintf('fmincon optimal: Kp=%.4f, Ki=%.4f, Kd=%.4f, Tf=%.4f, Cost=%.4f\n', ...
    K_fc(1), K_fc(2), K_fc(3), K_fc(4), J_fc);

%% --- Simulated Annealing ---
options_sa = optimoptions('simulannealbnd', ...
    'InitialTemperature', 100, ...
    'MaxIterations', 200, ...
    'Display', 'iter');
fprintf('Running Simulated Annealing...\n');
[K_sa, J_sa] = simulannealbnd(@objectiveFunction, x0, lb, ub, options_sa);
fprintf('SA optimal: Kp=%.4f, Ki=%.4f, Kd=%.4f, Tf=%.4f, Cost=%.4f\n', ...
    K_sa(1), K_sa(2), K_sa(3), K_sa(4), J_sa);

%% Compare step responses using best solution (use GA result)
K_opt = K_ga;   % or choose the best overall
Kp_opt = K_opt(1); Ki_opt = K_opt(2); Kd_opt = K_opt(3); Tf_opt = K_opt(4);
numC_opt = [Kp_opt*Tf_opt + Kd_opt,  Kp_opt + Ki_opt*Tf_opt,  Ki_opt];
denC_opt = [Tf_opt, 1, 0];
C_opt = tf(numC_opt, denC_opt);
CL_opt = feedback(C_opt * G, 1);

figure;
step(CL_opt, t); grid on;
title('Step Response with GA‑Optimized PIDF Controller');
xlabel('Time (s)'); ylabel('Output');
saveas(gcf, 'pidf_step_ga.png');

% Disturbance response: add step disturbance 0.5 at t=10 s at plant input
t_dist = 0:Ts:Tfinal;
d = (t_dist >= 10) * 0.5;
% Transfer from disturbance to output: S = G/(1 + C_opt*G)
S_dist = feedback(G, C_opt);          % disturbance enters at plant input
y_dist = lsim(S_dist, d, t_dist);
figure;
plot(t_dist, y_dist); grid on;
title('Response to Step Disturbance at Plant Input');
xlabel('Time (s)'); ylabel('Output');
saveas(gcf, 'pidf_disturbance.png');

% Robustness: vary time constant and delay ±20%
% Nominal: G_nom = 2*exp(-0.3s)/((4s+1)(3s+1))
pert_vals = [-20, 0, 20];  % percent
figure; hold on;
for p = pert_vals
    tau_p = 0.3 * (1 + p/100);
    G_pert = tf(2, conv([4 1], [3 1]), 'InputDelay', tau_p);
    CL_pert = feedback(C_opt * G_pert, 1);
    step(CL_pert, t);
end
legend('-20%', 'Nominal', '+20%');
title('Robustness to Delay Variation');
grid on; saveas(gcf, 'pidf_robustness_delay.png');

% Similarly for time constant variation (vary both time constants by same factor)
figure; hold on;
for p = pert_vals
    scale = 1 + p/100;
    G_pert = tf(2, conv([4*scale 1], [3*scale 1]), 'InputDelay', 0.3);
    CL_pert = feedback(C_opt * G_pert, 1);
    step(CL_pert, t);
end
legend('-20%', 'Nominal', '+20%');
title('Robustness to Time‑Constant Variation');
grid on; saveas(gcf, 'pidf_robustness_tc.png');

%% Effect of population size and generations (GA)
pop_sizes = [20, 50, 100];
gens = [30, 50, 80];
results = [];
for ps = pop_sizes
    for gn = gens
        opts = optimoptions('ga', 'PopulationSize', ps, 'MaxGenerations', gn, 'Display', 'off');
        [K_tmp, J_tmp] = ga(@objectiveFunction, 4, [], [], [], [], lb, ub, [], opts);
        results = [results; ps, gn, K_tmp, J_tmp];
    end
end
disp('Population size / Generations / Cost:');
for i = 1:size(results,1)
    fprintf('Pop=%d, Gen=%d, Cost=%.4f\n', results(i,1), results(i,2), results(i,end));
end

