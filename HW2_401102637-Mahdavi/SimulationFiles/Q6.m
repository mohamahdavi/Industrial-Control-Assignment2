clc; clear; close all;
s = tf('s');

%% Plant and desired polynomial
G = 1/(s+1);
% Desired closed-loop poles: -2±2j, extra at -5
% Δ_d(s) = (s+2-2j)(s+2+2j)(s+5)
Delta_d = (s+2-2i)*(s+2+2i)*(s+5);

%% Derived controller from part (c)
% Using augmented plant A'(s)=s^2*(s+1), B'(s)=1
% Solved: L(s)=1, P(s)=8s^2+28s+40
% Overall controller C(s) = P(s) / s^2
C = (8*s^2 + 28*s + 40) / s^2;

% Check characteristic polynomial
char_poly = minreal(1 + G*C); % should be s^3+9s^2+28s+40 / s^2(s+1)
% Get closed-loop denominator
[num_cl, den_cl] = tfdata(feedback(G*C,1), 'v');
roots(den_cl)   % display closed-loop poles

%% (f) Build closed-loop system and confirm poles
T = feedback(G*C, 1);
poles = pole(T);
fprintf('Closed-loop poles:\n');
disp(poles);
% Compare with roots of Δ_d(s)
desired_roots = roots([1 9 28 40]);
fprintf('Desired characteristic polynomial roots:\n');
disp(desired_roots);

%% (g) Ramp response with designed controller
t = 0:0.01:10;
r = t;                     % unit ramp
[y, ~] = lsim(T, r, t);
e = r - y;                 % tracking error

figure;
subplot(2,1,1);
plot(t, y, 'b', t, r, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Output');
legend('y(t)', 'r(t)=t', 'Location', 'northwest');
title('Ramp response with pole-placement controller');
grid on;

subplot(2,1,2);
plot(t, e, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Tracking error');
title('Tracking error e(t) = r(t)-y(t)');
grid on;
saveas(gcf, 'ramp_tracking.png');

%% (h) Ramp response with proportional controller
Kp = 1;
T_p = feedback(G*Kp, 1);
[y_p, ~] = lsim(T_p, r, t);
e_p = r - y_p;

figure;
subplot(2,1,1);
plot(t, y_p, 'b', t, r, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Output');
legend('y(t)', 'r(t)=t', 'Location', 'northwest');
title('Ramp response with proportional controller (K=1)');
grid on;

subplot(2,1,2);
plot(t, e_p, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Tracking error');
title('Tracking error (proportional)');
grid on;
saveas(gcf, 'ramp_tracking_proportional.png');

%% (i) Ramp disturbance at plant input (starting at t=3s)
% Disturbance signal d(t) = (t-3)*1(t-3)
d = (t-3).*((t-3)>=0);
% Transfer function from disturbance at plant input to output: S(s) = G/(1+G*C)
S_dist = feedback(G, C);   % sensitivity = G/(1+G*C) (if disturbance enters at plant input)
% Actually disturbance enters at input of plant: y = G*(u+d), u = C*(r-y) -> y = G*C*(r-y) + G*d => y*(1+G*C)=G*C*r+G*d => y = T*r + S_dist*d where S_dist = G/(1+G*C). Yes.

% Simulate disturbance response with zero reference
y_dist = lsim(S_dist, d, t);

figure;
plot(t, y_dist, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Output');
title('Response to ramp disturbance at plant input (starting at t=3s)');
grid on;
saveas(gcf, 'disturbance_response.png');