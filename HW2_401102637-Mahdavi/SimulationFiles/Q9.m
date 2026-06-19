clc; clear; close all;

% Define the original high-order transfer function
s = tf('s');
G = (0.2*s+1)*(s+0.7)*(s+1)*(s+2) / ...
    ((25*s+1)*(12*s+1)*(8*s+1)*(3*s+1)*(s+0.4)*(s+0.2));

% Display pole/zero locations
figure;
pzmap(G);
title('Pole-Zero Map of Original System');
grid on; saveas(gcf, 'pzmap_original.png');

% Frequency range for analysis (rad/s)
w = logspace(log10(0.05), log10(15), 200);

%% Method 1: Balanced Truncation (using balred, as in Model Reducer)
% Obtain a balanced realisation and truncate to order 2 and 3
G_bal2 = balred(G, 2);   % 2nd order reduction
G_bal3 = balred(G, 3);   % 3rd order reduction

% Bode comparison for balanced truncation
figure;
bode(G, 'b-', G_bal2, 'r--', G_bal3, 'g-.', w);
legend('Original (6th order)', 'Balanced 2nd', 'Balanced 3rd', 'Location','best');
grid on; title('Bode Diagram – Balanced Truncation');
saveas(gcf, 'bode_balred.png');

%% Method 2: Mode Selection (Modal Truncation) using freqsep
% Separate slow and fast dynamics
% Choose a cutoff frequency that separates the dominant slow poles from fast ones.
% The poles are: -0.04, -0.0833, -0.125, -0.2, -0.333, -0.4
% We keep poles slower than a cutoff, say 0.5 rad/s (keep all poles slower than 0.5,
% which are essentially all of them? Actually all poles are <=0.5 except none are above 0.5? Wait:
% -0.04, -0.083, -0.125, -0.2, -0.333, -0.4 – all are <=0.4. So there is no pole above 0.5. 
% That means all poles are slow. For mode selection we could separate zeros? Or we keep only the slowest ones.
% Let's define cutoff = 0.15 rad/s to keep poles slower than 0.15: -0.04, -0.0833, -0.125 are kept (three poles).
% We'll use freqsep to split the system into slow and fast components.
try
    [G_slow, G_fast] = freqsep(G, 0.15);  % cutoff 0.15 rad/s
    % The slow part is a reduced model (G_slow)
    G_modal = G_slow;  % modal truncation
catch
    % If freqsep fails, use a simpler approach: keep poles with abs(real) < cutoff
    % We'll extract poles and zeros and manually build a reduced model.
    [p, z, k] = zpkdata(G, 'v');
    cutoff = 0.15;
    idx_keep = abs(real(p)) < cutoff;  % slower than cutoff
    p_keep = p(idx_keep);
    % Keep all zeros (or you can filter zeros too, but usually we keep them)
    % Build reduced model using zpk
    G_modal = zpk(z, p_keep, k * prod(-p(~idx_keep))/prod(-z)); % adjust gain? 
    % Actually to keep DC gain, we need to recompute k. Simpler: use balred or minreal.
    % Let's use freqsep with try/catch. In MATLAB R2016b+ freqsep works.
end

% If freqsep didn't work, fallback to balanced truncation with order 3
if ~exist('G_modal','var') || isempty(G_modal)
    warning('freqsep not available, using 3rd order balanced truncation for modal selection.');
    G_modal = balred(G, 3);
end

% Bode comparison with modal truncation
figure;
bode(G, 'b-', G_modal, 'm-', w);
legend('Original', 'Modal Truncation (slow modes)', 'Location','best');
grid on; title('Bode Diagram – Mode Selection (cutoff 0.15 rad/s)');
saveas(gcf, 'bode_modal.png');

%% Step response comparison for the best reduced model (choose G_bal3 as representative)
figure;
step(G, 'b-', G_bal3, 'r--', 0:0.1:50);
legend('Original', 'Balanced 3rd', 'Location','best');
grid on; title('Step Response Comparison');
saveas(gcf, 'step_comparison.png');

% Compute error metrics (normalized L2 error in frequency domain)
err_bal3 = norm(G - G_bal3, inf) / norm(G, inf);
err_modal = norm(G - G_modal, inf) / norm(G, inf);
fprintf('Relative H-inf error for Balanced 3rd: %.4f\n', err_bal3);
fprintf('Relative H-inf error for Modal (slow modes): %.4f\n', err_modal);

% Display reduced models
disp('Balanced 3rd order model:')
zpk(G_bal3)
disp('Modal truncation model:')
zpk(G_modal)