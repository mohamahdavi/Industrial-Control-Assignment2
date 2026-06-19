clc; clear; close all;

%% General settings
Ts = 0.01;
T = 25;
t = 0:Ts:T;
N = length(t);

%% Part 1: Chirp input (with DC offset to ensure non-zero mean)
% Generate swept-frequency cosine signal + DC offset
f0 = 0.1;       % Hz at t=0
f1 = 50;        % Hz at t=T
u_chirp = chirp(t, f0, T, f1) + 0.1;   % DC offset added

% Define high-order systems
s = tf('s');
G1 = 20 * (s-3)^2 / ((s+1)^2 * (s+0.8)^4);
G2 = 120 / (s+1.6)^6 * exp(-s);   % exact delay: inputdelay

% Simulate output for chirp input
y1_chirp = lsim(G1, u_chirp, t);
y2_chirp = lsim(G2, u_chirp, t);

figure; plot(t, u_chirp); xlabel('Time (s)'); ylabel('Amplitude');
title('Chirp Input Signal'); grid on; saveas(gcf, 'chirp_input.png');

%% Function to compute moments of a signal
% M(k) = sum(t.^k .* signal * Ts)
mom = @(sig, k) sum(t(:).^k .* sig(:) * Ts);

%% Function to estimate model parameters from moments
function [G_fopdt, G_sopdt] = identify_from_moments(Mu, My, k_max)
    % Compute impulse response moments eta_k
    eta = zeros(k_max+1,1);
    eta(1) = My(1) / Mu(1);   % k=0
    if k_max >= 1
        eta(2) = (My(2) - eta(1)*Mu(2)) / Mu(1);
    end
    if k_max >= 2
        eta(3) = (My(3) - eta(1)*Mu(3) - 2*eta(2)*Mu(2)) / Mu(1);
    end
    K = eta(1);
    % FOPDT model: K*exp(-L*s)/(T*s+1)
    % eta0 = K, eta1 = K(T+L), eta2 = K(L^2 + 2*T*L + 2*T^2)
    % Solve for T, L (if possible)
    if K > 1e-6
        TL = eta(2)/K;
        T_sq = eta(3)/K - TL^2;   % T^2
        if T_sq > 0
            T_fopdt = sqrt(T_sq);
            L_fopdt = TL - T_fopdt;
            if L_fopdt < 0, L_fopdt = 0; T_fopdt = TL; end
        else
            T_fopdt = TL;
            L_fopdt = 0;
        end
    else
        T_fopdt = 1; L_fopdt = 0;
    end
    G_fopdt = tf(K, [T_fopdt 1], 'InputDelay', L_fopdt);
    
    % Second-order three-parameter model: K/((T1*s+1)*(T2*s+1))
    % eta0 = K, eta1 = K*(T1+T2), eta2 = 2*K*(T1^2+T2^2+T1*T2)
    if K > 1e-6
        S = eta(2)/K;
        P = S^2 - eta(3)/(2*K);
        disc = S^2 - 4*P;
        if disc >= 0
            T1 = (S + sqrt(disc))/2;
            T2 = (S - sqrt(disc))/2;
            if T1 < 0, T1 = 0; T2 = S; end
            if T2 < 0, T2 = 0; T1 = S; end
        else
            % complex roots -> use a double pole
            T1 = S/2; T2 = S/2;
        end
    else
        T1 = 1; T2 = 1;
    end
    G_sopdt = tf(K, conv([T1 1], [T2 1]));
end

% Compute moments (up to k=2)
Mu_chirp = [mom(u_chirp,0); mom(u_chirp,1); mom(u_chirp,2)];
My1_chirp = [mom(y1_chirp,0); mom(y1_chirp,1); mom(y1_chirp,2)];
My2_chirp = [mom(y2_chirp,0); mom(y2_chirp,1); mom(y2_chirp,2)];

% Identify models for G1
[G1_fopdt_chirp, G1_sopdt_chirp] = identify_from_moments(Mu_chirp, My1_chirp, 2);
fprintf('G1 - Chirp: FOPDT -> K=%.4f, T=%.4f, L=%.4f\n', dcgain(G1_fopdt_chirp), ...
    G1_fopdt_chirp.denominator{1}(2), G1_fopdt_chirp.InputDelay);
fprintf('G1 - Chirp: SOPDT -> K=%.4f, T1=%.4f, T2=%.4f\n', dcgain(G1_sopdt_chirp), ...
    G1_sopdt_chirp.denominator{1}(2), G1_sopdt_chirp.denominator{1}(3));

% Identify models for G2
[G2_fopdt_chirp, G2_sopdt_chirp] = identify_from_moments(Mu_chirp, My2_chirp, 2);
fprintf('G2 - Chirp: FOPDT -> K=%.4f, T=%.4f, L=%.4f\n', dcgain(G2_fopdt_chirp), ...
    G2_fopdt_chirp.denominator{1}(2), G2_fopdt_chirp.InputDelay);
fprintf('G2 - Chirp: SOPDT -> K=%.4f, T1=%.4f, T2=%.4f\n', dcgain(G2_sopdt_chirp), ...
    G2_sopdt_chirp.denominator{1}(2), G2_sopdt_chirp.denominator{1}(3));

% Plot step responses for chirp-based models
figure('Name','G1 - Chirp Identification');
step(G1, 'b-', G1_fopdt_chirp, 'g--', G1_sopdt_chirp, 'r-.', 0:0.1:20);
grid on; legend('True G1', 'FOPDT', 'SOPDT', 'Location','southeast');
title('Step Response Comparison for G1 (Chirp Input)');
saveas(gcf, 'G1_chirp.png');

figure('Name','G2 - Chirp Identification');
step(G2, 'b-', G2_fopdt_chirp, 'g--', G2_sopdt_chirp, 'r-.', 0:0.1:40);
grid on; legend('True G2', 'FOPDT', 'SOPDT', 'Location','southeast');
title('Step Response Comparison for G2 (Chirp Input)');
saveas(gcf, 'G2_chirp.png');

%% Part 2: Arbitrary rich input signal (random with non-zero mean)
rng(0);   % for reproducibility
u_rand = 0.5*randn(size(t)) + 0.5;   % random normal + DC

% Simulate outputs
y1_rand = lsim(G1, u_rand, t);
y2_rand = lsim(G2, u_rand, t);

% Moments
Mu_rand = [mom(u_rand,0); mom(u_rand,1); mom(u_rand,2)];
My1_rand = [mom(y1_rand,0); mom(y1_rand,1); mom(y1_rand,2)];
My2_rand = [mom(y2_rand,0); mom(y2_rand,1); mom(y2_rand,2)];

% Identify models
[G1_fopdt_rand, G1_sopdt_rand] = identify_from_moments(Mu_rand, My1_rand, 2);
fprintf('G1 - Rand: FOPDT -> K=%.4f, T=%.4f, L=%.4f\n', dcgain(G1_fopdt_rand), ...
    G1_fopdt_rand.denominator{1}(2), G1_fopdt_rand.InputDelay);
fprintf('G1 - Rand: SOPDT -> K=%.4f, T1=%.4f, T2=%.4f\n', dcgain(G1_sopdt_rand), ...
    G1_sopdt_rand.denominator{1}(2), G1_sopdt_rand.denominator{1}(3));

[G2_fopdt_rand, G2_sopdt_rand] = identify_from_moments(Mu_rand, My2_rand, 2);
fprintf('G2 - Rand: FOPDT -> K=%.4f, T=%.4f, L=%.4f\n', dcgain(G2_fopdt_rand), ...
    G2_fopdt_rand.denominator{1}(2), G2_fopdt_rand.InputDelay);
fprintf('G2 - Rand: SOPDT -> K=%.4f, T1=%.4f, T2=%.4f\n', dcgain(G2_sopdt_rand), ...
    G2_sopdt_rand.denominator{1}(2), G2_sopdt_rand.denominator{1}(3));

% Plot comparisons for random input
figure('Name','G1 - Random Input Identification');
step(G1, 'b-', G1_fopdt_rand, 'g--', G1_sopdt_rand, 'r-.', 0:0.1:20);
grid on; legend('True G1', 'FOPDT', 'SOPDT', 'Location','southeast');
title('Step Response Comparison for G1 (Random Input)');
saveas(gcf, 'G1_rand.png');

figure('Name','G2 - Random Input Identification');
step(G2, 'b-', G2_fopdt_rand, 'g--', G2_sopdt_rand, 'r-.', 0:0.1:40);
grid on; legend('True G2', 'FOPDT', 'SOPDT', 'Location','southeast');
title('Step Response Comparison for G2 (Random Input)');
saveas(gcf, 'G2_rand.png');

disp('All identifications completed.');