clc; clear; close all;

%% Load data and inspect variables
load('Ident.mat');
disp('Variables in Ident.mat:');
whos

% The file contains X (input) and Y (output), both 20000x1
u_raw = X;      % input
y_raw = Y;      % output

% The experiment is a step/impulse type. The toolbox requires that the initial
% portion of the data represents the system at equilibrium (steady state) before
% the step occurs. We therefore prepend 50 zero samples to both signals.
% This ensures that the model can capture the pre‑step steady state.
nx_pre = 50;    % number of extra zero samples (greater than max model order)

u = [zeros(nx_pre, 1); u_raw];
y = [zeros(nx_pre, 1); y_raw];

% Sampling and start time
Ts = 0.1;               % sample interval
Tstart = 1;             % original start time (now effectively shifted)

% Create iddata object
data = iddata(y, u, Ts, 'Tstart', Tstart);

%% Preprocess: remove mean
data = detrend(data, 0);   % 0 = remove constant (mean)

%% Estimate transfer function models (1..5 poles, 0..2 zeros)
best_fit = -Inf;
best_model = [];
best_order = [];
fits = [];

fprintf('\n=== Transfer Function Estimation Results ===\n');
fprintf('%-10s %-10s %-10s %-10s %-10s %-10s\n', ...
    'np', 'nz', 'Fit (%)', 'MSE', 'AIC', 'FPE');

for np = 1:5
    for nz = 0:min(2, np-1)   % number of zeros < poles for a proper model
        try
            sys = tfest(data, np, nz);
            [~, fit, ~] = compare(data, sys);
            % Residual analysis for MSE
            resid_data = resid(data, sys);
            e = resid_data.OutputData;
            N = length(e);
            V = (e' * e) / N;            % mean squared error
            n_par = np + nz + 1;        % number of estimated parameters
            AIC_val = log(V) + 2*n_par/N;
            FPE_val = V * (1 + 2*n_par/N);  % approximate FPE
            fprintf('%-10d %-10d %-10.2f %-10.4f %-10.4f %-10.4f\n', ...
                np, nz, fit, V, AIC_val, FPE_val);
            fits = [fits; np, nz, fit, V, AIC_val, FPE_val];
            if fit > best_fit
                best_fit = fit;
                best_model = sys;
                best_order = [np, nz];
            end
        catch
            % tfest may fail for some orders; simply skip
        end
    end
end

%% Display best model
if isempty(best_model)
    error('No model could be estimated. Check the data and prepended zeros.');
end
fprintf('\nBest model: np = %d, nz = %d, Fit = %.2f%%\n', ...
    best_order(1), best_order(2), best_fit);
disp('Transfer function:')
tf(best_model)

%% Plot results
figure;
compare(data, best_model);
title(sprintf('Best fit: %.2f%% (np=%d, nz=%d)', best_fit, best_order(1), best_order(2)));
saveas(gcf, 'ident_fit.png');

figure;
resid(data, best_model);
title('Residual Analysis');
saveas(gcf, 'ident_resid.png');

% Save the best model and performance metrics
save('estimated_model.mat', 'best_model', 'best_order', 'best_fit', 'fits');