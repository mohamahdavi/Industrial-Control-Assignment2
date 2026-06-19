    clc; clear; close all;
    
    %% Dataset filenames
    datasets = {'GrowthData1.csv', 'GrowthData2.csv', 'GrowthData3.csv', 'GrowthData4.csv'};
    nFiles = length(datasets);
    
    % Preallocate storage
    data = cell(1, nFiles);
    metrics = cell(nFiles, 4);   % {R2, RMSE} for each model
    fittedModels = cell(nFiles, 4); % store cfit objects (optional)
    
    %% Load CSV files
    for i = 1:nFiles
        T = readtable(datasets{i});
        if width(T) >= 2
            t_vec = T{:,1};
            y_vec = T{:,2};
        else
            error('Dataset %s must have at least two columns.', datasets{i});
        end
        data{i} = struct('t', t_vec, 'y', y_vec);
    end
    
    %% Define model equations as strings for fittype
    % Logistic: K / (1 + ((K-P0)/P0) * exp(-r*t))
    logistic_eq = 'K ./ (1 + ((K - P0) ./ P0) .* exp(-r .* t))';
    % Gompertz: K * exp(-exp(-r*(t - t0)))
    gompertz_eq = 'K .* exp(-exp(-r .* (t - t0)))';
    % Richards: K ./ ((1 + alpha .* exp(-r .* t)).^(1/nu))
    richards_eq = 'K ./ ((1 + alpha .* exp(-r .* t)).^(1/nu))';
    % Bertalanffy: K * (1 - exp(-r * t))^3
    bertalanffy_eq = 'K .* (1 - exp(-r .* t)).^3';
    
    % Create fittype objects
    ft_log = fittype(logistic_eq, 'independent', 't', 'coefficients', {'K', 'P0', 'r'});
    ft_gomp = fittype(gompertz_eq, 'independent', 't', 'coefficients', {'K', 't0', 'r'});
    ft_rich = fittype(richards_eq, 'independent', 't', 'coefficients', {'K', 'nu', 'r', 'alpha'});
    ft_bert = fittype(bertalanffy_eq, 'independent', 't', 'coefficients', {'K', 'r'});
    
    % Set bounds for parameters (all positive, and nu>0, P0>0, etc.)
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'Display', 'off', ...
                         'StartPoint', []);
    
    %% Fit each dataset
    for d = 1:nFiles
        t = data{d}.t;
        y = data{d}.y;
        n = length(y);
        
        % Initial guesses (will be set inside loop)
        K0 = max(y);
        P00 = y(1);          % Logistic: initial value
        r0 = 0.5;
        t0_guess = mean(t);
        nu0 = 1;
        alpha0 = 1;
        
        % --- Logistic ---
        opts_log = options;
        opts_log.StartPoint = [K0, P00, r0];
        opts_log.Lower = [0, 0, 0];
        [fit_log, gof_log] = fit(t, y, ft_log, opts_log);
        y_pred_log = feval(fit_log, t);
        metrics{d,1} = [gof_log.rsquare, gof_log.rmse];
        fittedModels{d,1} = fit_log;
        
        % --- Gompertz ---
        opts_gomp = options;
        opts_gomp.StartPoint = [K0, t0_guess, r0];
        opts_gomp.Lower = [0, 0, 0];   % t0 can be negative, so relax
        % Actually t0 can be negative; we only bound K and r positive
        opts_gomp.Lower = [0, -inf, 0];
        [fit_gomp, gof_gomp] = fit(t, y, ft_gomp, opts_gomp);
        y_pred_gomp = feval(fit_gomp, t);
        metrics{d,2} = [gof_gomp.rsquare, gof_gomp.rmse];
        fittedModels{d,2} = fit_gomp;
        
        % --- Richards ---
        opts_rich = options;
        opts_rich.StartPoint = [K0, nu0, r0, alpha0];
        opts_rich.Lower = [0, 0, 0, 0];   % all positive
        [fit_rich, gof_rich] = fit(t, y, ft_rich, opts_rich);
        y_pred_rich = feval(fit_rich, t);
        metrics{d,3} = [gof_rich.rsquare, gof_rich.rmse];
        fittedModels{d,3} = fit_rich;
        
        % --- Bertalanffy ---
        opts_bert = options;
        opts_bert.StartPoint = [K0, r0];
        opts_bert.Lower = [0, 0];
        [fit_bert, gof_bert] = fit(t, y, ft_bert, opts_bert);
        y_pred_bert = feval(fit_bert, t);
        metrics{d,4} = [gof_bert.rsquare, gof_bert.rmse];
        fittedModels{d,4} = fit_bert;
        
        % Plot results
        figure;
        plot(t, y, 'ko', 'MarkerSize', 6, 'DisplayName', 'Data'); hold on;
        plot(t, y_pred_log, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Logistic');
        plot(t, y_pred_gomp, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Gompertz');
        plot(t, y_pred_rich, 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Richards');
        plot(t, y_pred_bert, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Bertalanffy');
        hold off;
        xlabel('Time'); ylabel('Population');
        title(sprintf('Dataset %d – Growth Model Fits', d));
        legend('Location', 'best');
        grid on;
        saveas(gcf, sprintf('fit_dataset%d.png', d));
    end
    
    %% Print metrics
    fprintf('\n=== Goodness-of-Fit Metrics ===\n');
    fprintf('%-10s %-15s %-10s %-10s\n', 'Dataset', 'Model', 'R-squared', 'RMSE');
    for d = 1:nFiles
        fprintf('Dataset %d:\n', d);
        fprintf('  Logistic:    R^2 = %.4f, RMSE = %.4f\n', metrics{d,1}(1), metrics{d,1}(2));
        fprintf('  Gompertz:    R^2 = %.4f, RMSE = %.4f\n', metrics{d,2}(1), metrics{d,2}(2));
        fprintf('  Richards:    R^2 = %.4f, RMSE = %.4f\n', metrics{d,3}(1), metrics{d,3}(2));
        fprintf('  Bertalanffy: R^2 = %.4f, RMSE = %.4f\n', metrics{d,4}(1), metrics{d,4}(2));
    end
    
    %% Best model per dataset (highest R^2)
    modelNames = {'Logistic','Gompertz','Richards','Bertalanffy'};
    for d = 1:nFiles
        R2_vals = cellfun(@(x) x(1), metrics(d,:));
        [bestR2, bestIdx] = max(R2_vals);
        fprintf('Best model for Dataset %d: %s (R^2 = %.4f)\n', d, modelNames{bestIdx}, bestR2);
    end