
%% Standard parameterization: [kappa, beta, xi^2, I1, ..., IN]
classdef StandardParam < BaseParam
    methods
        function obj = StandardParam(N, tau)
            obj.N = N;
            obj.tau = tau;
            obj.name = 'Standard (κ, β, ξ², I)';
            obj.lower_bounds = [1e-9, -10, 1e-9, -100 * ones(1, N)];
            obj.upper_bounds = [100, 100, 100, 100 * ones(1, N)];
        end
        
        function params = compute_initial_guess(obj, data)
            % Standard initial guess computation
            [T, N] = size(data);
            tau = obj.tau;  % Use tau from object
            
            % Regression on sum process
            xmean = mean(data, 2);

            xsum = xmean;
            dxsum = diff(xsum);
            xsum = xsum(1:end-1);
            
            regressors_sum = [ones(length(xsum),1)*tau -xsum*tau];
            drift_sum = regressors_sum\dxsum;
            res_sum = regressors_sum * drift_sum - dxsum;
            var_sum = var(res_sum,1)/tau;
            
            % Regression on mean differences
            drift_diff = zeros(2,N);
            var_diff = zeros(1,N);
            for k = 1:N
                xdiffi = data(:,k) - xmean;
                dxdiffi = diff(xdiffi);
                xdiffi = xdiffi(1:end-1);
                regressors_diffi2 = [ones(length(xdiffi),1)*tau -xdiffi*tau];
                drift_diff(:, k) = regressors_diffi2\dxdiffi;
                res_diffi = regressors_diffi2 * drift_diff(:, k) - dxdiffi;
                var_diff(k) = var(res_diffi,1)/tau;
            end
            
            % Estimate parameters
            xi_squared = (sum(var_diff) + var_sum) / N;
            kb_diff_m = mean(drift_diff(2,:));
            kb_sum_m = drift_sum(2);
            
            kb_matrix = [1, -1; 1, N-1];
            kb_vector = [kb_diff_m; kb_sum_m];
            kb = kb_matrix \ kb_vector;
            
            kappa = kb(1);
            beta = kb(2);
            
            % Estimate individual inputs
            I_bar_i = drift_sum(1);
            I_i = drift_diff(1,:) + I_bar_i;
            
            params = [kappa, beta, xi_squared, I_i];
        end
        
        function [f, g, H] = compute_likelihood(obj, params, data)
            % Transform to internal parameterization and compute likelihood
            kappa = params(1);
            beta = params(2);
            
            % Transform to numerically robust internal parameters
            kb_diff = kappa - beta;
            kb_sum = kappa + (obj.N - 1) * beta;
            params_internal = [kb_diff, kb_sum, params(3:end)];
            
            % Jacobian for chain rule
            J = eye(length(params));
            J(1:2, 1:2) = [1, -1; 1, obj.N-1];
            
            % Call core function with tau
            if nargout == 1
                f = nlog_LCA(params_internal, data, obj.tau);
            elseif nargout == 2
                [f, g_internal] = nlog_LCA(params_internal, data, obj.tau);
                g = (J' * g_internal')';
            else
                [f, g_internal, H_internal] = nlog_LCA(params_internal, data, obj.tau);
                g = (J' * g_internal')';
                H = J' * H_internal * J;
            end
        end
        
        function [kappa, beta, xi_squared, I] = to_standard(obj, params)
            % Already in standard form
            kappa = params(1);
            beta = params(2);
            xi_squared = params(3);
            I = params(4:end);
        end
        
        function results = package_results(obj, params, stderr, log_lik, converged)
            results = struct();
            results.kappa = params(1);
            results.beta = params(2);
            results.xi_squared = params(3);
            results.xi = sqrt(params(3));
            results.I = params(4:end);
            results.stderr_kappa = stderr(1);
            results.stderr_beta = stderr(2);
            results.stderr_xi_squared = stderr(3);
            results.stderr_xi = stderr(3) / (2 * results.xi);
            results.stderr_I = stderr(4:end);
            results.log_likelihood = log_lik;
            results.converged = converged;
        end
        
        function display_results(obj, results)
            fprintf('  Leakage (κ):     %.4f ± %.4f\n', results.kappa, results.stderr_kappa);
            fprintf('  Inhibition (β):  %.4f ± %.4f\n', results.beta, results.stderr_beta);
            fprintf('  Noise (ξ):       %.4f ± %.4f\n', results.xi, results.stderr_xi);
            for i = 1:obj.N
                fprintf('  Input I_%d:       %.4f ± %.4f\n', i, results.I(i), results.stderr_I(i));
            end
        end
    end
end