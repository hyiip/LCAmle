
%% Adjusted input parameterization: [kappa, beta, xi^2, I_tilde1, ..., I_tildeN]
% where I_tilde_i = I_i - 0.5*xi^2
% Useful when I_tilde has a specific interpretation or for numerical reasons
classdef LogParam < BaseParam
    methods
        function obj = LogParam(N, tau)
            obj.N = N;
            obj.tau = tau;
            obj.name = 'Adjusted Input (κ, β, ξ², Ĩ)';
            obj.lower_bounds = [1e-9, -10, 1e-9, -100 * ones(1, N)];
            obj.upper_bounds = [100, 100, 100, 100 * ones(1, N)];
        end
        
        function params = compute_initial_guess(obj, data)
            % Get standard parameterization first
            std_param = StandardParam(obj.N, obj.tau);
            params_std = std_param.compute_initial_guess(data);
            
            % Transform to lognormal parameterization
            kappa = params_std(1);
            beta = params_std(2);
            xi_squared = params_std(3);
            I = params_std(4:end);
            
            I_tilde = I - 0.5 * xi_squared;
            params = [kappa, beta, xi_squared, I_tilde];
        end
        
        function [f, g, H] = compute_likelihood(obj, params, data)
            % Chain rule: LogParam -> StandardParam -> Internal
            
            % First transform to standard parameters
            [kappa, beta, xi_squared, I] = obj.to_standard(params);
            params_standard = [kappa, beta, xi_squared, I];
            
            % Jacobian from LogParam to Standard
            J_log_to_std = eye(length(params));
            J_log_to_std(4:end, 3) = 0.5;  % dI_i/dxi_squared = 0.5
            
            % Now use StandardParam to go to internal
            std_param = StandardParam(obj.N, obj.tau);
            
            if nargout == 1
                f = std_param.compute_likelihood(params_standard, data);
            elseif nargout == 2
                [f, g_standard] = std_param.compute_likelihood(params_standard, data);
                g = (J_log_to_std' * g_standard')';
            else
                [f, g_standard, H_standard] = std_param.compute_likelihood(params_standard, data);
                g = (J_log_to_std' * g_standard')';
                H = J_log_to_std' * H_standard * J_log_to_std;
            end
        end
        
        function [kappa, beta, xi_squared, I] = to_standard(obj, params)
            % Helper method to transform to standard parameterization
            kappa = params(1);
            beta = params(2);
            xi_squared = params(3);
            I_tilde = params(4:end);
            I = I_tilde + 0.5 * xi_squared;
        end
        
        function results = package_results(obj, params, stderr, log_lik, converged)
            results = struct();
            results.kappa = params(1);
            results.beta = params(2);
            results.xi_squared = params(3);
            results.xi = sqrt(params(3));
            results.I_tilde = params(4:end);
            results.I = results.I_tilde + 0.5 * results.xi_squared;
            results.stderr_kappa = stderr(1);
            results.stderr_beta = stderr(2);
            results.stderr_xi_squared = stderr(3);
            results.stderr_xi = stderr(3) / (2 * results.xi);
            results.stderr_I_tilde = stderr(4:end);
            results.log_likelihood = log_lik;
            results.converged = converged;
        end
        
        function display_results(obj, results)
            fprintf('  Leakage (κ):     %.4f ± %.4f\n', results.kappa, results.stderr_kappa);
            fprintf('  Inhibition (β):  %.4f ± %.4f\n', results.beta, results.stderr_beta);
            fprintf('  Noise (ξ):       %.4f ± %.4f\n', results.xi, results.stderr_xi);
            for i = 1:obj.N
                fprintf('  Ĩ_%d:             %.4f ± %.4f (I_%d = %.4f)\n', ...
                    i, results.I_tilde(i), results.stderr_I_tilde(i), i, results.I(i));
            end
        end
    end
end