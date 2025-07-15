
%% scaled parameterization: [kappa_bar, gamma, xi^2, mu1, ..., muN]
% where kappa_bar = kappa - beta, gamma = beta/kappa_bar, mu_i = I_i/kappa_bar
classdef ScaledParam < BaseParam
    methods
        function obj = ScaledParam(N, tau)
            obj.N = N;
            obj.tau = tau;
            obj.name = 'scaled (κ̄, γ, ξ², μ)';
            obj.lower_bounds = [1e-9, -10, 1e-9, -10 * ones(1, N)];
            obj.upper_bounds = [100, 10, 100, 10 * ones(1, N)];
        end
        
        function params = compute_initial_guess(obj, data)
            % Get standard parameterization first
            std_param = StandardParam(obj.N, obj.tau);
            params_std = std_param.compute_initial_guess(data);
            
            % Transform to scaled parameterization
            kappa = params_std(1);
            beta = params_std(2);
            xi_squared = params_std(3);
            I = params_std(4:end);
            
            kappa_bar = kappa - beta;
            if kappa_bar < 1e-6
                kappa_bar = 1e-6;  % Ensure positive
            end
            gamma = beta / kappa_bar;
            mu = I / kappa_bar;
            
            params = [kappa_bar, gamma, xi_squared, mu];
        end
        
        function [f, g, H] = compute_likelihood(obj, params, data)
            % Chain rule: ScaledParam -> StandardParam -> Internal
            
            % First transform to standard parameters
            [kappa, beta, xi_squared, I] = obj.to_standard(params);
            params_standard = [kappa, beta, xi_squared, I];
            
            % Jacobian from ScaledParam to Standard
            kappa_bar = params(1);
            gamma = params(2);
            mu = params(4:end);
            
            J_scaled_to_std = zeros(length(params));
            % d(kappa)/d(kappa_bar) = 1 + gamma
            J_scaled_to_std(1,1) = 1 + gamma;
            % d(kappa)/d(gamma) = kappa_bar
            J_scaled_to_std(1,2) = kappa_bar;
            % d(beta)/d(kappa_bar) = gamma
            J_scaled_to_std(2,1) = gamma;
            % d(beta)/d(gamma) = kappa_bar
            J_scaled_to_std(2,2) = kappa_bar;
            % xi_squared unchanged
            J_scaled_to_std(3,3) = 1;
            % d(I_i)/d(kappa_bar) = mu_i
            for i = 1:obj.N
                J_scaled_to_std(3+i,1) = mu(i);
                J_scaled_to_std(3+i,3+i) = kappa_bar;
            end
            
            % Now use StandardParam to go to internal
            std_param = StandardParam(obj.N, obj.tau);
            
            if nargout == 1
                f = std_param.compute_likelihood(params_standard, data);
            elseif nargout == 2
                [f, g_standard] = std_param.compute_likelihood(params_standard, data);
                g = (J_scaled_to_std' * g_standard')';
            else
                [f, g_standard, H_standard] = std_param.compute_likelihood(params_standard, data);
                g = (J_scaled_to_std' * g_standard')';
                H = J_scaled_to_std' * H_standard * J_scaled_to_std;
                
                % Add second-order terms (these come from nonlinear transformations)
                % d²(kappa)/d(kappa_bar)d(gamma) = 1
                H(1,2) = H(1,2) + g_standard(1);
                H(2,1) = H(2,1) + g_standard(1);
                
                % d²(beta)/d(kappa_bar)d(gamma) = 1
                H(1,2) = H(1,2) + g_standard(2);
                H(2,1) = H(2,1) + g_standard(2);
                
                % d²(I_i)/d(kappa_bar)d(mu_i) = 1
                for i = 1:obj.N
                    H(1,3+i) = H(1,3+i) + g_standard(3+i);
                    H(3+i,1) = H(3+i,1) + g_standard(3+i);
                end
            end
        end
        
        function [kappa, beta, xi_squared, I] = to_standard(obj, params)
            % Helper method to transform to standard parameterization
            kappa_bar = params(1);
            gamma = params(2);
            xi_squared = params(3);
            mu = params(4:end);
            
            beta = kappa_bar * gamma;
            kappa = kappa_bar + beta;
            I = kappa_bar * mu;
        end
        
        function results = package_results(obj, params, stderr, log_lik, converged)
            results = struct();
            results.kappa_bar = params(1);
            results.gamma = params(2);
            results.xi_squared = params(3);
            results.xi = sqrt(params(3));
            results.mu = params(4:end);
            results.stderr_kappa_bar = stderr(1);
            results.stderr_gamma = stderr(2);
            results.stderr_xi_squared = stderr(3);
            results.stderr_xi = stderr(3) / (2 * results.xi);
            results.stderr_mu = stderr(4:end);
            results.log_likelihood = log_lik;
            results.converged = converged;
            
            % Also compute derived quantities for reference
            results.beta = results.kappa_bar * results.gamma;
            results.kappa = results.kappa_bar + results.beta;
            results.I = results.kappa_bar * results.mu;
        end
        
        function display_results(obj, results)
            fprintf('  κ̄ (κ-β):         %.4f ± %.4f\n', results.kappa_bar, results.stderr_kappa_bar);
            fprintf('  γ (β/κ̄):         %.4f ± %.4f\n', results.gamma, results.stderr_gamma);
            fprintf('  Noise (ξ):       %.4f ± %.4f\n', results.xi, results.stderr_xi);
            for i = 1:obj.N
                fprintf('  μ_%d:             %.4f ± %.4f (I_%d/κ̄)\n', ...
                    i, results.mu(i), results.stderr_mu(i), i);
            end
            fprintf('\n  Derived: κ = %.4f, β = %.4f\n', results.kappa, results.beta);
        end
    end
end