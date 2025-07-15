
%% Base parameter handler class
classdef BaseParam < handle
    properties
        N           % Number of alternatives
        name        % Name of parameterization
        lower_bounds
        upper_bounds
        tau         % Time step - needed for likelihood computation
    end
    
    methods (Abstract)
        % Each parameterization must implement these
        params = compute_initial_guess(obj, data)
        [f, g, H] = compute_likelihood(obj, params, data)
        [kappa, beta, xi_squared, I] = to_standard(obj, params)
        results = package_results(obj, params, stderr, log_lik, converged)
        display_results(obj, results)
    end
    
    methods
        function stderr = compute_standard_errors(obj, hessian, params)
            % Compute standard errors from Hessian
            [R, p] = chol(hessian);
            if p == 0
                Rinv = inv(R);
                cov_matrix = Rinv * Rinv';
                stderr = sqrt(diag(cov_matrix));
            else
                warning('Hessian not positive definite');
                stderr = NaN(size(params));
            end
        end
    end
end