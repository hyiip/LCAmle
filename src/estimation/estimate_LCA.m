
function results = estimate_LCA(data, tau, param_handler)
    % ESTIMATE_LCA - Estimate LCA model parameters
    %
    % Inputs:
    %   data - T x N matrix of observations
    %   tau  - time step size
    %   param_handler - (optional) parameter handler object
    %
    % Outputs:
    %   results - structure with parameter estimates and diagnostics
    
    [T, N] = size(data);
    
    % Default to standard parameterization
    if nargin < 3
        param_handler = StandardParam(N, tau);
    end
    
    fprintf('Estimating LCA model: %d time points, %d alternatives\n', T, N);
    fprintf('Using parameterization: %s\n', param_handler.name);
    
    % Get initial guess
    initial_guess = param_handler.compute_initial_guess(data);
    
    % Optimization settings
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'trust-region-reflective', ...
        'SpecifyObjectiveGradient', true, ...
        'HessianFcn', 'objective', ...
        'OptimalityTolerance', 1e-6);
    
    % Objective function
    objective = @(x) param_handler.compute_likelihood(x, data);
    
    % Run optimization
    fprintf('\nStarting optimization...\n');
    [params_opt, fval, exitflag, output, ~, ~, hessian] = ...
        fmincon(objective, initial_guess, [], [], [], [], ...
                param_handler.lower_bounds, param_handler.upper_bounds, [], options);
    
    % Compute standard errors
    stderr = param_handler.compute_standard_errors(hessian, params_opt);
    
    % Package results (let param handler decide what to store)
    results = param_handler.package_results(params_opt, stderr, -fval, exitflag > 0);
    
    % Display results
    fprintf('\n========== Results ==========\n');
    fprintf('Parameterization: %s\n', param_handler.name);
    fprintf('Converged: %s\n', mat2str(results.converged));
    fprintf('Log-likelihood: %.2f\n', results.log_likelihood);
    fprintf('\nParameter estimates:\n');
    param_handler.display_results(results);
end