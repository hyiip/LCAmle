function data = simulate_LCA_data_taylor15(T, N, tau, kappa, beta, xi, I, x0)
    % SIMULATE_LCA_DATA - Simulate LCA using strong order 1.5 Taylor scheme
    %
    % This implements the strong order 1.5 Taylor scheme
    % for improved accuracy compared to Euler-Maruyama
    %
    % Inputs:
    %   T - number of time points
    %   N - number of alternatives
    %   tau - time step (dt)
    %   kappa - leakage parameter
    %   beta - inhibition parameter
    %   xi - noise parameter
    %   I - input vector (1 x N)
    %   x0 - initial condition (1 x N)
    %
    % Output:
    %   data - T x N matrix of simulated data
    
    % Initialize
    data = zeros(T, N);
    data(1, :) = x0;
    
    % Pre-compute constants
    sqrt_tau = sqrt(tau);
    sqrt_tau3 = sqrt(tau^3);
    half_tau2 = 0.5 * tau^2;
    
    % Main simulation loop
    for t = 2:T
        x = data(t-1, :);
        sumX = sum(x);
        
        % Generate random numbers for each alternative
        % ran1 for dW, ran2 for higher-order terms
        ran1 = randn(1, N);
        ran2 = randn(1, N);
        
        % Drift term: a = I - (κ-β)x - β*sum(x)
        a = I - (kappa - beta) * x - beta * sumX;
        
        % Brownian increments
        dW = sqrt_tau * ran1;
        dZ = 0.5 * sqrt_tau3 * (ran1 + ran2/sqrt(3));
        
        % Sum of drift for inhibition terms
        sum_a = sum(a);
        sum_dZ = sum(dZ);
        
        % Taylor 1.5 update for each alternative
        for j = 1:N
            % First-order terms
            dx_order1 = a(j) * tau + xi * dW(j);
            
            % Second-order deterministic term (from Taylor expansion of drift)
            dx_order2_det = -(kappa - beta) * a(j) - beta * sum_a;
            dx_order2_det = dx_order2_det * half_tau2;
            
            % Second-order stochastic term (Milstein correction)
            dx_order2_stoch = -xi * ((kappa - beta) * dZ(j) + beta * sum_dZ);
            
            % Total increment
            dx = dx_order1 + dx_order2_det + dx_order2_stoch;
            
            % Update
            data(t, j) = x(j) + dx;
        end
    end
end
