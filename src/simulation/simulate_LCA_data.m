function data = simulate_LCA_data(T, N, tau)
    % Simple wrapper with default parameters
    true_kappa = 1;
    true_beta = 0.5;
    true_xi = 0.1;
    true_I = [1.0, 0.8, 0.6];
    
    % Initial condition
    x0 = randn(1, N) * 0.1;
    
    % Call Taylor 1.5 implementation
    data = simulate_LCA_data_taylor15(T, N, tau, true_kappa, true_beta, ...
                                      true_xi, true_I, x0);
end
