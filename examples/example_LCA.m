%% example_LCA.m
% Example of using the LCA model estimation
%
% The LCA (Leaky Competing Accumulator) model describes how people
% accumulate evidence for different choices over time.

clear;
clc;

% Add paths to source directories
addpath('../src/estimation');
addpath('../src/simulation');
addpath('../src/params');

%% Load your data
% Your data should be a matrix with:
% - Rows = time points (e.g., measurements over time)
% - Columns = different choice alternatives

% Example: Load from CSV file
% data = readmatrix('my_experiment_data.csv');

% For this example, we'll use simulated data
T = 10000;  % 500 time points
N = 3;    % 3 choice alternatives
tau = 0.01;  % Time step (10 milliseconds)

% Simulate some data (in practice, use your real data)
rng(42);  % For reproducibility
data = simulate_LCA_data(T, N, tau);

%% Example 1: Standard parameterization (most common)
fprintf('=== Example 1: Standard Analysis ===\n');
results = estimate_LCA(data, tau);

% Access your results
fprintf('\nYour parameter estimates:\n');
fprintf('κ = %.3f, β = %.3f, ξ = %.3f\n', results.kappa, results.beta, results.xi);
fprintf('Inputs: [%.3f, %.3f, %.3f]\n\n', results.I);

%% Example 2: If you want to estimate adjusted inputs Ĩ
fprintf('\n=== Example 2: Adjusted Input Parameterization ===\n');
param_log = LogParam(N, tau);
results_log = estimate_LCA(data, tau, param_log);

% You directly estimate Ĩ = I - 0.5ξ²
fprintf('\nĨ values: [%.3f, %.3f, %.3f]\n', results_log.I_tilde);
fprintf('(Equivalent I values: [%.3f, %.3f, %.3f])\n', results_log.I);

%% Example 3: If you want scaled parameterization
fprintf('\n=== Example 3: scaled Parameterization ===\n');
param_scaled = ScaledParam(N, tau);
results_scaled = estimate_LCA(data, tau, param_scaled);

% You get κ̄ and γ directly
fprintf('\nκ̄ = %.3f, γ = %.3f\n', results_scaled.kappa_bar, results_scaled.gamma);
fprintf('μ values: [%.3f, %.3f, %.3f]\n', results_scaled.mu);
