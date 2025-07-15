# LCAmle - Leaky Competing Accumulator Maximum Likelihood Estimation

A MATLAB toolbox for maximum likelihood estimation of the Leaky Competing Accumulator (LCA) model parameters.

## Project Structure

```
LCAmle/
├── src/
│   ├── params/          # Parameter classes
│   │   ├── BaseParam.m      # Abstract base class
│   │   ├── StandardParam.m  # Standard LCA model ($\kappa$, $\beta$, $\xi^2$, $I$) parameterization
│   │   ├── LogParam.m       # Modified LCA model parameterization
│   │   └── ScaledParam.m     # Alternative scaled type ($\bar{\kappa}=\kappa-\beta$, $\gamma = N\beta/\bar{\kappa}$) parameterization
│   ├── simulation/      # Data simulation
│   │   ├── simulate_LCA_data.m           # Simple wrapper
│   │   └── simulate_LCA_data_taylor15.m  # strong order 1.5 Taylor scheme implementation
│   └── estimation/      # Parameter estimation
│       ├── estimate_LCA.m   # Main estimation interface
│       └── nlog_LCA.m       # Core likelihood computation
├── examples/
│   ├── example_LCA.m        # Usage examples
│   └── my_LCA_results.mat   # Example output
├── LICENSE
└── README.md
```

## Mathematical Model

The LCA model for N alternatives follows the stochastic differential equation:

$$
dx_i = [I_i - \kappa x_i - \beta \sum_{j\neq i}x_j]dt + \xi dW_i
$$

where:
- **$x_i(t)$**: Evidence for alternative i at time t
- **$I_i$**: Input/drift rate for alternative i
- **$\kappa$**: Leak/decay parameter ($\kappa > 0$)
- **$\beta$**: Lateral inhibition strength
- **$\xi$**: Noise intensity (diffusion coefficient)
- **$W_i$**: Independent Wiener processes

## Available Parameterizations

### 1. Standard LCA model Parameterization
**Parameters**: [$\kappa$, $\beta$, $\xi^2$, $I_1$, ..., $I_N$]

This is the most intuitive parameterization using the original model parameters.

```matlab
param = StandardParam(N, tau);
results = estimate_LCA(data, tau, param);
```

### 2. Modified LCA model parameterization
**Parameters**: [$\kappa$, $\beta$, $\xi^2$, $\tilde{I}_1$, ..., $\tilde{I}_N$]

where $\tilde{I}_i = I_i + \frac{1}{2}\xi^2$

This parameterization follows Lo and Ip (2021) modification to introduce a log-normality of the LCA model to ensure positive definiteness. The estimated parameter can be easily obtained by making adjustment to the drift. 

```matlab
param = LogParam(N, tau);
results = estimate_LCA(data, tau, param);
```

### 3. Scaled Parameterization
**Parameters**: [$\bar{\kappa}$, $\gamma$, $\xi^2$, $\mu_1$, ..., $\mu_N$]

where:
- $\bar{\kappa} = \kappa - \beta$
- $\gamma = \beta/\bar{\kappa} = \beta/(\kappa-\beta)$
- $\mu_i = \tilde{I}_i/\bar{\kappa}$

This parameterization:
- Separates timescaled $\bar{\kappa}$ from competition structure ($\gamma$)

```matlab
param = ScaledParam(N, tau);
results = estimate_LCA(data, tau, param);
```

## Usage Example

```matlab
% Load or simulate data
T = 10000;     % Time points
N = 3;         % Number of alternatives
tau = 0.01;    % Time step (seconds)
data = simulate_LCA_data(T, N, tau);

% Estimate parameters using standard parameterization
results = estimate_LCA(data, tau);

% Display results
fprintf('kappa = %.3f, beta = %.3f, xi = %.3f\n', results.kappa, results.beta, results.xi);
fprintf('Inputs: [%.3f, %.3f, %.3f]\n', results.I);
```


## Core Functions

### estimate_LCA - Main Estimation Function

The `estimate_LCA` function is the primary interface for parameter estimation:

```matlab
results = estimate_LCA(data, tau, param)
```

**Inputs:**
- `data`: T×N matrix of observations (T time points, N alternatives)
- `tau`: Time step size (scalar)
- `param`: (Optional) Parameter object specifying the parameterization (default: StandardParam)

**Process:**
1. **Initial Guess**: Calls `param.compute_initial_guess(data)` to get starting values
2. **Optimization**: Automatically configures based on what `compute_likelihood` provides:
   ```matlab
   % If compute_likelihood returns [f, g, H] - uses all analytical derivatives
   options = optimoptions('fmincon', ...
       'SpecifyObjectiveGradient', true, ...
       'HessianFcn', 'objective', ...
       'Algorithm', 'trust-region-reflective');
   
   % If compute_likelihood returns [f, g] - uses gradient only
   options = optimoptions('fmincon', ...
       'SpecifyObjectiveGradient', true, ...
       'Algorithm', 'trust-region-reflective');
   
   % If compute_likelihood returns only f - uses finite differences
   options = optimoptions('fmincon', ...
       'SpecifyObjectiveGradient', false, ...
       'Algorithm', 'sqp');  % or 'interior-point'
   ```
   Using analytical derivatives, we would have a much faster speed, more accurate result and better convergence.
3. **Results**: Packages results using `param.package_results()` and displays using `param.display_results()`

**Output:**
- `results`: Structure containing estimated parameters and diagnostics

### simulate_LCA_data_taylor15 - Data Simulation

Simulates LCA process trajectories using the strong order 1.5 Taylor scheme:

```matlab
data = simulate_LCA_data_taylor15(T, N, tau, kappa, beta, xi_squared, I, x0)
```

**Inputs:**
- `T`: Number of time points
- `N`: Number of alternatives
- `tau`: Time step size
- `kappa`: Leak parameter
- `beta`: Lateral inhibition strength
- `xi_squared`: Noise variance (ξ²)
- `I`: Input rates (N-vector)
- `x0`: (Optional) Initial conditions (default: zeros)

**Algorithm:**
The Taylor 1.5 scheme provides higher accuracy than Euler-Maruyama by including additional stochastic integrals:

$$
\begin{aligned}
x_{i}^{t+\Delta t} = & x_i^t + a_i \Delta t + \xi\Delta W_i - \frac{1}{2}\left\{(\kappa-\beta)a_i + \beta \sum_{j=1}^N a_j \right\} \left(\Delta t\right)^2 \\
&-\xi\left\{(\kappa-\beta)\Delta Z_i + \beta \sum_{j=1}^N \Delta Z_j \right\}
\end{aligned}
$$

where

- $a_i = \left\{(\kappa-\beta)x_i + \beta \displaystyle \sum_{j=1}^N x_j \right\}$
- $\Delta W_i = U_{i,1} \sqrt{\Delta t}$
- $\Delta Z_i = \dfrac{1}{2} \left(\Delta t\right)^{3/2}\left(U_{i,1} + \dfrac{1}{\sqrt{3}}U_{i,2}\right)$

$U_{i,1}$ and $U_{i,2}$ are uncorrelated random numbers drawn from a normal distribution with zero mean and unit variance whilst $\Delta Z_i$ is normally distributed with zero mean, variance $E((\Delta Z_i)^2) = \frac{1}{3}\left(\Delta t\right)^3$ and covariance $E(\Delta Z_i\Delta W_i) = \frac{1}{2}\left(\Delta t\right)^2$.

## Initial Parameter Estimation Method

The toolbox uses a regression-based approach to compute initial parameter guesses:

1. **Sum Process Regression**: For the sum $S(t) = \dfrac{1}{N}\sum_{i=1}^N x_i(t)$, $\bar{I} = \dfrac{1}{N}\sum_{i=1}^N I_i$

   $$
   dS \approx \left\{\bar{I} - \left[\kappa + \beta (N-1)\right] S\right\}dt + noise
   $$

   Linear regression gives estimates for $\bar{I}$ and $\kappa+(N-1)\beta$.

2. **Mean Difference Regression**: For differences $D_{i}(t) = x_i(t) - S(t)$

   $$
   dD_{i} \approx [(I_i - \bar{I}) - (\kappa - \beta)D_{i}]dt + noise
   $$

   Provides estimates for input differences $I_i - \bar{I}$ and $\kappa - \beta$.

3. **Noise Estimation**: Sum of variance of regression residuals estimates $\xi^2$.

4. **Parameter Recovery**: Solve linear system to extract $\kappa$, $\beta$, and individual $I_i$ values.

## Creating Custom Parameterizations

### 1. Define Your Parameterization Class

```matlab
classdef MyCustomParam < BaseParam
    methods
        function obj = MyCustomParam(N, tau)
            obj.N = N;
            obj.tau = tau;
            obj.name = 'My Custom Parameterization';
            % Define bounds for your parameters [θ₁, θ₂, ..., θₘ]
            obj.lower_bounds = [0, -10, 0, -100*ones(1,N)];
            obj.upper_bounds = [100, 100, 100, 100*ones(1,N)];
        end
    end
end
```

### 2. Implement Initial Guess (Example: Using Standard Estimates)

```matlab
function params = compute_initial_guess(obj, data)
    % Get standard parameterization estimates
    std_param = StandardParam(obj.N, obj.tau);
    std_guess = std_param.compute_initial_guess(data);
    
    % Transform to ratio parameterization: [κ, ρ=β/κ, ξ², μ=I/κ]
    kappa = std_guess(1);
    beta = std_guess(2);
    xi_squared = std_guess(3);
    I = std_guess(4:end);
    
    rho = beta / kappa;          % β/κ ratio
    mu = I / kappa;              % scaled inputs
    
    params = [kappa, rho, xi_squared, mu];
end
```

### 3. Chain Rule for Parameter Transformations

When creating custom parameterizations, you need to transform between your parameters $\theta$ and the internal parameters $\psi$ that `nlog_LCA` expects. The chain rule ensures correct gradient and Hessian computations.

#### Mathematical Foundation

For a transformation $\psi = T(\theta)$, where:
- $\theta$ = your custom parameters
- $\psi$ = internal parameters expected by `nlog_LCA`
- $f(\theta) = \ell(\psi(\theta))$ = negative log-likelihood

The chain rule gives:
- **Gradient**: $\nabla_\theta f = J^T \nabla_\psi \ell$
- **Hessian**: $H_\theta = J^T H_\psi J + \sum_{k=1}^{n_\psi} \frac{\partial f}{\partial \psi_k} \frac{\partial^2 \psi_k}{\partial \theta \partial \theta^T}$

where $J = \frac{\partial \psi}{\partial \theta}$ is the Jacobian matrix. The second term in the Hessian only appears for nonlinear transformations.

#### Example: Ratio Parameterization

Consider the parameterization $\theta = [\kappa, \rho, \xi^2, \mu_1, ..., \mu_N]$ where:
- $\rho = \beta/\kappa$ (inhibition-to-leak ratio)
- $\mu_i = I_i/\kappa$ (scaled inputs)

Since `nlog_LCA` internally expects $[\kappa-\beta, \kappa+(N-1)\beta, \xi^2, I_1, ..., I_N]$, we need two transformations:
1. From your parameters $\theta$ to standard parameters $[\kappa, \beta, \xi^2, I]$
2. From standard to internal (handled by StandardParam)

**Step 1: Transform to standard parameters**
$$\begin{align}
\kappa &= \theta_1 \\
\beta &= \kappa \cdot \rho = \theta_1 \cdot \theta_2 \\
\xi^2 &= \theta_3 \\
I_i &= \kappa \cdot \mu_i = \theta_1 \cdot \theta_{3+i} \quad \text{for } i = 1, ..., N
\end{align}$$

**Step 2: Compute the Jacobian**
$$J = \frac{\partial[\kappa, \beta, \xi^2, I]}{\partial\theta} = \begin{pmatrix}
1 & 0 & 0 & 0 & \cdots & 0 \\
\rho & \kappa & 0 & 0 & \cdots & 0 \\
0 & 0 & 1 & 0 & \cdots & 0 \\
\mu_1 & 0 & 0 & \kappa & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
\mu_N & 0 & 0 & 0 & \cdots & \kappa
\end{pmatrix}$$

**Step 3: Handle second-order derivatives for the Hessian**
For nonlinear transformations, the Hessian requires additional terms:
$$H_\theta = J^T H_\psi J + \sum_{k=1}^{n_\psi} g_{\psi,k} \frac{\partial^2 \psi_k}{\partial \theta \partial \theta^T}$$

where $g_{\psi,k}$ is the k-th component of the gradient $\nabla_\psi f$ (i.e., `g_standard(k)` in the code).

The non-zero second derivatives are:
- $\frac{\partial^2 \beta}{\partial \kappa \partial \rho} = 1$
- $\frac{\partial^2 I_i}{\partial \kappa \partial \mu_i} = 1$

All other second derivatives are zero.

### 4. Implement Likelihood with Chain Rule

```matlab
function [f, g, H] = compute_likelihood(obj, params, data)
    % Example: θ = [κ, ρ, ξ², μ₁, ..., μₙ]
    
    % Step 1: Transform to standard parameters [κ, β, ξ², I]
    kappa = params(1);
    rho = params(2);    % β/κ ratio
    xi_squared = params(3);
    mu = params(4:end); % I/κ ratios
    
    beta = kappa * rho;
    I = kappa * mu;
    params_standard = [kappa, beta, xi_squared, I];
    
    % Step 2: Build Jacobian for θ → standard transformation
    J = zeros(length(params));
    J(1,1) = 1;                   % ∂κ/∂κ
    J(2,1) = rho;                 % ∂β/∂κ
    J(2,2) = kappa;               % ∂β/∂ρ
    J(3,3) = 1;                   % ξ² unchanged
    for i = 1:obj.N
        J(3+i,1) = mu(i);         % ∂Iᵢ/∂κ
        J(3+i,3+i) = kappa;       % ∂Iᵢ/∂μᵢ
    end
    
    % Step 3: Use StandardParam (handles transformation to internal)
    std_param = StandardParam(obj.N, obj.tau);
    
    % Step 4: Apply chain rule
    if nargout == 1
        f = std_param.compute_likelihood(params_standard, data);
    elseif nargout == 2
        [f, g_standard] = std_param.compute_likelihood(params_standard, data);
        g = (J' * g_standard')';
    else
        [f, g_standard, H_standard] = std_param.compute_likelihood(params_standard, data);
        g = (J' * g_standard')';
        H = J' * H_standard * J;
        
        % Add second-order terms (only cross terms are non-zero)
        % ∂²β/∂κ∂ρ = 1, so we add g_β * 1
        H(1,2) = H(1,2) + g_standard(2) * 1;
        H(2,1) = H(2,1) + g_standard(2) * 1;
        
        % ∂²Iᵢ/∂κ∂μᵢ = 1, so we add g_Iᵢ * 1  
        for i = 1:obj.N
            H(1,3+i) = H(1,3+i) + g_standard(3+i) * 1;
            H(3+i,1) = H(3+i,1) + g_standard(3+i) * 1;
        end
    end
end
```

### 5. Package Results

```matlab
function results = package_results(obj, params, se, fval)
    results.kappa = params(1);
    results.rho = params(2);    % β/κ ratio
    results.xi_squared = params(3);
    results.mu = params(4:end); % I/κ ratios
    
    % Also include derived standard parameters
    results.beta = results.kappa * results.rho;
    results.I = results.kappa * results.mu;
    results.xi = sqrt(results.xi_squared);
    
    results.log_likelihood = -fval;
end
```

### 6. Display Results

```matlab
function display_results(obj, results)
    fprintf('\n%s Results:\n', obj.name);
    fprintf('κ = %.3f, ρ = %.3f (β/κ ratio)\n', results.kappa, results.rho);
    fprintf('Derived: β = %.3f, ξ = %.3f\n', results.beta, results.xi);
    fprintf('Scaled inputs μ: ');
    fprintf('%.3f ', results.mu);
    fprintf('\n');
end
```

### 7. Complete Implementation Example

```matlab
classdef RatioParam < BaseParam
    methods
        function obj = RatioParam(N, tau)
            obj.N = N;
            obj.tau = tau;
            obj.name = 'Ratio Parameterization (κ, ρ=β/κ, ξ², μ=I/κ)';
            % κ > 0, ρ can be any value, ξ² > 0, μ can be any values
            obj.lower_bounds = [1e-9, -10, 1e-9, -100*ones(1,N)];
            obj.upper_bounds = [100, 10, 100, 100*ones(1,N)];
        end
        
        function params = compute_initial_guess(obj, data)
            % Get standard parameterization estimates
            std_param = StandardParam(obj.N, obj.tau);
            std_guess = std_param.compute_initial_guess(data);
            
            % Transform to ratio parameterization: [κ, ρ=β/κ, ξ², μ=I/κ]
            kappa = std_guess(1);
            beta = std_guess(2);
            xi_squared = std_guess(3);
            I = std_guess(4:end);
            
            rho = beta / kappa;          % β/κ ratio
            mu = I / kappa;              % scaled inputs
            
            params = [kappa, rho, xi_squared, mu];
        end
        
        function [f, g, H] = compute_likelihood(obj, params, data)
            % Step 1: Transform to standard parameters [κ, β, ξ², I]
            kappa = params(1);
            rho = params(2);    % β/κ ratio
            xi_squared = params(3);
            mu = params(4:end); % I/κ ratios
            
            beta = kappa * rho;
            I = kappa * mu;
            params_standard = [kappa, beta, xi_squared, I];
            
            % Step 2: Build Jacobian for θ → standard transformation
            J = zeros(length(params));
            J(1,1) = 1;                   % ∂κ/∂κ
            J(2,1) = rho;                 % ∂β/∂κ
            J(2,2) = kappa;               % ∂β/∂ρ
            J(3,3) = 1;                   % ξ² unchanged
            for i = 1:obj.N
                J(3+i,1) = mu(i);         % ∂Iᵢ/∂κ
                J(3+i,3+i) = kappa;       % ∂Iᵢ/∂μᵢ
            end
            
            % Step 3: Use StandardParam (handles transformation to internal)
            std_param = StandardParam(obj.N, obj.tau);
            
            % Step 4: Apply chain rule
            if nargout == 1
                f = std_param.compute_likelihood(params_standard, data);
            elseif nargout == 2
                [f, g_standard] = std_param.compute_likelihood(params_standard, data);
                g = (J' * g_standard')';
            else
                [f, g_standard, H_standard] = std_param.compute_likelihood(params_standard, data);
                g = (J' * g_standard')';
                H = J' * H_standard * J;
                
                % Add second-order terms (only cross terms are non-zero)
                % ∂²β/∂κ∂ρ = 1, so we add g_β * 1
                H(1,2) = H(1,2) + g_standard(2) * 1;
                H(2,1) = H(2,1) + g_standard(2) * 1;
                
                % ∂²Iᵢ/∂κ∂μᵢ = 1, so we add g_Iᵢ * 1  
                for i = 1:obj.N
                    H(1,3+i) = H(1,3+i) + g_standard(3+i) * 1;
                    H(3+i,1) = H(3+i,1) + g_standard(3+i) * 1;
                end
            end
        end
        
        function results = package_results(obj, params, se, fval)
            results.kappa = params(1);
            results.rho = params(2);    % β/κ ratio
            results.xi_squared = params(3);
            results.mu = params(4:end); % I/κ ratios
            
            % Also include derived standard parameters
            results.beta = results.kappa * results.rho;
            results.I = results.kappa * results.mu;
            results.xi = sqrt(results.xi_squared);
            
            results.log_likelihood = -fval;
        end
        
        function display_results(obj, results)
            fprintf('\n%s Results:\n', obj.name);
            fprintf('κ = %.3f, ρ = %.3f (β/κ ratio)\n', results.kappa, results.rho);
            fprintf('Derived: β = %.3f, ξ = %.3f\n', results.beta, results.xi);
            fprintf('Scaled inputs μ: ');
            fprintf('%.3f ', results.mu);
            fprintf('\n');
        end
    end
end
```

## Bayesian Analysis

To use this toolbox with Bayesian methods, the log-likelihood, gradient, and Hessian function might be useful:

```matlab
% Get log-likelihood and derivatives
param_obj = StandardParam(N, tau);
[neg_log_lik, neg_grad, neg_hess] = param_obj.compute_likelihood(params, data);

% For Bayesian inference
log_likelihood = -neg_log_lik;
gradient = -neg_grad;
hessian = -neg_hess;
```

**Important**: The chain rule is already handled within each parameterization's `compute_likelihood` method. For custom parameterizations, ensure you properly implement the Jacobian transformation as shown in the examples above.

## Requirements

- MATLAB R2020b or later for `pagemtimes`
- Optimization Toolbox for `fmincon`

## License

See LICENSE file for details.