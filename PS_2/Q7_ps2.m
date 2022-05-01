clear 
% Set Seet: LGMV
randn('seed', 180618)

% Set global 
global D T W X Y Z

%check dummy file read.m
read;

[beta_T1, sigma_T1] = OLS_Est(W, T(:, 1));
[beta_T2, sigma_T2] = OLS_Est(W, T(:, 2));
[beta_T3, sigma_T3] = OLS_Est(W, T(:, 3));
[beta_T4, sigma_T4] = OLS_Est(W, T(:, 4));
[beta_D, sigma_D] = OLS_Est(Z, D); % Take sigma_D = 1, to simplify 
[beta_1, sigma_1] = OLS_Est(X, Y);
[beta_0, sigma_0] = OLS_Est(X, Y);

[alpha_T1, alpha_T2, alpha_T3, alpha_T4, alpha_I, alpha_1, alpha_0, sigma_theta] = deal(1);

% Initial guesses for 
guesses = [beta_0; beta_1; beta_D; beta_T1; beta_T2; beta_T3; beta_T4; ...
    sigma_0; sigma_1; sigma_T1; sigma_T2; sigma_T3; sigma_T4; ...
    alpha_0; alpha_1; alpha_I; alpha_T2; alpha_T3; alpha_T4; ...
    sigma_theta];

tic
options = optimoptions(@fminunc, 'Algorithm', 'quasi-newton', 'Display', 'iter',...
    'GradObj', 'off', 'HessUpdate', 'bfgs', 'UseParallel', false,...
    'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIter', 10000, 'MaxFunEvals', 10000);
[est_beta, est_F, exitflag, output, grad, hessian] = fminunc('loglikelihood_ps2', guesses, options);

runtime = toc;

se_Hess=sqrt(diag(inv(hessian)));
[est_beta se_Hess]