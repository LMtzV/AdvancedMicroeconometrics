function output = loglikelihood_ps2(init)


% init -> initial guess
global D T W X Y Z

colW = size(W, 2);
colX = size(X, 2);
colZ = size(Z, 2);

% Identification auxiliary
aux = [colX, ...
    2 * colX, ...
    2 * colX + colZ, ...
    2 * colX + colZ + colW, ...
    2 * colX + colZ + 2 * colW, ...
    2 * colX + colZ + 3 * colW, ...
    2 * colX + colZ + 4 * colW, ...
    2 * colX + colZ + 4 * colW + 1, ...
    2 * colX + colZ + 4 * colW + 2, ...
    2 * colX + colZ + 4 * colW + 3, ...
    2 * colX + colZ + 4 * colW + 4, ...
    2 * colX + colZ + 4 * colW + 5, ...
    2 * colX + colZ + 4 * colW + 6, ...
    2 * colX + colZ + 4 * colW + 7, ...
    2 * colX + colZ + 4 * colW + 8, ...
    2 * colX + colZ + 4 * colW + 9, ...
    2 * colX + colZ + 4 * colW + 10, ...
    2 * colX + colZ + 4 * colW + 11, ...
    2 * colX + colZ + 4 * colW + 12, ...
    2 * colX + colZ + 4 * colW + 13];

% Parameters
beta_0 = init(1:aux(1));
beta_1 = init(aux(1) + 1:aux(2));
beta_D = init(aux(2) + 1:aux(3));
omega_T1 = init(aux(3) + 1:aux(4));
omega_T2 = init(aux(4) + 1:aux(5));
omega_T3 = init(aux(5) + 1:aux(6));
omega_T4 = init(aux(6) + 1:aux(7));
sigma_0 = init(aux(7) + 1:aux(8));
sigma_1 = init(aux(8) + 1:aux(9));
sigma_T1 = init(aux(9) + 1:aux(10));
sigma_T2 = init(aux(10) + 1:aux(11));
sigma_T3 = init(aux(11) + 1:aux(12));
sigma_T4 = init(aux(12) + 1:aux(13));
alpha_0 = init(aux(13) + 1:aux(14));
alpha_1 = init(aux(14) + 1:aux(15));
alpha_I = init(aux(15) + 1:aux(16));
alpha_T2 = init(aux(16) + 1:aux(17));
alpha_T3 = init(aux(17) + 1:aux(18));
alpha_T4 = init(aux(18) + 1:aux(19));
sigma_theta = init(aux(19) + 1:aux(20));

MLfunc = @(theta) (normpdf(Y - X * beta_0 - theta * alpha_0, 0, exp(sigma_0)) ...
    .* normpdf(Y - X * beta_1 - theta * alpha_1, 0, exp(sigma_1)) ...
    .* normpdf(T(:, 1) - W * omega_T1 - theta, 0, exp(sigma_T1)) ...
    .* normpdf(T(:, 2) - W * omega_T2 - theta * alpha_T2, 0, exp(sigma_T2)) ...
    .* normpdf(T(:, 3) - W * omega_T3 - theta * alpha_T3, 0, exp(sigma_T3)) ...
    .* normpdf(T(:, 4) - W * omega_T4 - theta * alpha_T4, 0, exp(sigma_T4)) ...
    .* (1 - normcdf(Z * beta_D + theta * alpha_I, 0, 1)) .^ (1 - D) ...
    .* normcdf(Z * beta_D + theta * alpha_I, 0, 1) .^ D ...
    .* normpdf(theta, 0, exp(sigma_theta)));


% Can use Gauss- Hermit or Montecarlo 
% However there is a tradeoff between computing time and programming G-H
% integration.
sol = integral(MLfunc, -Inf, Inf, 'ArrayValued', true);

output = -sum(log(sol));

end
