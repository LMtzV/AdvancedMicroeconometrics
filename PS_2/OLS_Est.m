function [beta, sigma] = OLS_Est(X, Y)
%OLS regression
beta = inv(X' * X) * X' * Y;
eps = Y - X*beta;
sigma = sqrt((eps' * eps) / (size(Y, 1) - size(X, 2)));
end