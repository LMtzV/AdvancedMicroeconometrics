% Date: 23-04-2022
% Purpose: Optimize the log-likelihood for our extended Roy model

clear 
randn('seed',12345)
global y D X Z Q T


%%

% load data
data=csvread('sim_toestimate.csv');

ones=ones(size(data,1),1);

% t1 t2 t3 D Y Q X Z
t1=data(:,1);
t2=data(:,2);
t3=data(:,3);
T=[t1 t2 t3];
D=data(:,4);
y=data(:,5);
Q=[ones data(:,6)];
X=[ones data(:,7)];
Z=[ones data(:,8)];

%% Starting values for the betas

% Measuring system
beta_start_T1=inv(Q'*Q)*(Q'*T(:,1));
sigma_start_T1=sqrt(((T(:,1)-Q*beta_start_T1)'*(T(:,1)-Q*beta_start_T1))/(size(T,1)-size(Q,2)));

beta_start_T2=inv(Q'*Q)*(Q'*T(:,2));
sigma_start_T2=sqrt(((T(:,2)-Q*beta_start_T2)'*(T(:,2)-Q*beta_start_T2))/(size(T,1)-size(Q,2)));
delta_start_T2=1;

beta_start_T3=inv(Q'*Q)*(Q'*T(:,3));
sigma_start_T3=sqrt(((T(:,3)-Q*beta_start_T3)'*(T(:,3)-Q*beta_start_T3))/(size(T,1)-size(Q,2)));
delta_start_T3=1;

sigma_theta1_start=1;
sigma_theta2_start=1;


betas_start_1s=[beta_start_T1; beta_start_T2 ; delta_start_T2; beta_start_T3;...
    delta_start_T3; sigma_start_T1; sigma_start_T2; sigma_start_T3 ; sigma_theta1_start];

% Rest of the parameters
beta_start_D=inv(Z'*Z)*(Z'*D);
alphav_start=1;

beta_start_y1=inv(X'*X)*(X'*y);
sigma_start_y1=sqrt(((y-X*beta_start_y1)'*(y-X*beta_start_y1))/(size(y,1)-size(X,2)));
alpha1_start=1;

beta_start_y0=inv(X'*X)*(X'*y);
sigma_start_y0=sqrt(((y-X*beta_start_y0)'*(y-X*beta_start_y0))/(size(y,1)-size(X,2)));
alpha0_start=1;

betas_start_2s=[beta_start_D; alphav_start ;beta_start_y1 ;alpha1_start;...
    beta_start_y0;alpha0_start ; sigma_start_y1 ; sigma_start_y0];


%% Optimization - all at once
betas_start=[betas_start_1s ; betas_start_2s];

of=llik_roy(betas_start)

tic
options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display','iter',...
    'GradObj','off','HessUpdate','bfgs','UseParallel',false,...
    'TolFun',1e-6,'TolX',1e-6,'MaxIter',1e6,'MaxFunEvals',1e6);
[beta_est,fval,exitflag,output,grad,hessian] = fminunc('llik_roy',betas_start,options);
runtime=toc

se_H=sqrt(diag(inv(hessian)));

[beta_est se_H]



