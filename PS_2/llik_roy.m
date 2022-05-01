function f=llik_roy(betas)
% Date: 23-04-2022
% Purpose: computes the sample model for the session

global  y D X Z Q T

nQ=size(Q,2);
nX=size(Q,2);
nZ=size(Q,2);

theta1=betas(1:nQ,1);
theta2=betas(nQ+1:nQ+nQ,1);
delta2=betas(nQ+nQ+1,1);
theta3=betas(nQ+nQ+1+1:nQ+nQ+1+nQ,1);
delta3=betas(nQ+nQ+1+nQ+1,1);
sigma_t1=betas(nQ+nQ+1+nQ+1+1,1);
sigma_t2=betas(nQ+nQ+1+nQ+1+1+1,1);
sigma_t3=betas(nQ+nQ+1+nQ+1+1+1+1,1);
sigma_theta1=betas(nQ+nQ+1+nQ+1+1+1+1+1);

colZ=size(Z,2);
colX=size(X,2);

gamma=betas(3*nQ+6+1:3*nQ+6+colZ,1);
alphav=betas(3*nQ+6+colZ+1,1);
beta1=betas(3*nQ+6+colZ+2:3*nQ+6+colZ+1+colX,1);
alpha1=betas(3*nQ+6+colZ+1+colX+1,1);
beta0=betas(3*nQ+6+colZ+1+colX+2:3*nQ+6+colZ+1+colX+1+colX,1);
alpha0=betas(3*nQ+6+colZ+1+colX+1+colX+1,1);
sigma1=betas(3*nQ+6+colZ+1+colX+1+colX+1+1,1);
sigma0=betas(3*nQ+6+colZ+1+colX+1+colX+1+1+1,1);

fun=@(x) (normpdf(y- X*beta1 - x*alpha1,0,exp(sigma1)).*normcdf(Z*gamma + x*alphav ,0,1) ).^D...
    .*(normpdf(y- X*beta0 - x*alpha0,0,exp(sigma0)).*(1-normcdf(Z*gamma + x*alphav ,0,1)) ).^(1-D)...
    .*normpdf(T(:,1)- Q*theta1 - x,0,exp(sigma_t1))...
    .*normpdf(T(:,2)- Q*theta2 - x*delta2,0,exp(sigma_t2))...
    .*normpdf(T(:,3)- Q*theta3 - x*delta3,0,exp(sigma_t3))...
    .*normpdf(x,0,exp(sigma_theta1));

q = integral(fun,-Inf,Inf,'ArrayValued',true);

f=-sum(log(q));