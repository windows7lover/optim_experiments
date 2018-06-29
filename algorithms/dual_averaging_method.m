function [xplus, param] = dual_averaging_method(fun,x,k,param)

% Here, we  regularizer = 0.5||x-x0||

gx = fun.fp;

% Retrive old parameters
beta_hat = param.beta_hat;
sk = param.sk;
rho = param.rho; % = sqrt(2D)
% sigma = param.sigma; % = strong convexity regularizer, usually = 1
sigma = 1; % Here, we  regularizer = 0.5||x-x0||
% x0 = param.x0;
x0 = param.x0;

% Update parameters
gxx = gx(x);
sk = sk + gxx/norm(gxx);
beta_hat = beta_hat + 1/beta_hat;
beta = beta_hat/(rho*sqrt(sigma));

% Update point
xplus = x0-sk/beta;

% Update structure param
param.beta_hat = beta_hat;
param.sk = sk;