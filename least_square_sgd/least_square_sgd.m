
%%% Min square

%%

clear all
close all
clc

addpath(genpath('../algorithms'))
addpath(genpath('../fixed_point'))

%% Create Problem: min ||Ax-b||^2 + lambda ||x||^2

type = 'RandomStrongConvex';

opt.L = 1e1;
opt.mu = 1e-1; % could be zero!
lambda = 0; % could be zero too!
n = 300;
% rankA = round(sqrt(n));
rankA = n;
[Q1,~] = qr(rand(n));
[Q2,~] = qr(rand(n));
S = sqrt(diag([opt.mu + (opt.L-opt.mu) * rand(n-2,1) ; opt.mu ; opt.L]));
A = Q1 * S * Q2; 
b = ones(rankA,1);

param.A = A;
param.xstar = rand(n,1);
param.b = A*param.xstar;
param.x0 = 1000*rand(n,1);

finfo = getFunction('LeastSquare',param);

std_noise = 1000/sqrt(n);

% finfo.fp = @(x) finfo.fp(x) + std_noise*randn(n,1); % additive noise
finfo.fp = @(x) finfo.fp(x) + std_noise*(rand(n,1)-0.5)/sqrt(12); % additive noise
finfo.proxoperator.f = @(x) 0;

x0 = finfo.x0;
xstar = finfo.xstar;
L = finfo.L;

nIteMax = 10000;

%% RMPE param

rmpe_k = 20;

param_rmpe.doAdaptiveLambda = false; % Determine if lambda should change over time
param_rmpe.lambda = 0; % Determine if lambda should change over time
% param_rmpe.lambda = 1e-6; % Determine if lambda should change over time
param_rmpe.lambdamin = 1e-12; % Determine if lambda should change over time
param_rmpe.lambdaSVD = true; % Determine lambda in function of the SVD
param_rmpe.k = rmpe_k; % Number of iteration before extrapolation then restart
param_rmpe.rmpealgo=@ampe; % Extrapolation algorithm
param_rmpe.doLineSearch = false;
param_rmpe.forceDecrease = false;

%% SGD
x = x0;
err_sgd = zeros(nIteMax,1);
for idx = 1:nIteMax
    err_sgd(idx) = finfo.f(x) - finfo.f(xstar);
    x = x-(1/L)*finfo.fp(x);
end


%% aveSGD
x_ave = x0;
err_ave_sgd = zeros(nIteMax,1);
for idx = 1:nIteMax
    err_ave_sgd(idx) = finfo.f(x_ave) - finfo.f(xstar);
    x = x-(1/L)*finfo.fp(x);
    x_ave = x_ave*idx/(idx+1) + x*1/(idx+1);
end


%% Acc. SGD (Flammarion)

trC = std_noise^2*n;

x = x0;
xold = x0;
err_acc_sgd = zeros(nIteMax,1);
for idx = 1:nIteMax
    err_acc_sgd(idx) = finfo.f(x) - finfo.f(xstar);
    
    alpha = min(1/L , norm(x0-xstar)/(2*sqrt(trC)*idx^(3/2)));
    beta = min(idx*alpha,1/L);
    
    xgrad = x*idx*(alpha+beta)/(idx*alpha + beta) - xold*(idx-1)*beta/(idx*alpha + beta);
    
    xplus = x*2*idx/(idx+1) - xold*(idx-1)/(idx+1) - finfo.fp(xgrad) * (idx*alpha+beta)/(idx+1);
    xold = x;
    x = xplus;
end


%% RMPE
nIteRmpe = nIteMax/rmpe_k;
err_rmpe_sgd = zeros(nIteRmpe,1);
idx_vec_rmpe = zeros(nIteRmpe,1);

grad_call = 1;
x = x0;
warning('off')
for idx = 1:nIteRmpe
    
    err_rmpe_sgd(idx) = finfo.f(x) - finfo.f(xstar);
    idx_vec_rmpe(idx) = grad_call;
    
    y_vec = zeros(n,rmpe_k+1);
    y = x;
    y_vec(:,1) = y;
    for idx2 = 1:rmpe_k
        y = y-(1/L)*finfo.fp(y);
        grad_call = grad_call+1;
        y_vec(:,idx2+1) = y;
    end
    xnew = abstract_extrapolation_adaptive_lambda_ls_stepsize(finfo,y_vec,param_rmpe);
    [norm(x-xnew) norm(finfo.fp(x))/L]
    x = xnew;
end
warning('on')

%% Plot

figure
loglog(err_sgd,':b','linewidth',2)
hold on
loglog(err_ave_sgd,'--g','linewidth',2)
loglog(err_acc_sgd,'-.r','linewidth',2)
loglog(idx_vec_rmpe,err_rmpe_sgd,'c','linewidth',2)
legend({'sgd','avesgd','accsgd','rmpesgd'})

xlabel('xlabel','fontsize',16)
ylabel('valf','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
axis tight


