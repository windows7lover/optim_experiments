
%%% Min square

%%

clear all
close all
clc

addpath(genpath('/home/dscieur/Projects/toolboxes/cvx_free'))
addpath(genpath('../algorithms'))
addpath(genpath('../fixed_point'))

%% Create Problem: min ||Ax-b||^2 + lambda ||x||^2

type = 'RandomStrongConvex';

n = 10;
std_noise = 0.1;
var_noise = std_noise^2/n;
std_noise = sqrt(var_noise);
dist_xstar = 0.1;
norm_xstar = 1e3;

L = 1e1;
mu = 1e0;

opt.L = sqrt(L);
opt.mu = sqrt(mu); % could be zero!
lambda = 0; % could be zero too!

rankA = n;
A = rand(rankA,n);
A = (opt.L-opt.mu)*A/norm(A) + opt.mu*eye(n);
b = ones(rankA,1);

param.A = A;
param.xstar = norm_xstar*rand(n,1);
param.b = A*param.xstar;
random_vec = rand(n,1)-0.5;
random_vec = random_vec/norm(random_vec);

param.x0 = param.xstar + dist_xstar*random_vec;

finfo = getFunction('LeastSquare',param);

finfo.fp = @(x) finfo.fp(x) + std_noise*randn(n,1);

nIteMax = 10;


%% Run algo

x_grad = zeros(n,nIteMax);
x_mean = zeros(n,nIteMax);
x_cheby = zeros(n,nIteMax);

x_grad(:,1) = finfo.x0;
x_mean(:,1) = finfo.x0;
x_cheby(:,1) = finfo.x0;

error_x_grad = zeros(1,nIteMax);
error_x_mean = zeros(1,nIteMax);
error_x_cheby = zeros(1,nIteMax);

error_fun = @(x) norm(x-param.xstar)^2;

%%
for i=1:nIteMax
    
    display(i)
    
    if( i >= 2)
        x_grad(:,i) = x_grad(:,i-1) - (1/L) * finfo.fp(x_grad(:,i-1));
        x_mean(:,i) = (x_mean(:,i-1)*i + x_grad(:,i) ) / (i+1);

        warning off
        [p, max_val_poly, max_val_var, max_val_problem] = stoch_cheby_poly(i-1,n*var_noise/(dist_xstar^2),mu/L,100,1000);
        warning on
        x_cheby(:,i) = x_grad(:,1:i)*flipud(p);
%         x_mean(:,i) = (x_mean(:,i-1)*i + x_cheby(:,i) ) / (i+1);
    end
    
    error_x_grad(i) = error_fun(x_grad(:,i));
    error_x_mean(i) = error_fun(x_mean(:,i));
    error_x_cheby(i) = error_fun(x_cheby(:,i));
    
    
end

save(['nite : ' num2str(nIteMax)])


%%

legendcell = {};
itevec = 1:nIteMax;

figure
semilogy(itevec,error_x_grad); legendcell = [legendcell, 'SGD'];
hold on
semilogy(itevec,error_x_mean); legendcell = [legendcell, 'SGD ave'];
semilogy(itevec,error_x_cheby); legendcell = [legendcell, 'SGD cheby'];
legend(legendcell)



