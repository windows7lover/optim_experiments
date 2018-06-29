% Basic use of RMPE
%
% The purpose of this file is to show how one should plug the RMPE
% algorithm to the gradient method. Moreover, a comparison is done between
% RMPE (regularized acceleration) and AMPE (extrapolation without 
% regularization).
% Since this file generates randomly the dataset, do not hesitate to launch
% it several times to see the difference between the raw extrapolation and
% the regularized extrapolation.
%
% Note also that, in theory, the regularization parameter should decrease
% over time. It means that the line-search on the regularizer is mandatory,
% but for simplicity this line-search is not coded here. Fore more
% information about a proper implementation, see the "advanced folder".


%% Clear everything

clear all
close all
clc

%% Parameters

k_rmpe = 5;

% Regularization parameter: Even if it may works, setting reg_rmpe = 0
% leads sometimes to an unstable and/or divergent process
% Moreover, fixing reg_rmpe is not a good idea since reg_rmpe should
% decrease over time (reg_rmpe = O(perturbation) in the theory ), but it is
% a toy experiment.
reg_rmpe = 1e-12;   % even if it appears small, keep in mind that it is in 
                    % fact the square of the regularization parameter

nIte = k_rmpe*1e2; % Should be a multiple of k_rmpe
lambda_logistic_ratio = 1e-6 ;

% Define some random dataset
nPoints = 1e3;
nFeatures = 1e1;
sigma = rand(nFeatures);
sigma = (sigma'*sigma)^2; % kind of covariance matrix, to make X badly conditionned
X = sigma*(rand(nFeatures,nPoints)-0.5);
y = sign(rand(nPoints,1)-0.5);

lambda_logistic = lambda_logistic_ratio * norm(X'*X)/4;


%% Define Logistic regression   

yXt = (X*spdiags(y,0,nPoints,nPoints))';
    
f = @(w) sum(  log(1+exp(-(yXt*w) )  )) + (lambda_logistic/2)*norm(w)^2;
grad = @(w) -sum( X*spdiags(y./(1+exp(yXt*w)),0,nPoints,nPoints) ,2) + lambda_logistic*w;

L = norm(X*X')/4+lambda_logistic ;
mu = lambda_logistic;


%% Mem alloc
x0 = zeros(nFeatures,1);

tolf_grad = [f(x0),zeros(1,nIte)];
tolf_ampe = [f(x0),zeros(1,nIte/k_rmpe)];
tolf_rmpe = [f(x0),zeros(1,nIte/k_rmpe)];

fstar_approx = -inf;
fstar_approx_handle = @(x) f(x) - (1/(2*mu))*norm(grad(x))^2; % Duality grap for strongly convex functions

%% Regular Gradient method
xnew = x0;
for i=1:nIte
    xold = xnew;
    xnew = xold - grad(xold)*2/(mu+L); % Gradient step for strongly convex objective functions
    tolf_grad(i+1) = f(xnew);
    
    % Update lower bound on fstar
    new_fstar_approx = fstar_approx_handle(xnew);
    if(new_fstar_approx > fstar_approx && new_fstar_approx~=inf)
        fstar_approx = new_fstar_approx ;
    end
end


%% Approximate MPE (AMPE) on gradient method (no regularization)


warning('off','MATLAB:nearlySingularMatrix')

x_temp = zeros(nFeatures,1+k_rmpe);
onesk = ones(k_rmpe,1);
xnew = x0;
for i=1:(nIte/k_rmpe)
    
    xold = xnew;
    
    % Do k steps of gradient method
    x_temp(:,1) = xold;
    for j=1:k_rmpe
        xold = x_temp(:,j);
        x_temp(:,j+1) = xold - grad(xold)*2/(mu+L);
    end
    
    % Extrapolate using AMPE algorithm
    U = diff(x_temp,1,2);
    UU = U'*U;
    UU = UU/norm(UU); % Normalization of U'U
    z = UU\onesk; % We do not use regularization here
    c = z/sum(z);
    x_extrapolated = x_temp(:,2:end)*c; % Extrapolation is a linear combination of prevous iterates
    
    % "Restart" by setting x0 = x_extrapolated
    xnew = x_extrapolated;
    tolf_ampe(i+1) = f(xnew);
    
    
    % Update lower bound on fstar
    new_fstar_approx = fstar_approx_handle(xnew);
    if(new_fstar_approx > fstar_approx && new_fstar_approx~=inf)
        fstar_approx = new_fstar_approx ;
    end
end


warning('on','MATLAB:nearlySingularMatrix')



%% Regularized MPE (RMPE) on gradient method


warning('off','MATLAB:nearlySingularMatrix')

x_temp = zeros(nFeatures,1+k_rmpe);
reg_rmpe_eyek = reg_rmpe*eye(k_rmpe);
onesk = ones(k_rmpe,1);
xnew = x0;
for i=1:(nIte/k_rmpe)
    
    xold = xnew;
    
    % Do k steps of gradient method
    x_temp(:,1) = xold;
    for j=1:k_rmpe
        xold = x_temp(:,j);
        x_temp(:,j+1) = xold - grad(xold)*2/(mu+L);
    end
    
    % Extrapolate using RMPE algorithm
    U = diff(x_temp,1,2);
    UU = U'*U;
    UU = UU/norm(UU); % Normalization of U'U, so that lambda scales better
    z = (UU+reg_rmpe_eyek)\onesk;
    c = z/sum(z);
    x_extrapolated = x_temp(:,2:end)*c; % Extrapolation is a linear combination of prevous iterates
    
    % "Restart" by setting x0 = x_extrapolated
    xnew = x_extrapolated;
    tolf_rmpe(i+1) = f(xnew);
    
    
    % Update lower bound on fstar
    new_fstar_approx = fstar_approx_handle(xnew);
    if(new_fstar_approx > fstar_approx && new_fstar_approx~=inf)
        fstar_approx = new_fstar_approx ;
    end
end


warning('on','MATLAB:nearlySingularMatrix')


%% Update tolerance vectors

tolf_grad_updated = abs(tolf_grad-fstar_approx);
tolf_ampe_updated = abs(tolf_ampe-fstar_approx);
tolf_rmpe_updated = abs(tolf_rmpe-fstar_approx);

%% Plot

figure
semilogy(0:nIte,tolf_grad_updated,'r','linewidth',2)
hold on
semilogy(0:k_rmpe:nIte,tolf_ampe_updated,'g','linewidth',2)
semilogy(0:k_rmpe:nIte,tolf_rmpe_updated,'c','linewidth',2)
xlabel('Gradient oracle calls')
ylabel('Duality gap (estimated)')

legend({'Gradient method','AMPE-k on gradient','RMPE-k on gradient'})
