%% Example of use of algorithm "rmpe_blackbox" algorithm
% 
% Here, we generate a quadratic minimization problem, then define the
% gradient method. Finally, we use the black box rmpe algorithm to
% accelerate the rate of convergence of gradient method using only its
% iterates and the objective function.

%% Clear everything

clear all
close all
clc
RmpeParam = [];

%% Parameter of RMPE
% These are the default parameters of RMPE.

% RmpeParam.doAdaptiveLambda = true; 
% % Determine if lambda should change over time

% RmpeParam.doLineSearch = true; 
% % Determine if we should perform a line search at the end

% RmpeParam.lambda = 1; 
% RmpeParam.lambdamin = 1e-10;
% % These two values determine the range of the grid search [lambda,lambdamin]

% RmpeParam.forceDecrease = true; 
% % Optionnal, check if the extrapolated value is smaller than the last iterate of gradient method

%% Sequence size

k = 5; 
% This means k oracle calls, so the length of the sequence is k+1.
% => Define the number of points used in the extrapolation algorithm.
% A good value for k is, usually, between 5 and 10


%% Generation of the problem
% here we minimize 0.5*||Ax-b||^2.
% for simplicity, we will use gradient method

dim = 100; % Dim. of the problem
L = 1e3; % norm of A'*A

x0 = zeros(dim,1);
xstar = rand(dim,1);
A = rand(dim,dim);
AA = A'*A;

% Scaling so that the function is L-smooth
normAA = norm(AA);
A = A*sqrt(L/normAA);
AA = AA*(L/normAA);

b = A*xstar;
Ab = A'*b;

fx = @(x) 0.5*norm(A*x-b)^2;
error_fun = @(x) norm(A*x-b); % Error function

% Deterministic gradient
gradx = @(x) (AA*x-Ab);


% % This corresponds to a (kind of) stochastic gradient
% std_grad = sqrt(dim);
% gradx = @(x) (AA*x-Ab).*(1+std_grad*randn(dim,1)); 


% % Another noise for gradient (uniform)
% std_grad = sqrt(dim);
% gradx = @(x) (AA*x-Ab).*(1-std_grad+2*std_grad*rand(dim,1)); 

%% Definition of the algorithm

algorithm = @(x) x-(1/L)*gradx(x); % the used algorithm here is the basic gradient descend

%% The two loops

memory = zeros(dim,k+1); % contains the sequence to accelerate

nite_mainloop = 500; % Number of call of extrapolation algorithm
% Remark : the number of oracle calls is nite_mainloop*k.

error_vec = zeros(1,1+nite_mainloop); % for the plot

warning('off','MATLAB:nearlySingularMatrix')

x_exrapolated = x0;
error_vec(1) = error_fun(x_exrapolated);
for i=1:nite_mainloop
    memory(:,1) = x_exrapolated; % We start at the previous extrapolated point
    for j=1:k
        memory(:,j+1) = algorithm(memory(:,j)); % memory contains k+1 iterates of gradient descend
    end
    x_exrapolated = rmpe_blackbox(fx,memory,RmpeParam); % Extrapolation of k+1 iterates
    error_vec(i+1) = error_fun(x_exrapolated); % square norm of the residue
end

warning('on','MATLAB:nearlySingularMatrix')

%% Comparison with classic gradient method

error_grad = zeros(1,1+nite_mainloop);
x_grad = x0;
error_grad(1) = error_fun(x_grad);
for i=1:nite_mainloop
    for j=1:k
        x_grad = algorithm(x_grad);
    end
    error_grad(i+1) = error_fun(x_grad);
end

%% Plot
figure
semilogy(error_grad,'r')
hold on
semilogy(error_vec,'b')
legend('Gradient', 'Acceleration of gradient')

