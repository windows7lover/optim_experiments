%% example.m

% create a random quadratic function
d = 100;
Q = rand(d);
Q = Q'*Q;

xstar = randn(d,1);
x0 = randn(d,1);

fx = @(x) (x-xstar)'*Q*(x-xstar)/2;
nablafx = @(x) Q*(x-xstar);
ls = @(x,p) -p'*Q*(x-xstar)/(p'*Q*p);

% Define gradient descent
L = norm(Q); % Stepsize
g = @(x) x-nablafx(x)/L;

% define offline and online RNA
N = 10; % Memory size
lambda = 1e-8;
beta = 1;

% Choose if we use line search
% param.linesearch = NaN; % No line search
param.linesearch = @(x,p) ls(x,p);  % with line search (override beta)


offlineRNA = RNA(N,lambda,beta,param);
onlineRNA = RNA(N,lambda,beta,param);

% Storage for function values over iterations
nIter = 100;
f_grad = zeros(1,nIter);
f_rna_offline = zeros(1,nIter);
f_rna_online = zeros(1,nIter);



% Gradient descent
x=x0;
y=x0;
for i=1:nIter
    f_grad(i) = fx(y);
    x = g(y);
    y=x;
end


% Gradient descent + offline RNA
x=x0;
y=x0;
yextr = x0;
for i=1:nIter
    f_rna_offline(i) = fx(yextr);
    x = g(y);
    offlineRNA = offlineRNA.store(x,y);
    y=x;
    yextr = offlineRNA.accelerate();
end


% Gradient descent + online RNA
x=x0;
y=x0;
for i=1:nIter
    f_rna_online(i) = fx(y);
    x = g(y);
    onlineRNA = onlineRNA.store(x,y);
    y = onlineRNA.accelerate();
end

%%
semilogy(f_grad,'linewidth',2)
hold on
semilogy(f_rna_offline,'linewidth',2)
semilogy(f_rna_online,'linewidth',2)
legend({'Gradient Method', 'Gradient + offline RNA', 'Gradient + online RNA'},'location','SW','box','off')
set(gca,'fontsize',16)
