% 

L = 10;
mu = 1;
n = 10;
type = 'Random';

k = 5;
x0 = rand(n,1);
x0 = x0/norm(x0);

sigma = 1-mu/L;

opt.L = sigma;

G = generateRandMatrix(n,type,opt);
% G = eye(n)*sigma;
 
X = linear_fixed_point(G,0*x0,x0,k);
Xbar = X - repmat(mean(X,2),1,size(X,2));

display([norm(X) norm(Xbar)])