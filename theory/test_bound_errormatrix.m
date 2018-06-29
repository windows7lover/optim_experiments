% 

sigmavec = linspace(0.01,10,25);
Sigma = diag(sigmavec);
Sigmasquare = diag(sigmavec.^2);
Sigmainv = diag(1./sigmavec);

lambda = 1;
normP = 0.1;

n = length(sigmavec);

I = eye(n);

cvx_begin

cvx_precision high

variable M(n,n)

    minimize ( norm(((Sigmasquare+M)*Sigmainv)) )
subject to
    (M + Sigmasquare) - lambda*I == semidefinite(n);
    norm(M) <= lambda+normP;

cvx_end

realbound = 1/cvx_optval;

bound1 = max(sigmavec/lambda);
bound2 = sqrt(2)/(sqrt(lambda+normP));
[bound1, bound2, min(bound1,bound2), realbound]