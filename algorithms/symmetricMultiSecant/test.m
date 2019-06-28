
d = 100;
nIter = 100;
N = 10;

A = rand(d);
A = A'*A;
A = A/norm(A);

x0 = rand(d,1);
xstar = rand(d,1);

%%

f = @(x) (x-xstar)'*A*(x-xstar)/2;
g = @(x) A*(x-xstar);
accum_bfgs = Accumulator(d,N);
accum_bfgs_quad = Accumulator(d,N);
accum_bfgs_quad_precond = Accumulator(d,N);
accum_nlSymQN = Accumulator(d,N);
accum_rna = Accumulator(d,N);


err_fun = @(x) f(x);
% err_fun = @(x) norm(g(x));

err_grad = zeros(nIter,1);
err_bfgs = zeros(nIter,1);
err_bfgs_quad = zeros(nIter,1);
err_bfgs_quad_precond = zeros(nIter,1);
err_bfgs_nlSymQN = zeros(nIter,1);
err_rna = zeros(nIter,1);


% Gradient
y=x0;
for i=1:nIter
    x = y-g(y);
    err_grad(i) = err_fun(y) ;
    y=x;
end


% BFGS
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_bfgs = accum_bfgs.store(y,x-y);
    qNstep = getSymHessian(accum_bfgs,-1,g(y));
    err_bfgs(i) = err_fun(y) ;
    y = y+qNstep;
end



% BFGS_quad
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_bfgs_quad = accum_bfgs_quad.store(y,x-y);
    qNstep = quad_bfgs(accum_bfgs_quad,-1,g(y));
    err_bfgs_quad(i) = err_fun(y) ;
    y = y+qNstep;
end


% BFGS_quad_precond
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_bfgs_quad_precond = accum_bfgs_quad_precond.store(y,x-y);
    qNstep = quad_bfgs_precond(accum_bfgs_quad_precond,-1,g(y));
    err_bfgs_quad_precond(i) = err_fun(y) ;
    y = y+qNstep;
end


% BFGS nlSymQN
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_nlSymQN = accum_nlSymQN.store(y,x-y);
    [qNstep] = nonLinearPrecondQN(accum_nlSymQN,-1,g(y));
    err_bfgs_nlSymQN(i) = err_fun(y) ;
    y = y+qNstep;
end


% RNA
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_rna = accum_rna.store(y,x-y);
    err_rna(i) = err_fun(y) ;
    y = rna(accum_rna);
end



semilogy(err_grad)
hold on
semilogy(err_bfgs,'linewidth',2)
semilogy(err_bfgs_quad,'linewidth',2)
semilogy(err_bfgs_quad_precond,'linewidth',2)
semilogy(err_bfgs_nlSymQN,'.','markersize',25)
semilogy(err_rna)

legend({'gradient','multisecant BFGS','multisecant BFGS quadratic','multisecant BFGS quadratic preconditionned','multisecant nlSymQN','rna'})