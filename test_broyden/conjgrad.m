function [tolx,tolgrad,tolf] = conjgrad(finfo,iter)
% CONJGRAD  Conjugate Gradient Method.
%   X = CONJGRAD(A,B) attemps to solve the system of linear equations A*X=B
%   for X. The N-by-N coefficient matrix A must be symmetric and the right
%   hand side column vector B must have length N.
%
%   X = CONJGRAD(A,B,TOL) specifies the tolerance of the method. The
%   default is 1e-10.
%
% Example (highlight lines between %{ and %}, then press F9):
%{
  n = 6000;
  m = 8000;
  A = randn(n,m);
  A = A * A';
  b = randn(n,1);
  tic, x = conjgrad(A,b); toc
  norm(A*x-b)
%}
% By Yi Cao at Cranfield University, 18 December 2008
% Updated on 6 Feb 2014.
%
A = finfo.A;
b = finfo.b;
x = finfo.x0;
r = b - A*x;
d = r;

tolx = zeros(iter+1,1);
tolgrad = zeros(iter+1,1);
tolf = zeros(iter+1,1);
k=0;
tolx(k+1) = norm(x-finfo.xstar);
tolgrad(k+1) = norm(finfo.fp(x));
tolf(k+1) = finfo.f(x)-finfo.fstar;

for k = 1:iter
    
    alpha = r'*r/(d'*A*d);
    x = x+alpha*d;
    r_new = -finfo.fp(x);
    beta = r_new'*r_new/(r'*r);
    r = r_new;
    d = r + beta*d;
    
    tolx(k+1) = norm(x-finfo.xstar);
    tolgrad(k+1) = norm(finfo.fp(x));
    tolf(k+1) = finfo.f(x)-finfo.fstar;
end
end