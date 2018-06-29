function [xstar,ystar,fstar_primal,fstar_dual] = solve_maxcut_sdp(L)

n = size(L,1);

% % Primal problem
% cvx_begin
% 
%     variable X(n,n)
%     minimize -trace(L*X)
%     subject to 
%     X == semidefinite(n,n);
%     diag(X) == ones(n,1)/n;
%     trace(X) == 1;
% 
% cvx_end

% xstar = X;
% fstar_primal = cvx_optval;
xstar = [];
fstar_primal = [];


% Dual problem
cvx_begin

    variables y(n,1)
    minimize -sum(y)/n + lambda_max(L+diag(y))
    subject to
    L-diag(y) == semidefinite(n,n);
    
cvx_end

ystar = y;
fstar_dual = cvx_optval;