function [ xout, c ] = rmpe( x, UU_norm, lambda )
%	Solve the problem 
%
%	c^* = argmin_c |Uc|^2 + lambda*|c|^2
%	s.t. sum(c) == 1
%
% where U(:,i) = x(:,i+1) - x(:,i)
% and UU_norm = U'*U/norm(U'*U).
%
% output the value sum_i=1^k c^*_i x_i.

if(nargin<2)
    lambda = 0;
end

if(size(UU_norm,2)==0)
    xout = x;
    return;
end

% U = diff(x,1,2);
k = size(UU_norm,2);

matrix = (UU_norm + eye(k)*lambda);

c = matrix\ones(k,1);
c = c/sum(c);

xout=x(:,1:end-1)*c;
