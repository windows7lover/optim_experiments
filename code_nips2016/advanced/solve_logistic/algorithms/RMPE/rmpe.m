function [ xout, c ] = rmpe( x, UU_norm, lambda )

if(nargin<2)
    lambda = 0;
end

% In case we cannot compute U
if(size(x,2)==1)
    xout = x;
    return;
end

% U = diff(x,1,2);
k = size(UU_norm,2);

% % Compute minimal polynomial coefficients
% UU = U'*U;

% if(norm(UU)>100*eps)
%     UU = UU/norm(UU);
% else
%     UU = 0*UU;
% end

matrix = (UU_norm + eye(k)*lambda);

c = matrix\ones(k,1);
c = c/sum(c);

xout=x(:,1:end-1)*c;
