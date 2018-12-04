function [ finfo ] = getFunction( type , param)
%GETFUNCTION Summary of this function goes here
%   Detailed explanation goes here

if strcmpi(type, 'Quadratic')
    n = param.n;
    A = rand(n);
    AA = A'*A;
    
    finfo.L = eigs(AA,1);
    finfo.mu = eigs(AA,1,'sm'); % sm for smallest eigenvalues
    
    finfo.f = @(x) x'*AA*x/2;
    finfo.fp = @(x) AA*x;
    finfo.fpp = @(x) AA;
    finfo.n = n;
    finfo.x0 = param.x0;
    
    finfo.xstar = zeros(n,1);
    finfo.fstar = finfo.f(finfo.xstar);
end


if strcmpi(type, 'LeastSquare')
    
    A = param.A;
    xstar = param.xstar;
    b = A*xstar;
    
    AA = A'*A;
    Ab = A'*b;
    I = eye(size(AA));
    n = size(AA,1);
    
    finfo.L = eigs(AA,1);
    finfo.mu = eigs(AA,1,'sm'); % sm for smallest eigenvalues
    
    finfo.fp = @(x) AA*x-Ab;
%     finfo.f = @(x) 0.5*(norm(A*x-b)^2 + lambda*norm(x)^2);
    finfo.f = @(x) max(abs((A*x-b))); % instead of real f, we compute the residual
    finfo.fpp = @(x) AA ;
    finfo.n = n;
    finfo.x0 = param.x0;
    
    finfo.xstar = param.xstar;
    finfo.fstar = 0;
end



if strcmpi(type, 'LeastSquare2')
    
    A = param.X';
    b = param.y;
    lambda = param.lambda;
    
    AA = A'*A;
    I = eye(size(AA));
    AA = AA + lambda*I;
    
    Ab = A'*b;
    finfo.xstar = AA\Ab;
    n = size(AA,1);
    
    finfo.L = norm(AA);
    finfo.mu = lambda; % sm for smallest eigenvalues
    
    finfo.fp = @(x) AA*x-Ab;
    finfo.A = AA;
    finfo.b = Ab;
%     finfo.fp = @(x) A'*(A*x-b)+lambda*x;
%     finfo.f = @(x) 0.5*(norm(A*x-b)^2 + lambda*norm(x)^2);
    finfo.f = @(x) norm(A*x-b)^2 + lambda*norm(x)^2; % instead of real f, we compute the residual
    finfo.fpp = @(x) AA ;
    finfo.n = n;
    finfo.x0 = param.x0;
    
    finfo.lscoeffun = @(x,p) (p'*finfo.fp(x))/(p'*AA*p);
    finfo.lsfun = @(x,p) x - finfo.lscoeffun(x,p)*p;
    
    finfo.lscoefnormx = @(x,p) (p'*(x-finfo.xstar))/(p'*p);
    
    finfo.lscoefnormgrad = @(x,p) (p'*AA*finfo.fp(x))/(norm(AA*p))^2;
    finfo.lsnormgrad = @(x,p) x - finfo.lscoefnormgrad(x,p)*p;
    finfo.fstar = finfo.f(finfo.xstar);
end


if(strcmpi(type,'normx4'))
    n = param.n;
    A = rand(n);
    A = A'*A;
    
    finfo.xstar = param.xstar;
    finfo.x0    = param.x0;
    
    xstar = param.xstar;
    finfo.f     = @(x)  0.5*((x-xstar)'*A*(x-xstar))^2;
    finfo.fp    = @(x)  ((x-xstar)'*A*(x-xstar))*A*(x-xstar);
    
    finfo.fpp   = @(x)  ((x-xstar)'*A*(x-xstar))*A + 2*(A*(x-xstar))*(A*(x-xstar))';
    
    finfo.fstar = finfo.f(finfo.xstar);
    finfo.n = n;
    finfo.L = norm(finfo.fpp(finfo.x0));
    finfo.mu = 0;
    
end

if(strcmpi(type,'logistic'))
    
    X = param.X;
    y = param.y;
    lambda = param.lambda;
    [nfeat,npt] = size(X);
    yXt = (X*spdiags(y,0,npt,npt))';
    
    finfo.f = @(w) (1/npt)*sum(  log(1+exp(-(yXt*w) )  )) + (lambda/2)*norm(w)^2;
    finfo.fp = @(w) -(1/npt)*sum( X*spdiags(y./(1+exp(yXt*w)),0,npt,npt) ,2) + lambda*w;
    
    
    % We multiply everything by npt to make it consistent with the notation
    % (1/npt) * sum_i f_i(x)
    
    finfo.frand = @(w,i) (log(1+exp(-( y(i)*X(:,i)'*w ) )  ) + (lambda/2)*norm(w)^2); 
%     finfo.fprand = @(w,i) npt*(( - X(:,i)*y(i)/(1+exp(y(i)*X(:,i)'*w)) ) + (lambda/npt)*w); 
    % Now can handle batches
    finfo.fprand = @(w,i) (1/length(i))*(  -sum( X(:,i)*spdiags(y(i)./(1+exp(yXt(i,:)*w)),0,length(i),length(i)) ,2)  + lambda*w); 
    
    if(~isfield(param,'approxL') || ~param.approxL)
        
        Hessianbound = X*X';
        
    else    
        
        Hessianbound = zeros(nfeat,1);
        for j=1:npt
            Hessianbound = Hessianbound + X(:,j).^2;
        end
        Hessianbound = max(Hessianbound)*nfeat;
        
    end
    
%     Hessianbound_indiv = 0;
%     for j=1:npt
%         Hessianbound_indiv = max(Hessianbound_indiv, norm(X(:,j))^2);
%     end
    
    finfo.L = (1/npt)*norm(full(Hessianbound))/4+lambda;
    finfo.L_max_sample = compute_Lmax_logistic(X);
    finfo.mu = lambda;
    
    finfo.xstar = nan(size(param.x0));
    finfo.fstar = 0;
    finfo.n = size(X,1);
    finfo.x0 = param.x0;
    finfo.nterms = npt;
    
end

if strcmpi(type, 'QuadPreDefined')
    A = param.A;
    
    finfo.xstar = param.xstar;
    finfo.x0 = param.x0;
    finfo.n = size(A,1);
    
    finfo.L = param.L;
    finfo.mu = param.mu; 
    
    finfo.f = @(x) (x-finfo.xstar)'*A*(x-finfo.xstar)/2;
    Axstar = A*finfo.xstar;
    finfo.fp = @(x) A*x-Axstar;
%     finfo.fp = @(x) A*(x-finfo.xstar);
    finfo.fpp = @(x) A;
    
    finfo.fstar = finfo.f(finfo.xstar);
end

if strcmpi(type, 'Eye')
    n = param.n;
    
    AA = eye(n);
    
    AA(1,1) = param.m*AA(1,1);
    
    finfo.L = eigs(AA,1);
    finfo.mu = eigs(AA,1,'sm'); % sm for smallest eigenvalues
    
    finfo.f = @(x) x'*AA*x/2;
    finfo.fp = @(x) AA*x;
    finfo.fpp = @(x) AA;
    finfo.n = n;
    finfo.x0 = param.x0;
    
    finfo.xstar = zeros(n,1);
    finfo.fstar = finfo.f(finfo.xstar);
end

if strcmpi(type, 'Quadsine')
    
    mu = param.mu;
    sigma = param.sigma;
    L=param.L;
    
    h1 = 2*(sigma-mu)/L;
    h2 = (sigma-mu)/2;
    
    f = @(x) -(h2/h1^2)*sin(h1*x)+(h2+mu)*(x.^2)/2+h2/h1;
    fp = @(x) -(h2/h1)*cos(h1*x)+(h2+mu)*x+h2/h1;
    fpp = @(x) h2*sin(h1*x)+h2+mu;
    
%     finfo.f = @(x) f(x(1)) + f(x(2));
%     finfo.fp = @(x) [fp(x(1)) ; fp(x(2))];
%     finfo.fpp = @(x) [fpp(x(1)) 0 ; 0 fpp(x(2))];
    
    finfo.f = @(x) f(x);
    finfo.fp = @(x) fp(x);
    finfo.fpp = @(x) fpp(x);
    
    finfo.n = 1;
    finfo.L = sigma;
    finfo.mu = mu;
    finfo.x0 = param.x0(1:finfo.n);
    
    finfo.xstar = zeros(finfo.n,1);
    finfo.fstar = finfo.f(finfo.xstar);
end

if(strcmpi(type, 'Lasso'))
    A = param.A;
    b = param.b;
    lambda = param.lambda;
    
    finfo.L = norm(A)^2; % i.e. norm (A'A)
    finfo.mu = 0;
    
    finfo.f = @(x) 0.5*norm(A*x-b)^2;
    finfo.fp = @(x) A'*(A*x-b);
    finfo.fpp = @(x) nan;
    
    finfo.proxoperator.f = @(x) lambda*norm(x,1);
    finfo.proxoperator.fp = @(x) lambda*sign(x,1);
    finfo.proxoperator.prox = @(x,h) sign(x).*max(abs(x)-h*lambda,0);
    finfo.proxoperator.mu = 0;
    
    finfo.n = size(param.A,2);
    finfo.x0 = param.x0;
    
    finfo.xstar = param.xstar;
    finfo.fstar = finfo.f(finfo.xstar);
    finfo.proxoperator.fstar = finfo.proxoperator.f(finfo.xstar);
end


if(strcmpi(type, 'DualSvm'))
    X = param.X;
    y = param.y;
    
    yX = diag(y)*X;
    
    finfo.L = norm(yX)^2; % i.e. norm (A'A)
    finfo.mu = 0;
    
    
    finfo.proxoperator.f = @(z) 0;
    finfo.proxoperator.fp = @(z) nan;
    finfo.proxoperator.prox = @(z,h) min(max(z,0),1);
    finfo.proxoperator.mu = 0;
    
    finfo.f = @(z) 0.5*norm(yX'*(finfo.proxoperator.prox(z)))^2-sum(finfo.proxoperator.prox(z));
    finfo.fp = @(z) yX*(yX'*z)-1;
    finfo.fpp = @(z) nan;
    
    finfo.n = length(y);
    finfo.x0 = param.x0;
    
    finfo.xstar = param.xstar;
    finfo.fstar = finfo.f(finfo.xstar);
    finfo.proxoperator.fstar = finfo.proxoperator.f(finfo.xstar);
end

if(strcmpi(type,'MaxCut'))
    
    L = param.laplacian;
    n = size(L,1);
    
    finfo.xstar = param.xstar;
    finfo.fstar = param.fstar;
    finfo.x0 = param.x0;
    
    finfo.proxoperator.f = @(z) 0;
    finfo.proxoperator.fp = @(z) nan;
    finfo.proxoperator.prox = @(z,h) z*0;
    finfo.proxoperator.mu = 0;
    finfo.proxoperator.fstar = 0;
    
    finfo.L = nan;
    finfo.mu = 0;
    finfo.n=n;
    
    finfo.f = @(z) -sum(z)/n+eigs(L+diag(z),1,'la');
    finfo.fp = @(z) -ones(size(z))/n+maximal_eigenvector(L+diag(z)).^2; % differential of lambda max is the maximal eigenvector
    finfo.fpp = @(z) nan;
    
end

if strcmpi(type, 'Smooth-abs')
    
    n = param.n;
%     D = abs(diag(rand(n)));
    D = ones(n,1);
    
    maxd = max(D);
    
    finfo.L = maxd^2;
    
    finfo.f = @(x) sum(log(exp(D'*x)+exp(-D'*x)));
    finfo.fp = @(x) D.*(exp(2*D.*x)-1)./(exp(2*D.*x)+1);
    finfo.fpp = @(x) diag((4*D.^2).*exp(2*D.*x)./((1+exp(2*D.*x)).^2));
    finfo.n = n;
    finfo.x0 = param.x0;
    
    finfo.xstar = zeros(n,1);
    finfo.fstar = finfo.f(finfo.xstar);
    
end


if strcmpi(type, 'logsumexp')
    n = param.n;  
    [finfo.f, finfo.fp, finfo.fpp, finfo.L, finfo.xstar] = logsumexp_diag(n);
    finfo.n = n;
    finfo.x0 = param.x0;
    finfo.mu = 0;
    finfo.fstar = finfo.f(finfo.xstar);
end

if strcmpi(type, 'Self-concordant')
    
    x0 = param.x0;
    mu = param.mu;
    n = length(x0);
    finfo.n = n;
    
    % box constraint: x <= 2*x0 ; -x <= 2*x0 ;
    Abox = [eye(n) ; -eye(n)];
    bbox = 1.0001*[abs(x0) ; abs(x0)];
    
    % Random constraints
%     A = rand(n);
%     b = rand(n,1);
%     
%     % makes constraints feasible for x0
%     unfeasConstr = A*x0-b > 0; % Ax shoudl be <= b
%     A(unfeasConstr,:) = -A(unfeasConstr,:);
%     b(unfeasConstr) = -b(unfeasConstr);
    A = [];
    b = [];
    
    AConstr = [A ; Abox];
    bConstr = [b ; bbox];
    % Check if everything is ok
    if(any(AConstr*x0-bConstr>0))
        error('Problem: x0 not feasible')
    end
    
    
    finfo.f = @(x) mu*sum(log(bConstr-AConstr*x));
    finfo.fp = @(x) mu*(AConstr'*(1./(bConstr-AConstr*x)));
    finfo.fpp = @(x) mu*fpp_selfconcordant(x,AConstr,bConstr);
    
    finfo.L = norm(finfo.fpp(x0));
    finfo.x0 = param.x0;
    
    finfo.xstar = zeros(n,1);
    finfo.fstar = finfo.f(finfo.xstar);
end

%% Additionnal functions

    function Accu = fpp_selfconcordant(x,A,b)
        
        Accu = zeros(size(A,2));
        for i=1:length(b)
            Accu = Accu + A(i,:)'*A(i,:)*1/((b(i)-A(i,:)*x)^2);
        end
    end

    function out = proj_cone_definite_positive(S)
        % S is assumed symmetric
        [EigenVec,EigenVal] = eig(S);
        EigenVal = max(EigenVal,0);
        out = EigenVec*EigenVal*EigenVec';
        
    end

    function out = minimal_eigenvector(x)
        [V,~] = eigs(x,1,'sa');
        out = V(:,1);
        out = out/norm(out);
    end

    function out = maximal_eigenvector(x)
        [V,~] = eigs(x,1,'la');
        out = V(:,1);
        out = out/norm(out);
    end

end

