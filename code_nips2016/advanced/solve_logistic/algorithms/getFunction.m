function [ finfo ] = getFunction( type , param)

if(strcmpi(type,'logistic'))
    
    X = param.X;
    y = param.y;
    lambda = param.lambda;
    [nfeat,npt] = size(X);
    yXt = (X*spdiags(y,0,npt,npt))';
    
    finfo.f = @(w) sum(  log(1+exp(-(yXt*w) )  )) + (lambda/2)*norm(w)^2;
    finfo.fp = @(w) -sum( X*spdiags(y./(1+exp(yXt*w)),0,npt,npt) ,2) +lambda*w;
    
    if(~isfield(param,'approxL') || ~param.approxL)
        
        Hessianbound = X*X';
        
    else
        % Another way to compute L, leads to similar performances but does
        % not require too much memory
        Hessianbound = zeros(nfeat,1);
        for j=1:npt
            Hessianbound = Hessianbound + X(:,j).^2;
        end
        Hessianbound = max(Hessianbound)*nfeat;
        
    end
    finfo.L = norm(full(Hessianbound))/4+lambda;
    finfo.mu = lambda;
    
    finfo.xstar = nan(size(param.x0));
    finfo.fstar = 0;
    finfo.n = size(X,1);
    finfo.x0 = param.x0;
    
end

end