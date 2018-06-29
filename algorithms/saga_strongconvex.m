function [xplus, param] = saga_strongconvex(finfo,x,k,param)

n = finfo.nterms;

if(k == 2) % compute full gradient
%     param.meangrad = finfo.fp(x);
    param.meangrad = zeros(length(x),1);
    param.table = zeros(length(x),n);
    for i=1:n
        param.table(:,i) = finfo.fprand(x,i);
        param.meangrad = param.meangrad + param.table(:,i)/n;
    end
end

randIdx = randi(n,1);
gRand = finfo.fprand(x,randIdx);

% gamma = 1/(3*(finfo.mu+finfo.L)); % Stepsize
gamma = 1/(3*finfo.L_max_sample); % Stepsize

difference = gRand-param.table(:,randIdx);
xplus = x-gamma*( difference + param.meangrad );


param.meangrad = param.meangrad + difference/n;
param.table(:,randIdx) = gRand;

% fpx = finfo.fp(x);
% display( norm( param.meangrad - fpx)/norm(fpx) )
