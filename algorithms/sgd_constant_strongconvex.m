function [xplus, param] = sgd_constant_strongconvex(finfo,x,k,param)

n = finfo.nterms;

if(~isfield(param,'batchsize'))
    param.batchsize = 1;
end

% randIdx = randi(n,param.batchsize);
if(param.batchsize ~= 1)
    randIdx = datasample(1:n,param.batchsize,'Replace',false);
else
    randIdx = randi(n,1);
end

gRand = finfo.fprand(x,randIdx(1));

for i=2:length(randIdx)
    gRand = gRand+finfo.fprand(x,randIdx(i));
end

gRand = gRand/param.batchsize;

gamma = 1/finfo.L; % Stepsize
gamma = gamma*(param.batchsize/n);

xplus = x-gamma * gRand ;

