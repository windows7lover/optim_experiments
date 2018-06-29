function [xout,param,fxplus] = rna_first_order(finfo,x,k,param)

mu = finfo.mu;
fp = finfo.fp;
yold = x;

% init
if(~isfield(param,'online'))
    param.online = true;
end
if(~isfield(param,'lambda'))
    param.lambda = 1e-8;
end
if(~isfield(param,'U'))
    param.U = []; %zeros(length(x),param.window_size);
end
if(~isfield(param,'x_seq'))
    param.x_seq = x; %zeros(length(x),param.window_size+1);
end
if(~isfield(param,'L'))
    param.L = finfo.L;
end

if(~isfield(param,'dorna'))
    param.dorna = 1;
end


if(~isfield(param,'alpha_k'))
    param.alpha_k = 1;
end

if(param.strong_convex)
    beta = (sqrt(finfo.L)-sqrt(mu))/(sqrt(finfo.L)+sqrt(mu));
else
    alpha_new = roots([1, param.alpha_k^2, -param.alpha_k^2]);
    alpha_new = alpha_new(alpha_new>0 & alpha_new<1);
    beta = param.alpha_k*(1-param.alpha_k)/(param.alpha_k^2 + alpha_new);
    param.alpha_k = alpha_new;
end

if(~param.accelerated)
    beta = 0;
end

if(k == 2)
    param.xold = yold;
end

xold = param.xold;
g = fp(yold);
if(param.backtracking)
    fold = finfo.f(yold);
    normg = norm(g);
    xplus = yold - g/param.L;
    fxplus = finfo.f(xplus);
    while(fxplus > fold - (1/(2*param.L))*normg^2) % descent condition
        param.L = param.L*2;
        xplus = yold - g/param.L;
        fxplus = finfo.f(xplus);
    end
else
    xplus = yold - g/param.L;
    fxplus = finfo.f(xplus);
end

new_rna_point = yold - g/finfo.L;
if( size(param.x_seq,2) < param.window_size+1)
    param.x_seq = [param.x_seq, new_rna_point];
    param.U = [param.U g];
else
    param.x_seq = [param.x_seq(:,2:end), new_rna_point];
    param.U = [param.U(:,2:end) g];
end
UU = param.U'*param.U;
UU = UU/norm(UU);
% UU = UU/(norm(UU)^(3/2));



if(param.dorna)
    
    if(size(param.x_seq,2) <= 1)
        x_rna = xplus;
    else
        [ x_rna, c ] = ampe( param.x_seq, UU, param.lambda );
    end
    x_ls = (x_rna+beta*xold)/(1+beta);
    
    if(param.online)
        fls = finfo.f(x_ls);
        if( (fls < fxplus) )
            z = x_ls;
            fxplus = finfo.f(x_rna);
        else
            z = xplus;
        end
    else
        z = xplus;
        fxplus = finfo.f(x_rna);
    end
else
    z = xplus;
end

yplus = (1+beta)*z - beta*xold;
param.xold = z;

xout = yplus;

if(param.backtracking)
    param.L = param.L/2;
end