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
if(~isfield(param,'steplength'))
    param.steplength = 1/finfo.L; %zeros(length(x),param.window_size+1);
end
if(~isfield(param,'L'))
    param.L = finfo.L;
end
% if(~isfield(param,'Lstack'))
%     param.Lstack = param.L;
% end
if(~isfield(param,'y_seq'))
    param.y_seq = [];
end

if(~isfield(param,'dorna'))
    param.dorna = 1;
end


if(~isfield(param,'alpha_k'))
%     if(param.strong_convex)
%         param.alpha_k = sqrt(mu/finfo.L);
%     else
%         param.alpha_k = 1;
%     end
    param.alpha_k = 1;
end

% if(param.strong_convex)
% %     beta = (sqrt(finfo.L)-sqrt(mu))/(sqrt(finfo.L)+sqrt(mu));
%     q = mu/finfo.L;
% else
%     alpha_new = roots([1, param.alpha_k^2, -param.alpha_k^2]);
% end
if(param.strong_convex)
    q = mu/finfo.L;
else
    q = 0;
end

alpha_new = roots([1, param.alpha_k^2-q, -param.alpha_k^2]);
alpha_new = alpha_new(alpha_new>0 & alpha_new<1);
beta = param.alpha_k*(1-param.alpha_k)/(param.alpha_k^2 + alpha_new);
param.alpha_k = alpha_new;

if(~param.accelerated)
    beta = 0;
end

if(k == 2)
    param.xold = yold;
end

% param.Lstack = param.L;

% backtracking line search
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
        if(param.L > finfo.L)
            param.L = finfo.L;
            break
        end
    end
%     param.Lstack(end+1) = param.L;
else
    xplus = yold - g/param.L;
    fxplus = finfo.f(xplus);
end

new_rna_point = yold;% - g/finfo.L;
if( size(param.x_seq,2) < param.window_size+1)
    param.steplength = [param.steplength, 1/param.L];
    param.x_seq = [param.x_seq, new_rna_point];
    param.U = [param.U g];
    param.y_seq = [param.y_seq yold];
else
    param.steplength = [param.steplength(:,2:end), 1/param.L];
    param.x_seq = [param.x_seq(:,2:end), new_rna_point];
    param.U = [param.U(:,2:end) g];
    param.y_seq = [param.y_seq(:,2:end) yold];
end

UU = param.U'*param.U;
% UU = param.U'*param.y_seq;
% UU = UU + UU';
UU = UU/norm(UU);


if(param.dorna)
    
    if(size(param.x_seq,2) <= 1)
        x_rna_out = xplus;
    else
        [ x_rna_out, c ] = ampe( param.x_seq, UU, param.lambda );
    end
    
    steplength = 1/finfo.L; % worst case choise
%     steplength = max(param.steplength); % pessimist choice
%     steplength = max(param.steplength); % optimist
    
    Uc = param.U*c;
    x_rna = x_rna_out - steplength*Uc;
    x_ls = (x_rna+beta*xold)/(1+beta);
    
    if(param.online)
        fls = finfo.f(x_ls);
        if( (fls < fxplus) )
            % agressive line search
            
%             fls_old = inf;
%             while(fls < fls_old)
%                 fls_old = fls;
%                 steplength = steplength*2;
%                 x_rna = x_rna_out - steplength*Uc;
%                 x_ls = (x_rna+beta*xold)/(1+beta);
%                 fls = finfo.f(x_rna);
% %                 break;
%             end
%             steplength = steplength/2;
            x_rna = x_rna_out - steplength*Uc;
            x_ls = (x_rna+beta*xold)/(1+beta);
            
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