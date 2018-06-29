function [solution,algo_tolf,algo_tolx,algoparam,time_ite] = dokiter_minfunc(algoparam,~,nIteMax,finfo,tolfun,nIterTol)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Init do k iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isfield(finfo,'proxoperator'))
    finfo.proxoperator.f = @(x) 0;
    finfo.proxoperator.fp = @(x) 0;
    finfo.proxoperator.mu = 0;
    finfo.proxoperator.fstar = 0;
end
valf_handle = @(x) finfo.f(x) + finfo.proxoperator.f(x);
valfp_handle = @(x) finfo.fp(x) + finfo.proxoperator.fp(x);

if(nargin<5)
    tolfun = struct();
end

% We always compute the first one
if(nargin<6)
    nIterTol = nIteMax;
end

temp = linspace(1,nIteMax+1,nIterTol+1);
itertol = floor(temp);

if(any(temp~=itertol))
    itertol
    temp
    warning('nIteMax/nIterTol not an integer, rounding...')
    [nIteMax nIterTol]
end

nTol = length(itertol); % we compute x0

if( isnan(finfo.fstar) )
    error('cannot approximate fstar here')
end
fstar = finfo.fstar;

nIteMax = nIteMax+1; % "iteration" 1 = x0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Use minFunc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if( ~isfield(algoparam,'minFuncOpt'))
    error('No algo specified');
end

if( isfield(algoparam.minFuncOpt,'Method'))
    options.Method = algoparam.minFuncOpt.Method; % example : 'lbfgs'
end

if( isfield(algoparam.minFuncOpt,'LS_init'))
    options.LS_init = algoparam.minFuncOpt.LS_init;
end

if( isfield(algoparam.minFuncOpt,'LS_interp'))
    options.LS_interp = algoparam.minFuncOpt.LS_interp;
end

if( isfield(algoparam.minFuncOpt,'LS'))
    options.LS = algoparam.minFuncOpt.LS;
end


if( isfield(algoparam.minFuncOpt,'Corr'))
    options.Corr = algoparam.minFuncOpt.Corr; % memory size, default 100
end


options.Display = 'off';
options.MaxIter = nIteMax;
options.maxFunEvals = inf;
options.optTol = 0;
options.progTol = 0;

algoparam.minFuncOpt = options;
algoparam.trace = algoparam;


funObj = @(x) convert_finfo(x);
x0 = finfo.x0;


t_beg = cputime;
[solution,~,exitflag,output] = minFunc(funObj,x0,options);
display(exitflag)
t_end = cputime;

algoparam.output = output;

total_time = t_end-t_beg;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trace = output.trace;
n = length(trace.fval);
time_ite = linspace(0,total_time,n);
algo_tolf = trace.fval';
if(~isnan(fstar))
    algo_tolf = algo_tolf - fstar;
end
algo_tolx = NaN(n,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Conversion of finfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [obj,grad] = convert_finfo(x)
        
        valf_handle = @(x) finfo.f(x) + finfo.proxoperator.f(x);
        valfp_handle = @(x) finfo.fp(x) + finfo.proxoperator.fp(x);
        
        obj = valf_handle(x);
        grad = valfp_handle(x);
%         hess = finfo.fpp(x); % should be something else
        
    end

end