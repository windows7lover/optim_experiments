function [algo_x,algo_tolf,algo_tolx,algoparam,time_ite] = do_k_iterations(algoparam,algohandle,nIteMax,finfo,tolfun,nIterTol)

if(~isfield(finfo,'proxoperator'))
    finfo.proxoperator.f = @(x) 0;
    finfo.proxoperator.fp = @(x) 0;
    finfo.proxoperator.mu = 0;    
    finfo.proxoperator.fstar = 0;
end
valf_handle = @(x) finfo.f(x) + finfo.proxoperator.f(x);
valfp_handle = @(x) finfo.fp(x) + finfo.proxoperator.fp(x);
mu = finfo.mu+finfo.proxoperator.mu;


    
compute_tol = true;
if(nargin<5)
    tolfun = struct();
end

% We always compute the first one
if(nargin<6 || isempty(nIterTol) || isnan(nIterTol) )
    nIterTol = nIteMax;
end

temp = linspace(1,nIteMax+1,nIterTol+1);
itertol = floor(temp);

if(any(temp~=itertol))
    warning('nIteMax/nIterTol not an integer, rounding...')
    display('nIteMax/nIterTol not an integer, rounding...')
%     itertol
%     temp
%     [nIteMax nIterTol]
end

nTol = length(itertol); % we compute x0

if( (~isstruct(tolfun)) && ((any(isnan(tolfun)) || isempty(tolfun))) )
    compute_tol = false;
    need_approx_fstar = false;
else
    if(isfield(tolfun,'tolf'))
        need_approx_fstar = false;
        tolf_handle = tolfun.tolf;
    else
        need_approx_fstar = false;
        if(~isnan(finfo.fstar))
            tolf_handle = @(x) valf_handle(x) - finfo.fstar - finfo.proxoperator.fstar ;
        elseif(isnan(finfo.fstar) && mu>0)
            tolf_handle = @(x) valf_handle(x);
            need_approx_fstar = true;
            approx_fstar_handle = @(x) valf_handle(x)  - (1/(2*mu))*norm(valfp_handle(x))^2;
            approx_fstar = approx_fstar_handle(finfo.x0) ;
        end
    end
    if(isfield(tolfun,'tolx'))
        tolx_handle = tolfun.tolx;
    else
        tolx_handle = @(x) norm(x-finfo.xstar);
    end
end

nIteMax = nIteMax+1; % "iteration" 1 = x0

if(compute_tol)
    algo_tolf = zeros(1,nTol);
    algo_tolx = zeros(1,nTol);
    time_ite = zeros(1,nTol);
    

    algo_tolf(1) = tolf_handle(finfo.x0);
    algo_tolx(1) = tolx_handle(finfo.x0);
else
    algo_tolf = NaN;
    algo_tolx = NaN;
    time_ite = NaN;
end

algo_x = zeros(finfo.n,nTol);
algo_x(:,1) = finfo.x0;

x_new = finfo.x0;

if(compute_tol)
    totaltime_begin = clock;
end

theitertol = 2; % start at 2, because first entry is 1

t_beg = cputime;
for k=2:nIteMax
    [x_new,algoparam,fval] = algohandle(finfo,x_new,k,algoparam);
    
    if( all( [theitertol <= length(itertol) , itertol(theitertol) == k] ) )
        
        algo_x(:,theitertol) = x_new;
    
        if( all( [compute_tol , theitertol <= length(itertol) , itertol(theitertol) == k] ))

            t_end = cputime;        
            time_ite(theitertol) = time_ite(theitertol-1) + (t_end-t_beg);
            
            if( ~isnan(fval) )
                algo_tolf(k) = fval  - finfo.fstar - finfo.proxoperator.fstar ;
            else
                algo_tolf(k)  = tolf_handle(x_new);
            end

            algo_tolx(:,theitertol) = tolx_handle(x_new);

            t_beg = cputime; % reset timer
        end
        theitertol = theitertol+1;
    end
    
    
    if(need_approx_fstar)
        approx_fstar = max(approx_fstar,approx_fstar_handle(algo_x(:,k)) );
    end
end

if(compute_tol) % cputime not really accurate, so we normalize using the totaltime
    totaltime_end = clock;
    factor = min( etime(totaltime_end,totaltime_begin)/time_ite(end) , 1);
    time_ite = time_ite*factor;
end


if(need_approx_fstar)
    algo_tolf = algo_tolf-approx_fstar;
end
