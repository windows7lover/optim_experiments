function [algo_x,algo_tolf,algo_tolx,algoparam] = do_k_iterations(algoparam,algohandle,nIteMax,finfo,tolfun)

% DO_K_ITERATION
% Perform k steps of some user-defined algorithm algohandle
%
% Inputs: algoparam,algohandle,nIteMax,finfo,tolfun
%   - algoparam:    (structure) Parameters of the user-defined algorithm
%   - algohandle:   (function handle) User-defined algorithm. At each
%                   step, the function calls
%                   "[x,algoparam_new] = algohandle(finfo,algo_x,k,algoparam)"
%                   where x is the next iterate and algoparam_new is an
%                   updated version of algoparam.
%   - nIteMax:      (interger >0) Desired number of iterations
%   - fInfo:        (struture) The function to optimize. The structure of finfo
%                   should at least contains
%                       * finfo.f(x) (handle) the objective function value
%                       * finfo.fp(x) (handle) the gradient of the 
%                       objective function
%                       * finfo.x0 (vector) the initial value for the
%                       algorithm algohandle
%                   fInfo may also constain
%                       * finfo.fstar (float) Optimal value of f(x). 
%                       * finfo.mu (float) the strong convexity parameter,
%                       used to approximate fstar if fstar is NaN or not
%                       defined
%   - tolfun        (structure, optionnal) define the way the tolerances
%                   should be computed. The structure should contain
%                       * tolfun.tolf(x) (handle) information about
%                       f(x)-f(xstar)
%                       * tolfun.tolx(x) (handle) information about |x-xstar|
%                   by default, tolf(x) = f(x)-f(xstar) and 
%                   tolx = norm(x-xstar).
%                   if tolfun is NaN, then the tolerances are not computed.
%                   
% Outputs: algo_x,algo_tolf,algo_tolx,algoparam
%   - algo_x:       ( n*(nIteMax+1) matrix )  Contains all iterates produced 
%                   by algohandle
%   - algo_tolf:    (vector of size (nIteMax+1)) Contains all values of
%                   tolfun.tolf(algo_x(:,k)) for k=1:(nIteMax+1)
%   - algo_tolx:    (vector of size (nIteMax+1)) Contains all values of
%                   tolfun.tolx(algo_x(:,k)) for k=1:(nIteMax+1)
%   - algoparam:    (structure) last value of algoparam



% Some shortcuts
valf_handle = @(x) finfo.f(x);
valfp_handle = @(x) finfo.fp(x);


% Check all fields of finfo
if ~isfield(finfo,'mu')
    warning('Function do_k_iterations: no "mu" given, so finfo.mu = nan')
    mu = nan;
else
    mu = finfo.mu;
end

if ~isfield(finfo,'xstar')
    warning('Function do_k_iterations: no "xstar" given, so finfo.xstar = nan')
    finfo.xstar = nan(size(finfo.x0));
end

if ~isfield(finfo,'fstar')
    warning('Function do_k_iterations: no "fstar" given, so finfo.fstar = nan')
    finfo.fstar = nan;
end


if ~isfield(finfo,'n')
    warning('Function do_k_iterations: no "n" given, so finfo.n = length(finfo.x0)')
    finfo.n = length(finfo.x0);
end


% Check if we use the default tolerances computation

compute_tol= true;
if(nargin<5)
    tolfun = struct();
end


if( (~isstruct(tolfun)) && isnan(tolfun)) 
    % Here, we do not compute any kind of tolerance
    compute_tol = false;
    need_approx_fstar = false;
else
    % Determine how to compute tolerance on the function
    if(isfield(tolfun,'tolf'))
        % tolf is defined by the user
        need_approx_fstar = false;
        tolf_handle = tolfun.tolf;
    else
        % tolf is not defined
        need_approx_fstar = false;
        
        if(~isnan(finfo.fstar))
            % if fstar is defined, then tolf = f(x) -fstar
            tolf_handle = @(x) valf_handle(x) - finfo.fstar;
            
        elseif(isnan(finfo.fstar) && mu>0)
            % Value f_star not available, then we approximate the dual gap 
            % using fstar >= f(x) - (||f'(x)||^2) / (2*mu)
            tolf_handle = @(x) valf_handle(x);
            need_approx_fstar = true;
            approx_fstar_handle = @(x) valf_handle(x)  - (1/(2*mu))*norm(valfp_handle(x))^2;
            approx_fstar = approx_fstar_handle(finfo.x0) ;
            % we will use the best approx_fstar so far
        end
    end
    
    % Determine how to compute tolerance on the iterates
    if(isfield(tolfun,'tolx'))
        % tolx is defined by the user
        tolx_handle = tolfun.tolx;
    else
        % tolx is not defined by the user, so we use norm(x-xstar)
        tolx_handle = @(x) norm(x-finfo.xstar);
    end
end

nIteMax = nIteMax+1; % "iteration" 1 = x0

% Memory allocation
algo_x = zeros(finfo.n,nIteMax);
algo_x(:,1) = finfo.x0;

if(compute_tol)
    algo_tolf = zeros(1,nIteMax);
    algo_tolx = zeros(1,nIteMax);
    

    algo_tolf(1) = tolf_handle(finfo.x0);
    algo_tolx(1) = tolx_handle(finfo.x0);
else
    algo_tolf = NaN;
    algo_tolx = NaN;
end


% Main loop
for k=2:nIteMax
    [algo_x(:,k),algoparam] = algohandle(finfo,algo_x,k,algoparam);
    
    if(compute_tol)
        algo_tolx(:,k) = tolx_handle(algo_x(:,k));
        algo_tolf(k)  = tolf_handle(algo_x(:,k));
    end
    
    % Update the best value of fstar so far
    if(need_approx_fstar)
        approx_fstar = max(approx_fstar,approx_fstar_handle(algo_x(:,k)) );
    end
end

% if fstar is approximated, then we update algo_tolf with the best approx
% (i.e. the lowest value) of approx_fstar
if(need_approx_fstar)
    algo_tolf = algo_tolf-approx_fstar;
end
