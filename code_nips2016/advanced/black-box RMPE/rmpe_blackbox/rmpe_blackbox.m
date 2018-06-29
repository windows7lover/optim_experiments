function xout = rmpe_blackbox(funval,iterations,param)
% x_extrapolated = RMPE_BLACKBOX(funval,iterations,param)
%
% Apply the RMPE algorithm on the given iterations to extrapolate a new 
% point.
% The regularization parameter is automatically found using a grid search.
%
%   List of inputs
%   --------------
%
% 	- funval: a function handle, returning the value of the objective 
%     function
% 	- iterations: Matrix of size d x (k+1) , where d is the dimention of 
%     the space and k+1 the number of points to extrapolate
% 	- param: a structure which contains the options of the acceleration 
%     algorithm
%
%   Description of the parameters
%   -----------------------------
%
% * param.doAdaptiveLambda = true; 
% 		Determine if lambda should change over time
%
% * param.doLineSearch = true; 
% 		Determine if we should perform a line search on the stepsize at the
%       end
%
% * param.lambda = 1; 
%   param.lambdamin = 1e-10;
%       >   (if param.doAdaptiveLambda == true)
%       These two values determine the range of the grid search 
%       [lambda,lambdamin]
%       >   (if param.doAdaptiveLambda == false)
%		The value param.lambda fixes the value of the 
%       regularization.
%
% * param.forceDecrease = true; 
%       Optionnal, check if the extrapolated value is smaller than the last
%       iterate of gradient method


%% Inputs check

if(~isfield(param,'doLineSearch'))
    param.doLineSearch = true;
end

k = size(iterations,2)-1;

if(~isfield(param,'lambda'))
    param.lambda = 1;
end

if(~isfield(param,'lambdamin'))
    param.lambdamin = 1e-10;
end

lambda = param.lambda; % we use lambda as "lambda max"
lb = param.lambdamin; % "lambda min"

if(~isfield(param,'doAdaptiveLambda'))
    param.doAdaptiveLambda = true;
end

if(~isfield(param,'forceDecrease'))
    param.forceDecrease = true;
end

if(param.doAdaptiveLambda == true)
    lambdavec = [0, logspace(log10(lb),log10(lambda),k)];
else
    lambdavec = lambda;
end

%% Initializations

fvalvec = zeros(size(lambdavec)); % for the grid search
algo_x = iterations; % sequence
UU = diff(algo_x,1,2);
UU = UU'*UU;
UU_norm = UU/norm(UU); % normalized matrix U

%% Grid search on lambda

for i=1:length(lambdavec) % length(lambdavec) = 1 if not adaptive
    xout = rmpe(algo_x, UU_norm, lambdavec(i)); 
    % extrapolation using differents values of lambda
    
    if(param.doAdaptiveLambda)
        fvalvec(i) = funval(xout);
    else
        fvalvec = 0;
    end
end

if( param.forceDecrease && min(fvalvec) > funval(algo_x(:,end)) ) % Force decrease
    xout = algo_x(:,end); % If bad extrapolation, return the last point of the algorithm
else
    [~, idx_min] = min(fvalvec);
    lambdamin = lambdavec(idx_min);
    xout = rmpe(algo_x, UU_norm, lambdamin);
end

% Line search on the stepsize
if(param.doLineSearch)
    % Find a good stepsize, i.e a good alpha such that
    %           f( x_0 - alpha*(step) )
    % is small, where step = (x_rmpe-x_0).
    x = iterations(:,1);
    step = xout-x;
    sizestep = 1;
    fold = funval(x+sizestep*step);
    sizestep = 2*sizestep;

    while(true)
        fnew = funval(x+sizestep*step);
        if(fold>fnew)
            fold = fnew;
            sizestep = 2*sizestep;
        else
            break;
        end
    end
    sizestep = sizestep/2;
    xout = x+sizestep*step;
end
