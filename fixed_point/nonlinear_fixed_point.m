function X = nonlinear_fixed_point(f,x0,algo,k,param)
% Try top optimize f with 'algo'. Returns the results after k iterations 
% (i.e. X = [x0,x1,...,xk])
% We call at each steps algo(f,xHist,k,param), where at iteration "i" we 
% have xHist = [x0,x1,...,xi,0,0,...,0] where length(xHist) = k+1 at each
% iteration.


x = x0;
X = [x, zeros(length(x),k)];

for i=1:k
    x = algo(f,X,i+1,param);
    X(:,i+1) = x;
end