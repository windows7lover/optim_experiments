function [maxval,p,normp,optval] = reg_minimax(sigma,lambda,degree,npoints,npointstest)



points = linspace(0,sigma,npoints);
pointstest = linspace(0,sigma,npointstest);
pointPower = vander(points);
pointPower = pointPower(:,end-degree:end);

cvx_begin quiet
    cvx_precision high
    variable c(degree+1,1)
    variable t
    
    minimize( t^2 + lambda*c'*c)
    
    subject to
        sum(c) == 1;
        abs(pointPower*c) <= ones(npoints,1)*t;
cvx_end

optval = cvx_optval;

p = c;
maxval = max(abs(polyval(p,pointstest)));
normp = norm(c);
