function [p,maxval,normp] = reg_minimax(sigma,lambda,degree,npoints)



points = linspace(0,sigma,npoints);
pointPower = vander(points);
pointPower = pointPower(:,end-degree:end);

cvx_begin quiet
    variable c(degree+1,1)
    variable t
    
    minimize( t^2 + lambda*c'*c)
    
    subject to
        sum(c) == 1;
        abs(pointPower*c) <= ones(npoints,1)*t;
cvx_end

p = c;
maxval = t;
normp = norm(c);