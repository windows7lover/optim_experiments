function [maxval,p,normp,optval] = const_minimax(sigma,tau,degree,npoints,npointstest)



points = linspace(0,sigma,npoints);
pointstest = linspace(0,sigma,npointstest);
pointPower = vander(points);
pointPower = pointPower(:,end-degree:end);

cvx_begin quiet
    cvx_precision high
    variable c(degree+1,1)
    variable t
    
    minimize( t )
    
    subject to
        sum(c) == 1;
        abs(pointPower*c) <= ones(npoints,1)*t;
        c'*c <= (1+tau)/(sqrt(degree+1));
cvx_end

optval = cvx_optval;

p = c;
maxval = max(abs(polyval(p,pointstest)));
normp = norm(c);
