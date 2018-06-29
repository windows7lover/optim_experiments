function [p,maxval,normp,optval] = reg_minimax_U(sigma,lambda,degree,npoints1,npoints2)

if(nargin<5)
    npoints2 = npoints1;
end

points = linspace(0,sigma,npoints1);
pointsTest = [linspace(0,sigma,npoints2), points];
pointPower = vander(points);
pointPower = pointPower(:,end-(degree+1):end); % p(x)

pointsPoly = pointPower(:,1:end-1)-pointPower(:,2:end); % (1-x)p(x)

cvx_begin quiet
    cvx_precision best
    variable c(degree+1,1)
    variable t
    
    minimize( t^2 + lambda*(c'*c))
    
    subject to
        sum(c) == 1;
        abs(pointsPoly*c) <= ones(npoints1,1)*t;
cvx_end

p = c;
maxval = max(abs(polyval(p,pointsTest)));
maxval2 = max(abs(polyval([0 ; p]-[p ; 0],pointsTest))); % maxval of (1-x)*p
optval = maxval2^2+lambda*(p'*p);
normp = norm(p);