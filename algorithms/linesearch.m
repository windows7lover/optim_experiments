function alpha = linesearch(alphaUB,x,gradf,tol)

lb = 0;
ub = alphaUB;

d = -gradf(x);

gradAlpha = @(alpha) gradf(x+alpha*d)'*d;

if(gradAlpha(ub)<0) % too close
    alpha = ub;
    return
end

while((ub-lb)/alphaUB > tol)
    trial = (ub+lb)/2;
    if(gradAlpha(trial)>0) % too far
        ub = trial;
    else
        lb = trial;
    end
end
alpha = (ub+lb)/2;
