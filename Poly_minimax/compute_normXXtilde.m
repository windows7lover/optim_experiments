function out = compute_normXXtilde(r,L,mu,M,normX0,k)

e = zeros(k,1);

r1 = 1-mu/L;
r2 = r^2;

for i=2:k
    e(i) = r1*e(i-1) + (M/(2*L))*(r2^i)*normX0^2;
end

out = sum(e);