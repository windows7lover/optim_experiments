function [Hv, H] = getSymHessian_compare(accumulator,H0,v)

% NOT OPTIMAL: Should return matrix-vector vector product instead


[YC,RC] = accumulator.getC();

[n,k] = size(YC);
k=k+1;
[U,S,V1] = svd(RC','econ');
I = eye(n);

s = diag(S);
Sinv = diag(1./s);
Sinv(isnan(Sinv)) = 0;
B1 = (U'*YC')*V1;
Z1 = zeros(k-1,k-1);
for i=1:k-1
    for j=1:k-1
        if((s(i)^2+s(j)^2) > 0)
            Z1(i,j) = (s(i)*B1(i,j) + s(j)*B1(j,i))/(s(i)^2+s(j)^2);
        end
    end
end

H1 = V1*Z1*V1';
H2 = (YC-V1*(V1'*YC))*(U*Sinv)*V1'; % Order of operations matters! Check latter
H2 = H2+H2';
H3 = (I-V1*V1')*H0*(I-V1*V1'); % Order of operations matters! Check latter

H = H3 + H1 + H2;
Hv = H*v;

end