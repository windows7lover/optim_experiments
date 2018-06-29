function [Q, p, A, b] = transform_svm_dual(c,X,y)

[npoints,~] = size(X);
yX = diag(y)*X;

Q = yX*yX';
p = -ones(npoints,1);

A = [-eye(npoints) ; eye(npoints)];
b = [zeros(npoints,1) ; c*ones(npoints,1)];