function [xnew,c] = block_broyden(Y,X,G,lambda)

c = [];

dY = diff(Y,1,2);
dG = diff(X-Y,1,2);

H = dY*pinv(dG);

xnew = Y(:,end)-H*(X(:,end)-Y(:,end));