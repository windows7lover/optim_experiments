function out = rand_posdef(n)

A = rand(n);
[Q,~] = qr(A);
randvec = rand(n,1);
sigma = diag(randvec);
out = Q'*sigma*Q;
out = out/norm(out);