% test woodbury

d = 10;
N = 5;

W = rand(d);
W = W'*W;

Jref = rand(d);
Jref = (Jref'*Jref);

G = rand(d);
G = G'*G;
GmI = (G-eye(d));

Y = rand(d,N);
xstar = randn(d,1);
XSTAR = repmat(xstar,1,N);
R = GmI*(Y-XSTAR);

YC = diff(Y,1,2);
RC = diff(R,1,2);

J = (YC'*W*YC)\(YC'*W);
J = Jref + (RC-Jref*YC)*J;

Jrefinv = inv(Jref);
JrefinvRC = (Jref\RC);

Jinv = (YC'*W*Jrefinv*RC)\(YC'*W*Jrefinv);
Jinv = Jrefinv + (YC-JrefinvRC)*Jinv;

norm(inv(J) - Jinv)
[norm(Y - Jinv*R - XSTAR) norm(Y - R - XSTAR)]
