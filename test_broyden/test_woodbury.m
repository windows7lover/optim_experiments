% test woodbury

d = 10;
N = 5;

I = eye(d);


beta = rand(1);
J0 = rand(d);
% J0 = (1/beta)*(J0'*J0);
J0 = (1/beta)*I;
J0inv = inv(J0);

G = rand(d);
G = G'*G;
G = 0.99*G/norm(G);
GmI = (G-I);

% W = rand(d);
% W = W'*W;
W = (GmI);

Y = rand(d,N);
xstar = randn(d,1);
XSTAR = repmat(xstar,1,N);
R = GmI*(Y-XSTAR);

YC = diff(Y,1,2);
RC = diff(R,1,2);

YCWP = (YC'*W*YC)\(YC'*W);
PW = YC*YCWP;

J = (RC*YCWP) + (RC*YCWP)'*(I-PW) + (I-PW)'*J0*(I-PW);
T = (RC*YCWP)' + (I-PW)'*J0;
Tinv = (RC'*J0inv*W*YC)\(YC-J0inv*RC)';
Tinv = J0inv + J0inv*W*(YC)*Tinv;
Jinv = (YC'*W*Tinv*RC)\(YC'*W*Tinv);
Jinv = Tinv + (YC-Tinv*RC)*Jinv;

[norm(T*Tinv-I), norm(J*Jinv-I), norm(inv(J)-Jinv), norm(Jinv-Jinv'), norm(Jinv*RC-YC)]

temp1 = (RC'*RC)\RC';
temp2 = (YC'*RC)\YC';
Jinvsimple = beta*(I-RC*temp1) + YC*temp2;
norm(Jinvsimple - Jinv)

% SR1 test

H0 = -1000*J0inv;
sr1 = ((YC-H0*RC)'*RC)\(YC-H0*RC)';
sr1 = H0 + (YC-H0*RC)*sr1;

o = ones(N,1);
gRand = rand(N);
gRand = gRand/sum(gRand);
mat1 = R'*(inv(GmI)-H0)*R;
mat2 = Y'*R-R'*H0*R;
g1 = mat1\o; g1 = g1/sum(g1);
g2 = mat2\o; g2 = g2/sum(g2);
[ (Y-sr1*R)*gRand, (Y-H0*R)*g1,  (Y-H0*R)*g2 ]
