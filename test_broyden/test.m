%% Init

clear all
close all

d = 20;
N = 10;

G = rand(d,d);
G = G'*G;
I = eye(d);

Y = rand(d,N);
xstar = 10*rand(d,1);
XSTAR = repmat(xstar,1,N);
R = (G-I)*(Y-XSTAR);
C = eye(N,N-1) - [ zeros(1,N-1); eye(N-1,N-1)];
RC = R*C;
YC = Y*C;

beta = 0.5;
Href = beta*I;
Jref = (1/beta)*I;
% W = 0.1*eye(d);
W = inv(G-I);
Winv = inv(W);
Whalf = sqrtm(W);
Winvhalf = sqrtm(Winv);
gamma = zeros(N,1); gamma(end) = 1;


%% Check bad broyden formula

cvx_begin quiet
    variable H(d,d)
    
    minimize(  norm( Winvhalf*(H-Href)*Winvhalf ,'fro') )
    subject to
        H*RC == YC
    
cvx_end
Hcvx = H;
Hformula = Href + (YC-Href*RC)*((RC'*W*RC)\RC')*W;

% Check good broyden formula

cvx_begin quiet
    variable J(d,d)
    
    minimize(  norm( Whalf*(J-Jref)*Whalf ,'fro') )
    subject to
        RC == J*YC
    
cvx_end

Jcvx = J;
Jformula = Jref + (RC-Jref*YC)*((YC'*Winv*YC)\YC')*Winv;
Jinvformula = Href + (YC-Href*RC)*((YC'*Winv*Href*RC)\YC')*Winv*Href;


[norm(Hformula-Hcvx,'fro'), norm(Jformula-Jcvx,'fro'), norm(Jinvformula - inv(Jformula),'fro')]
[norm((Y-(Hformula-Jinvformula)*R)*gamma)]


%% Test anderson

onevec = ones(N,1);
cW = (R'*W*R)\onevec;
cW = cW/sum(cW);
cY = (Y'*R)\onevec;
cY = cY/sum(cY);
norm(cW-cY)

