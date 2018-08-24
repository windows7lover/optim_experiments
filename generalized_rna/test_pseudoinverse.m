addpath(genpath('../cvx'))

d = 20;


Href = rand(d);
Href = Href+Href';

% W = eye(d);

W = rand(d);
W = W*W' + 0.1*eye(d);
Whalf = sqrtm(W);


% V = W;
V = rand(d);
V = V*V' + 0.1*eye(d);

k = 4;
S = rand(d,k)';
% S = [S;S]
% Y = rand(d,k)';
Y = S*V;


cvx_begin
    variable H(d,d) symmetric
    minimize norm( Whalf* (H-Href) *Whalf, 'fro' )
    subject to
        Y*H==S
cvx_end
H_cvx = H;

% Api = pinv(A');
I = eye(d);

Ypi = (Y*inv(W))'*pinv(Y*inv(W)*Y');
P = I-Ypi*Y;
H_pi = Ypi*S + P*(Ypi*S)' + P*Href*P';

% norm(H_pi,'fro') - norm(Ypi*S,'fro') -norm(P*(Ypi*S)','fro')-norm(P*Href*P','fro')
% norm(H_pi,'fro') - norm(Ypi*S+P*(Ypi*S)'+P*Href*P','fro')

[norm( Whalf*(H_pi-Href)*Whalf ,'fro') norm(H_pi-H_cvx,'fro')]

%%
A = Y;
B = S-A*Href;

Aplus = pinv(A);
mat1 = Aplus*B + (I-Aplus*A)*(Aplus*B)';
mat2 = B'*pinv(A*B')*B;

s1 = svd(mat1);
s2 = svd(mat2);
[s1 s2]'



