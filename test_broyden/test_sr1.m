% test B

clearvars
close all

d = 15;
N = 5;

A = rand(d);
A = A'*A+1e-1*norm(A);
A = A/norm(A);

B0 = eye(d);
B = B0;
x0 = rand(d,1);
xstar = rand(d,1);


g = @(x) -A*(x-xstar);

x=x0;
R = zeros(d,N);
X = zeros(d,N);
Y = zeros(d,N);

error_x = zeros(1,N);
error_x_multi = zeros(1,N);

for i=1:N
    Y(:,i) = x;
    if(i == 1)
        grad = g(x);
    else
        oldgrad = grad;
        grad = g(x);
        y = grad-oldgrad;
        s = x-oldx;
        
        % SR1
%         res = s-B*y;
%         B = B + res*res'/(res'*y);

        % BFGS
%         rho = y'*s;
%         I = eye(d);
%         B = (I-s*y'/rho)*B*(I-y*s'/rho) + s*s'/(y'*s);

%         % Broyden Type-II
%         B = B + (s-B*y)*y'/(y'*y);
        
        
        % Broyden Type-I
        B = B + (s-B*y)*s'*B/(s'*B*y);
    end
    oldx = x;
    x = x-B*grad;
    X(:,i) = x;
    R(:,i)=grad;
    
    error_x(i) = norm(g(x));
end



C = diff(eye(N),1,2);
RC = R*C;
YC=Y*C;
XC=X*C;

norm(B*RC-YC)
norm(B*RC-XC)
norm(-B*A-eye(d))/norm(A)


x=x0;
R = zeros(d,N);
X = zeros(d,N);
Y = zeros(d,N);
Bmulti = B0;
for i=1:N
    Y(:,i) = x;
    if(i > 1)
        C = diff(eye(i),1,2);
        RC = R(:,1:i)*C;
        YC=Y(:,1:i)*C;
        XC=X(:,1:i)*C;
        
        B0 = Bmulti;
        
        % Broyden Type-I
        Bmulti = (YC'*B0*RC)\(YC'*B0);
        Bmulti = B0 + (YC-B0*RC)*(Bmulti);
    end
    grad = g(x);
    
    R(:,i) = grad;
    
%     gamma = rand(i,1);
    gamma = (Y(:,1:i)'*B0*R(:,1:i))\ones(i,1);
    gamma = gamma/sum(gamma);
    
%     x = (Y(:,1:i)-Bmulti*R(:,1:i))*gamma;
    x = (Y(:,1:i)-B0*R(:,1:i))*gamma;
    X(:,i) = x;
    
    error_x_multi(i) = norm(g(x));
end

figure
semilogy(error_x)
hold on
semilogy(error_x_multi)

