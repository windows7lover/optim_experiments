function  [f,fp,fpp,L,xstar] = logsumexp_diag(n)

d = rand(n,1);

dx = @(x) d.*x;

h = @(x) exp(dx(x)) + exp(-dx(x));
hp = @(x) d.*(exp(dx(x)) - exp(-dx(x)));
hpp = @(x) (d.^2).*h(x);

f = @(x) log(sum(h(x)));
fp = @(x) (1/sum(h(x))) * hp(x) ;
fpp = @(x) -hp(x)*hp(x)'/(sum(h(x))^2) + diag(hpp(x))/sum(h(x));

L = 2*norm(fpp(0)); % security
xstar = zeros(n,1);