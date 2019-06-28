function [ f_seq ] = adapaccel( x0, finfo, kmax )

f = finfo.f;
d = finfo.n;
fstar = finfo.fstar;

mk = f(x0);
Ak = 0;
Gk = zeros(d,1);
lambdak = 1;
Lk = finfo.L;

f_seq = zeros(1,kmax+1);
f_seq(1,1) = f(x0) - fstar;

xk = x0;
vk = x0;

for i = 2:kmax+1
    % lambdak_plus = lambdak/2; 
    lambdak_plus = lambdak;
    
    ratio = (lambdak/lambdak_plus);
    zk = ratio * vk + (1-ratio)*x0;
    polynomial = [ -1/(2*lambdak_plus), 1/Lk, Ak/Lk ];
    alphak_plus = max(roots(polynomial));
    Ak_plus = Ak+alphak_plus;
    wk_plus = (Ak*xk + alphak_plus*zk)/Ak_plus;
    
    [xk_plus, Lk_plus, fxk_plus, gxk_plus] = damped_gradient(wk_plus, finfo, Lk);
    
    mk_plus = mk + alphak_plus*(fxk_plus + gxk_plus'*(x0-xk_plus));
    Gk_plus = Gk + alphak_plus*gxk_plus;
    vk_plus = x0-Gk_plus/lambdak_plus;
    
    if( Ak_plus*fxk_plus > mk_plus - norm(Gk_plus)^2/(2*lambdak_plus))
        
        display('WTF????')
        
        polynomial = [ -1/(2*lambdak), 1/Lk, Ak/Lk ];
        % roots(polynomial)
        alphak = max(roots(polynomial));
        wk = (Ak*xk + alphak*vk)/(Ak+alphak);
        Ak = Ak+alphak;

        [xk, Lk, fxk, gxk] = damped_gradient(wk, finfo, Lk);

        mk = mk + alphak*(fxk + gxk'*(x0-xk));
        Gk = Gk + alphak*gxk;
        vk = x0-Gk/lambdak;
        
        if( Ak*fxk > mk - norm(Gk)^2/(2*lambdak))
            display('lfpwerl WTF????')
        end
        
    else
        Ak = Ak_plus;
        mk = mk_plus;
        Gk = Gk_plus;
        vk = vk_plus;
        
        xk = xk_plus;
        Lk = Lk_plus;
        lambdak = lambdak_plus;
    end
    lambdak
    f_seq(1,i) = f(xk)-fstar;
end

end

