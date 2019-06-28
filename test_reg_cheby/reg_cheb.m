function poly = reg_cheb(N,alpha)

if(N == 0)
    poly = [1];
    return
end

if(N == 1)
    poly = (alpha)*[1 0] + (1-alpha)*[1 1]/2;
    return
end

tm2 = [0 1];
tm1 = [1 0];

for i=2:N
    tm = alpha*(2*[tm1 0] - [0 tm2]) + (1-alpha)*ones(1,i+1)/(i+1);
    tm2 = [0 tm1];
    tm1 = tm;
end
poly = tm;