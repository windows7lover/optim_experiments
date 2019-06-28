function xextr = rna( accum )

[Y,R] = accum.getYR();

k = size(Y,2);
onevec = ones(k,1);

c = (R'*R)\onevec;
c = c/sum(c);

xextr = (Y+R)*c;

end

