function Hv = quad_bfgs_precond(accumulator,H0,v)


[YC,RC] = accumulator.getC();
RCplus = (RC'*YC)\YC';

RCplusv = RCplus*v;
Pv = RC*(RCplusv);
IPv = v-Pv;

H1 = YC*(RCplusv);
H2 = RCplus'*(YC'*IPv);
H3 = H0*IPv - RC*(RCplus*H0*IPv);

Hv = H1+H2+H3;