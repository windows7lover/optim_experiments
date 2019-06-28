function Hv = quad_bfgs(accumulator,H0,v)


[YC,RC] = accumulator.getC();
RCplus = (RC'*RC)\RC';

RCplusv = RCplus*v;
Pv = RC*(RCplusv);
IPv = v-Pv;

H1 = YC*(RCplusv);
H2 = RCplus'*(YC'*IPv);
H3 = H0*IPv - RC*(RCplus*IPv);

Hv = H1+H2+H3;