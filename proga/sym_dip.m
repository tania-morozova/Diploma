syms x1 x2 x3 x4 r1 r2 r3 r4 b1 b2 b3 c2 c3 c4 u1 u2

dx1 = x1*(r1+u1-b1*x2);
dx2 = x2*(-r2 - b2*x3 + c2*x1);
dx3 = x3*(-r3 + u2 - b3*x4 + c3*x2);
dx4 = x4*(-r4 + c4*x3);

P1 = (r2*c4 + b2*r4)/(c4*c2);
P2 = (r1 + u1)/b1;
P3 = r4/c4;
P4 = (c3*P2 - r3 + u2)/b3;

K = x1 - P1*log(x1) + ...
    b1/c2*(x2 - P2*log(x2)) + ...
    (b1*b2)/(c2*c3)*(x3 - P3*log(x3)) + ...
    (b1*b2*b3)/(c2*c3*c4)*(x4 - P4*log(x4));
K1 = subs(subs(subs(subs(K,x1,P1),x2,P2),x3,P3),x4,P4);
KK = subs(subs(K,x1,P1),x2,P2)
%%
K2 = (K - K1);
syms u10 u20
K0=subs(subs(K,u2,u20),u1,u10);
dK0 = simplify(diff(K0,x1)*dx1+diff(K0,x2)*dx2+diff(K0,x3)*dx3+diff(K0,x4)*dx4);
d2K0 = simplify(diff(dK0,x1)*dx1+diff(dK0,x2)*dx2+diff(dK0,x3)*dx3+diff(dK0,x4)*dx4)