clear
clc
%%
syms x1 x2 x3 x4 r1 r2 r3 r4 b1 b2 b3 c2 c3 c4 u1 u2

dx1 = x1*(r1+u1-b1*x2);
dx2 = x2*(-r2 - b2*x3 + c2*x1);
dx3 = x3*(-r3 + u2 - b3*x4 + c3*x2);
dx4 = x4*(-r4 + c4*x3);

P1 = (r2*c4 + b2*r4)/(c4*c2);
P2 = (r1 + u1)/b1;
P3 = r4/c4;
P4 = (c3*P2 - r3 + u2)/b3;

%%
syms xi1 xi2 xi3 xi4 u_01 u_02 v1 v2

xi = [xi1; xi2; xi3; xi4]
x = [x1 x2 x3 x4].'
dx = [dx1 dx2 dx3 dx4].'
P = [P1 P2 P3 P4].'
u = [u1 u2].'
v = [v1 v2].'
u0 = [u_01 u_02].'

%%

%u = u0 + v
dxi = subs(dx,x,xi + P);
dxi = subs(subs(dxi,u1,u_01+v1),u2,u_02 + v2)

%%
% z1 = xi2
% z2 = dxi(2)
% z3 = xi4
% z4 = dxi(4)

syms z1 z2 z3 z4

S = solve([z2 == dxi(2), z4 == dxi(4)],[xi1,xi3],'ReturnConditions',true)
%%

% xi_1 = subs(subs([simplify(S.xi1);z2;simplify(S.xi3);z4],xi2,z1),xi4,z3)
% 
% dxi_1 = simplify(subs(subs(subs(subs(subs(subs(dxi,xi1,xi_1(1)),xi2,xi_1(2)),xi3,xi_1(3)),xi4,xi_1(4)),z2,dxi(2)),z4, dxi(4)))
