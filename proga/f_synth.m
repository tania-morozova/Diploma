function dx = f_synth(t,x,u_min,u_max,f,P_curr)
if x(1) >= P_curr(1)
    u_curr = u_min;
else
    u_curr = u_max;
end
dx=f(t,x,u_curr);
end