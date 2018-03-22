function dx = f_synth(t,x,u_min,u_max,f,P,r,b)
    eps = 10^(-6);
    P_min = P(u_min);
    P_max = P(u_max);

    if (abs(x(1) - P_min(1)) < eps)&&(x(2) < P_max(2))&&(x(2) > P_min(2))
        alpha = (r(1) + u_max - b(1)*x(2))/(u_max - u_min);
        u_curr = alpha*u_min + (1-alpha)*u_max;
    elseif x(1) >= P_min(1)
        u_curr = u_min;
    else
        u_curr = u_max;
    end
    dx=f(t,x,u_curr);
end