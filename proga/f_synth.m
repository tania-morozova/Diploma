function dx = f_synth(t,x,u_min,u_max,f,P_curr)
    u_curr = [0,0];
    if x(1) >= P_curr(1)
        u_curr(1) = u_min(1);
    else
        u_curr(1) = u_max(1);
    end
    if x(3) >= P_curr(3)
        u_curr(2) = u_min(2);
    else
        u_curr(2) = u_max(2);
    end
    dx=f(t,x,u_curr);
end