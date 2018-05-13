function [value,isterminal,direction] = events_func3(t,x,c,b,r,P,u_min,u_max,eps)
    P_min = P(u_min);
    P_max = P(u_max);
    d_min = (c(3)*x(2) + u_min(2) - r(3))/(b(3));
    d_max = (c(3)*x(2) + u_max(2) - r(3))/(b(3));
    if (abs(x(1) - P_min(1)) <= eps)&&(x(2) >= P_min(2) - eps)&&(x(2) < P_max(2) + eps)&&(abs(x(3) - P_min(3)) <= eps)&&(x(4) > d_min - eps)&&(x(4) < d_max + eps)
        isterminal = 1;
        value = 0;
    else
        value = 1;
        isterminal = 0;
    end
    direction = 0;
end