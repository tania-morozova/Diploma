function u = f_u(t,x,u_min,u_max,P,r,b,c,eps)
    P_min = P(u_min);
    P_max = P(u_max);
    if (abs(x(1) - P_min(1)) < eps)&&(x(2) < P_max(2))&&(x(2) > P_min(2))
        u_curr(1) = b(1)*x(2) - r(1);
    elseif x(1) >= P_min(1)
        u_curr(1) = u_min(1);
    else
        u_curr(1) = u_max(1);
    end
    
    if (abs(x(3) - P_min(3)) < eps)&&(x(4) < (c(3)*x(2) - r(3) + u_max(2))/b(3))...
            &&(x(4) > (c(3)*x(2) - r(3) + u_min(2))/b(3))
        u_curr(2) = r(3) + b(3)*x(4) - c(3)*x(2);
    elseif x(3) <= P_min(3)
        u_curr(2) = u_max(2);
    else
        u_curr(2) = u_min(2);
    end
    
    u = u_curr;
end