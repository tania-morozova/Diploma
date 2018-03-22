function [value,isterminal,direction] = events_func(t,x,P,u2)
    P_min = P(u2(1));
    P_max = P(u2(2));
    
    value = [x(1) - P_min(1), (x(2) - P_min(2))*(x(2) - P_max(2))];
    isterminal = [1,0];
    direction = [0,-1];
end