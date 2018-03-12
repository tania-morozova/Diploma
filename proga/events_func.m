function [value,isterminal,direction] = events_func(t,x,P)
    value = x(1) - P(1);
    isterminal = 1;
    direction = 0;
end