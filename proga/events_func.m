function [value,isterminal,direction] = events_func2(t,x,P)
    value = [x(1) - P(1), x(3) - P(3)];
    isterminal = [1,1];
    direction = [0,0];
end