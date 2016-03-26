function [value, isterminal, direction] = y_crossing(t,state,~)

value=0-state(2);
isterminal=1;
direction=0;
end