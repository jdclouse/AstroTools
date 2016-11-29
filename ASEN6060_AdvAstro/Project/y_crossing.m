function [value, isterminal, direction] = y_crossing(t,state,~)

value=state(2);
isterminal=1;
direction=1;
end