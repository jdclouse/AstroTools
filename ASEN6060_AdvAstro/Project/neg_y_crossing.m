function [value, isterminal, direction] = neg_y_crossing(t,state,~)

value=state(2);
isterminal=1;
direction=-1;
end