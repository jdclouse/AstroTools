function [value, isterminal, direction] = neg_y_crossing(t,state,~)

value=norm([state(1) state(2)])-1;
isterminal=1;
direction=-1;
end