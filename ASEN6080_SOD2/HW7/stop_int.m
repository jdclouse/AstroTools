function [value,isterminal,direction] = stop_int(t,state,xxx)
% Locate the time when height passes through zero in a decreasingdirection
% and stop integration.
value = norm(state(1:3)) - 3*925000; % detect y-1/2 = 0
isterminal = 1; % stop the integration
direction = -1; % negative direction