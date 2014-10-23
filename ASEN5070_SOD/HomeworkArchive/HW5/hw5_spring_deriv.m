function state_dot = hw5_spring_deriv( times, state )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

w2 = (2.5+3.7)/1.5;
state_dot = zeros(6,1);
state_dot(1) = state(2);
state_dot(2) = -w2*state(1);

STM = reshape(state(3:end),2,2); %2x2 for one dimension
A = [0 1; -w2 0];
STM_dot = A*STM;

state_dot(3:6) = reshape(STM_dot,4,1);

end

