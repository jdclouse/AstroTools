function state_dot = hw5_deriv(times, state)

num_states = 4;
x = state(1);
y = state(2);
r = sqrt(x*x + y*y);
r3 = r*r*r;
state_dot = [state(3); state(4); -x/(r3); -y/(r3)];

if length(state) > num_states %Need to integrate the STM
    STM = reshape(state(num_states+1:end),num_states,num_states);

    r5=r3*r*r;
    A = [
        0, 0, 1, 0;
        0, 0, 0, 1;
        (-1/r3+3*x*x/r5), (3*x*y/r5), 0, 0;
        (3*x*y/r5), (-1/r3+3*y*y/r5), 0, 0];

    STM_dot = A*STM;
    state_dot = [state_dot; reshape(STM_dot, num_states*num_states, 1)];
end