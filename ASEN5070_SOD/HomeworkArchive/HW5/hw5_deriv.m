function state_dot = hw5_deriv(times, state)


x = state(1);
y = state(2);
r = sqrt(x*x + y*y);

state_dot = [state(3); state(4); -x/(r*r*r); -y/(r*r*r)];