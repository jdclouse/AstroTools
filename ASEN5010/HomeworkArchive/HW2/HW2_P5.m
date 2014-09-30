%% Problem 5: S&J, Problem 3.13
fprintf('Problem 5: S&J, Problem 3.13\n')
clearvars -except function_list pub_opt
fprintf('Initial Euler Angles (3-1-3):\n')
euler_angles = [-30; 40; 20]; % degrees
printVector(euler_angles, 'degrees')

PRV = Euler2PRV('313', euler_angles * pi/180); % PHI * e
phi = norm(PRV);
e = PRV/phi;
fprintf('a) Principle rotation axis: ')%\n[%f\n %f\n %f]\n', e1, e2, e3)
printVector(e, '');
fprintf('b) Principle rotation angle Phi = %f degrees\n', phi * 180/pi)
fprintf('   Principle rotation angle Phi Prime = %f degrees\n', ...
    (phi - 2*pi) * 180/pi)
EP = Euler2EP('313', euler_angles * pi/180);
fprintf('c) Euler parameters (quaternions): \n[%f\n %f\n %f\n %f]\n', ...
    EP(1), EP(2), EP(3), EP(4))
CRP = zeros(3,1);
CRP(1) = EP(2)/EP(1);
CRP(2) = EP(3)/EP(1);
CRP(3) = EP(4)/EP(1);
fprintf('d) Classical Rodriques parameters: \n[%f\n %f\n %f]\n', ...
    CRP(1), CRP(2), CRP(3))
MRP = zeros(3,1);
MRP(1) = EP(2)/(1+EP(1));
MRP(2) = EP(3)/(1+EP(1));
MRP(3) = EP(4)/(1+EP(1));
fprintf('e) Modified Rodriques parameters: \n[%f\n %f\n %f]\n', ...
    MRP(1), MRP(2), MRP(3))