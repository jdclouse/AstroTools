%% HW 7
%% a)
clear
A = [-10 0 -10 0;
    0 -0.7 9 0;
    0 -1 -0.7 0;
    1 0 0 0];
B = [20 3; 0 0 ; 0 0; 0 0];

% Calculate the eigenvalues:
[V, e_vals] = eig(A);
fprintf('Eigenvalues:\n')
disp(diag(e_vals))

% Not Hurwitz. But the eigenvalues are all distinct, so we can calculate
% B_tilde and figure out if it's reachable. 

B_tilde = inv(V)*B
% Not all rows are reachable, since there is a pole of zero. The 
% unreachable modes are the complex conjugate pair. However, the 
% unreachable modes have real values on the LHP. Therefore, A is 
% stabilizable. 

%% b)
clear
A = [-5 1 0; 0 -5 0; 0 0 -10];
B = [0;1;1];

% Calculate the eigenvalues:
e_vals = eig(A);
fprintf('Eigenvalues:\n')
disp((e_vals))

% A is Hurwitz. Compute the Grammian to see if its rank matches rank(A)
sys = ss(A,B,zeros(0,length(A)),0);
G = gram(sys,'c')
fprintf('rank(A): ')
disp(rank(A))
fprintf('rank(G): ')
disp(rank(G))
% The system is reachable. For the controllability effort, get the min
% eigen value of G:
effort = min(eig(G));
fprintf('Reachability effort is %f.\n', 1/effort)

%% c)
A = [-3 1 0 0;0 -3 0 0;0 0 -2 1; 1 0 0 -2];
B = [0;0.001;0;1];

% Calculate the eigenvalues:
[V, e_vals] = eig(A);
fprintf('Eigenvalues:\n')
disp(diag(e_vals))

% A is Hurwitz. Compute the Grammian to see if its rank matches rank(A)
sys = ss(A,B,zeros(0,length(A)),0);
G = gram(sys,'c')
fprintf('rank(A): ')
disp(rank(A))
fprintf('rank(G): ')
disp(rank(G))
% The system is reachable. For the controllability effort, get the min
% eigen value of G:
effort = min(eig(G));
fprintf('Reachability effort is %f.\n', 1/effort)

%% d)
A = [5 -1 -3; 0 5.5 0; 0 0 -6];
B = [1;0;5];

% Calculate the eigenvalues:
[V, e_vals] = eig(A);
fprintf('Eigenvalues:\n')
disp(diag(e_vals))

% Not Hurwitz. But the eigenvalues are all distinct, so we can calculate
% B_tilde and figure out if it's reachable. 
B_tilde = inv(V)*B

% The system is not reachable, the second mode does not have a non-zero
% value for its control input. It's not stabilizable, because that
% eigenvalue has a positive real part. 

%% Conclusion
% For the reachable systems, the system in b) took much less effort to
% control than that of part c).  You can tell that the Grammian in part c)
% was close to being singular in the first two columns, compared to part
% b)'s Grammian.