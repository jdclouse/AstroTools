%% Problem 12: S&J 4.3
fprintf('Problem 12: S&J 4.3\n');
clearvars -except function_list pub_opt
close all

b1=[0;1;0];
b2=[0;0;1];
b3=[1;0;0];

C_NB = [b1'; b2'; b3'];

I = [15 0 0 ; 0 11 5; 0 5 16];
[CT, lambda] = eigs(I);
if abs(det(CT)+1) < 1e-10
    CT(:,3) = cross(CT(:,1),CT(:,2));
end
C_FB = CT';
C_FB*I*C_FB' ;

f1_in_f = [1;0;0];
f2_in_f = [0;1;0];
f3_in_f = [0;0;1];

f1_in_n = C_NB*C_FB'*f1_in_f;
f2_in_n = C_NB*C_FB'*f2_in_f;
f3_in_n = C_NB*C_FB'*f3_in_f;

fprintf('a) Principle Inertias:\n');
printVector(max(lambda),'');

fprintf('b) Rotation Matrix FB:\n');
C_FB

fprintf('c) Principle body axes expressed in N-frame:\n');
fprintf('f1:\n');
printVector(f1_in_n,'');
fprintf('f2:\n');
printVector(f2_in_n,'');
fprintf('f3:\n');
printVector(f3_in_n,'');
