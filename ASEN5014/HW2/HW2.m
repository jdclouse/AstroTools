%% HW 2 , Problem 2b and 2c

A = [0 1;-25 -4];
B = [1 1;0 1];
C = eye(2);
D = zeros(2,2);
[n,d] = ss2tf(A,B,C,D,1)
[n,d] = ss2tf(A,B,C,D,2)
figure('OuterPosition', [0 50 hw_pub.figWidth hw_pub.figHeight])
impulse(ss(A,B,C,D),5)