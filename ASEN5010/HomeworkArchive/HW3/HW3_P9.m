%% Problem 9: S&J 2.12
fprintf('Problem 9: S&J 2.12\n');
clearvars -except function_list pub_opt
close all

M = [1,1,2,2];
R(:,1) = [1;-1;2];
R(:,2) = [-1;-3;2];
R(:,3) = [2;-1;-1];
R(:,4) = [3;-1;-2];


Rdot(:,1) = [2;1;1];
Rdot(:,2) = [0;-1;1];
Rdot(:,3) = [3;2;-1];
Rdot(:,4) = [0;0;1];

Rc = [dot(R(1,:),M); dot(R(2,:),M); dot(R(3,:),M)]/sum(M);

Rcdot = [dot(Rdot(1,:),M); dot(Rdot(2,:),M); dot(Rdot(3,:),M)]/sum(M);

for i = 1:length(M)
    r(:,i) = R(:,i)-Rc;
    rdot(:,i) = Rdot(:,i)-Rcdot;
end

KE_trans = 0.5*sum(M)*dot(Rcdot, Rcdot);
KE_rot = 0;
for i = 1:length(M)
    KE_rot = KE_rot + 0.5*M(i)*dot(rdot(i), rdot(i));
end
fprintf('a)\n')
fprintf('Translational Kinetic Energy: %f mass units * velocity units\n',...
    KE_trans)
fprintf('Rotational Kinetic Energy: %f mass units * velocity units\n',...
    KE_rot)
fprintf('b)\n')

% eqn 2.65
origin_sigmas = R;
COM_sigmas = r;
origin_sigma_dot = Rdot;
COM_sigma_dot = rdot;

% The following fills out eqn 2.66. 
H_origin = 0;
H_COM = 0;
for i = 1:length(M)
    H_origin = H_origin + ...
        cross(origin_sigmas(:,i),M(i)*origin_sigma_dot(:,i));
    H_COM = H_COM + cross(COM_sigmas(:,i),M(i)*COM_sigma_dot(:,i));
end
fprintf('Angular Momentum Vector about origin:\n')
printVector(H_origin, 'mas*pos*vel');
fprintf('Angular Momentum Vector about COM:\n')
printVector(H_COM, 'mas*pos*vel');