function DCM_BI = QUESTDCM( body_vecs, inrtl_vecs, weights )
fcnPrintQueue(mfilename('fullpath'))
%Assumes input vectors are unitized!

% # columns = # measurements
[rows, cols] = size(body_vecs);

B = zeros(3);
for i = 1:cols
    B = B + weights(i)*(body_vecs(:,i)*inrtl_vecs(:,i)');
end

S = B + B';
sigma = trace(B);
Z = [B(2,3)-B(3,2);...
    B(3,1)-B(1,3);...
    B(1,2)-B(2,1)];

K = [sigma, Z'; Z, S-sigma*eye(3)];

% Newton's method to find max eigenvalue
% First guess = sum of the weights
eigenval = sum(weights);
f_eigen = det(K-eigenval*eye(4));
tol = 1e-15;
while f_eigen > tol
    eigenval = eigenval - f_eigen/det(eye(4));
    f_eigen = det(K-eigenval*eye(4));
end
    
CRP = inv((eigenval+sigma)*eye(3) - S)*Z;

DCM_BI = CRP2DCM(CRP);
end

