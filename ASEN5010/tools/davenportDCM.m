function DCM_BI = davenportDCM( body_vecs, inrtl_vecs, weights )
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

[eigenvecs, eigenvals] = eigs(K);

DCM_BI = EulerParam2DCM(eigenvecs(:,1));
end

