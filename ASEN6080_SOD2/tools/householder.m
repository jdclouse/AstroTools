function [ out_A ] = householder( A )
%householder Transformation for StatOD
%   This only householders out to cols-1!

[rows, cols] = size(A);
n = cols - 1;
for kk = 1:n
    sig = sign(A(kk,kk))*sqrt(sum(A(kk:end,kk).*A(kk:end,kk)));
    uk = A(kk,kk) + sig;
    A(kk,kk) = -sig;
    u = A(kk+1:end,kk);
    beta = 1/(sig*uk);
    for jj = (kk+1):(n+1)
        gamma = beta*(uk*A(kk,jj) + sum(u.*A(kk+1:end,jj)));
%         A(kk+1:end,jj) = A(kk+1:end,jj) - gamma*u;
        A(kk:end,jj) = A(kk:end,jj) - gamma*[uk; u];
    end
    A(kk+1:end,kk) = zeros(rows-(kk),1);
end

out_A = A;
end

