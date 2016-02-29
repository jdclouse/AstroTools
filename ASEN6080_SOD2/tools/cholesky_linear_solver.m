function x = cholesky_linear_solver(Y,N)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

vlen = length(N);
z = zeros(vlen,1);
x = zeros(vlen,1);
U = chol(Y); % Upper
UT = U';
for ii = 1:vlen
    tmp = N(ii);
    for jj = 1:(ii-1)
        tmp = tmp - UT(ii,jj)*z(jj);
    end
    z(ii) = tmp/UT(ii,ii);
end

for ii = vlen:-1:1
    tmp = z(ii);
    for jj = (ii+1):vlen
        tmp = tmp - U(ii,jj)*x(jj);
    end
    x(ii) = tmp/U(ii,ii);
end
