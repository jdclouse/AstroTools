function x = backsub( R,b,n )
%backsub Summary of this function goes here
%   Detailed explanation goes here

x = zeros(n,1);
x(n) = b(n)/R(n,n);
for ii = n-1:-1:1
    x(ii) = (b(ii) - sum(R(ii,ii+1:end).*x(ii+1:end)'))/R(ii,ii);
end

end

