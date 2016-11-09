function A = gram_schmidt(B)

the_size = size(B)

A = zeros(the_size);
A(:,1) = B(:,1)/norm(B(:,1));

for ii = 2:the_size(2)
    e = B(:,ii);
    for jj = ii-1:1
        e = e-dot(A(:,jj),B(:,ii))*A(:,jj);
    end
    e_unit = e/norm(e);
    A(:,ii) = e_unit;
end