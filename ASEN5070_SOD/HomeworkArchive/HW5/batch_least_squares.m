function x_hat = batch_least_squares( x, y, W, H, x_bar, W_bar )
%batch_least_squares Prototype of BLS code.
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

STM = eye(length(x));
lam = W_bar;
N=lam*x_bar;

lam = lam + H'*W*H;
N = N + H'*W*y;

x_hat = lam\N;
end

