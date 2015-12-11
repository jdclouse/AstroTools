close all
% clear all

G1 = BitShiftRegister(10, [3,10]);
getG1 = [];
for ii = 1:1023
    getG1 = [getG1 G1.update()];
end
ns = 10; nrep = 10;
y = getG1;
yy = repmat(y-mean(y),ns,nrep);
len = prod(size(yy)); 
fstep = 1.023e6*ns/len;
yy = reshape(yy,len,1); f = fstep*[0:1:len-1]';
ps = 20*log10(abs(fft(yy)));
figure
plot(f,ps,'.')
xlim([0 max(f)/2]), ylim([0 100])
title('G1 Code (Maximal)') 