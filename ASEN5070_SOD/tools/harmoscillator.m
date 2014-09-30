function dx = harmoscillator( t, x, kmratio )
%harmoscillator  Return output for 
fcnPrintQueue(mfilename('fullpath')) % Add this code to code appendix

dx = zeros(2,1);
if length(x) ~= 2
    return
end

dx(1) = x(2);
dx(2) = -kmratio*x(1);