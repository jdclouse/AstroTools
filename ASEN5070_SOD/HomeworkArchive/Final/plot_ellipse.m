function C=plot_ellipse(a,b,angle,x,y,color)
fcnPrintQueue(mfilename('fullpath')); % Add this code to code app 
if nargin > 5
    c = color;
else
    c = 'k';
end

tmp=(0:pi/100:2*pi)';
x_ell_t=a * cos(tmp);
y_ell_t=b * sin(tmp);

x_ell = x_ell_t * cos(angle) - y_ell_t * sin(angle) + x;
y_ell = x_ell_t * sin(angle) + y_ell_t * cos(angle) + y;
C=[cos(angle) -sin(angle); sin(angle) cos(angle)];

plot(x_ell,y_ell,c), axis equal