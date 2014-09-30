function C=drawellipse(a,b,th,xc,yc,lt)
%function C=drawellipse(a,b,th,xc,yc,lt)
% 
% Input: a=semimajor axis length, b=semiminor axis length
%        th=angle between semimajor axis and positive x axis.
%        (xc,yc) are optional inputs specifying the ellipse center
%        lt is an optional input specifying the line type
% Output: rotation matrix
% Author: P. Axelrad, University of Colorado 2/94

E=0:pi/100:2*pi;
x=a*cos(E);
y=b*sin(E);

C=[cos(th) -sin(th); sin(th) cos(th)];

z=C*[x;y];

if (nargin == 4)
   ltt=xc;
elseif (nargin == 6)
   ltt=lt;
else
   ltt='g';
end

if (nargin > 4)
   z=[z(1,:)+xc; z(2,:)+yc];
end

plot(z(1,:),z(2,:),ltt)
axis('equal')
