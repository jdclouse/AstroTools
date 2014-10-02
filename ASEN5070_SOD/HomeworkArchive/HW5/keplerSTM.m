function Phi = keplerSTM(x0,dt,mu)
% keplerSTM - Generate the state transition matrix for keplerian orbits.
%   Phi = keplerSTM(x0,dt,mu) will return the state transition matrix Phi
%   for a set of objects described by initial positions and velocities (x0)
%   with gravitational parameters mu, propogated along Keplerian orbits for
%   time dt.
% 
%   x0 must be a row vector of length 6n where n is the number of planets. 
%   For each planet, x0 must include the 3 position and 3 velocity values
%   in the order [r1,r2,r3,v1,v2,v3].' These positions and velocities are 
%   measured in a cartesian coordinate system with the central object 
%   (i.e. star) at the origin.  For multiple planets, simply stack
%   several of these vectors on top of each other.  mu must contain n
%   values equal to G(m+ms) where G is the gravitational constant, m is the
%   mass of the orbiting object, and ms is the mass of the central object.
%   dt is a scalar value of time to propogate. 
%   NOTE: All units should be consistent.  If the positions are in AU and
%   velocities are in AU/day, then dt must be in days, and mu must be in
%   AU^3/day^2.
% 
%   This is calculated using the algorithm described in Shepperd, 1984
%   which employs Goodyear's universal variables and solves the Kepler
%   problem using continued fractions.
% 
%   Example:
%       %propogate a planet 10 days along its orbit:
%       x0 = [-0.8684, -0.6637, 0.7207, 0.0039, 0.0056, 0.0120].';
%       x1 = keplerSTM(x0,10,2.6935e-04)*x0;
% 
%   Written by Dmitry Savransky, 22 April 2008
%   17 November 2008 - Fixed several bugs associated with calculations for
%   more than one planet

%determine number of planets and validate input
nplanets = length(x0)/6;
if nplanets - floor(nplanets) > 0
    disp('The length of the input vector must be a multiple of 6.');
    Phi = -1;
    return;
end
if length(mu) ~= nplanets
    disp('The mu vector must contain the same number of values as the length of the input vector divided by 6.');
    Phi = -2;
    return;
end

%orient mu properly
mu = mu(:).';

%create position and velocity matrices
x0 = reshape(x0,6,nplanets);
r0 = x0(1:3,:);
v0 = x0(4:6,:);

%constants
r0norm = sqrt(sum(r0.^2));
nu0 = dot(r0,v0);
beta = 2*mu./r0norm - dot(v0,v0);

%initialization
u = zeros(1,nplanets);

%For elliptic orbits, calculate period effects
deltaU = zeros(1,length(beta));
eorbs = beta > 0;
if sum(eorbs) > 0
    P = 2*pi*mu(eorbs).*beta(eorbs).^(-3/2);
    n = floor((dt + P/2 - 2*nu0(eorbs)./beta(eorbs))./P);
    deltaU(eorbs) = 2*pi*n.*beta(eorbs).^(-5/2);
end

%continued fraction constants
a = 5;
b = 0;
c = 5/2;
k = 1 - 2*(a-b);
l = 2*(c-1);
d = 4*c*(c-1);
n = 4*b*(c-a);

%kepler iteration loop
t = zeros(1,nplanets);
counter = 0;
%loop until convergence of the time array to the time step
while sum(abs((t-dt)/dt*100)) > 1e-3 && counter < 1000;
    q = beta.*u.^2./(1+beta.*u.^2);
    U0w2 = 1 - 2*q;
    U1w2 = 2*(1-q).*u;
    U = 16/15*U1w2.^5 .* contFrac(q) + deltaU;
    U0 = 2*U0w2.^2-1;
    U1 = 2*U0w2.*U1w2;
    U2 = 2*U1w2.^2;
    U3 = beta.*U + U1.*U2/3;
    r = r0norm.*U0 + nu0.*U1 + mu.*U2;
    t = r0norm.*U1 + nu0.*U2 + mu.*U3;
    u = u - (t-dt)./(4*(1-q).*r);
    counter = counter+1;
end
if counter == 1000
    disp('Failed to converge on t');
end

%kepler solution
f = 1 - mu./r0norm.*U2;
g = r0norm.*U1 + nu0.*U2;
F = -mu.*U1./r./r0norm;
G = 1 - mu./r.*U2;

%construct the state transition matrix
Phi = zeros(6*nplanets);
for j=1:nplanets
    st = (j-1)*6+1;
    Phi(st:st+5,st:st+5) = [eye(3)*f(j)  eye(3)*g(j);eye(3)*F(j) eye(3)*G(j)];
end

    %calculate continued fraction
    function G = contFrac(x)
        
        %initialize arrays
        A = zeros(1,length(x))+1;
        B = zeros(1,length(x))+1;
        G = zeros(1,length(x))+1;

        Gprev = zeros(1,length(x))+2;
        counter = 0;
        %loop until convergence of continued fraction
        while sum(abs((G-Gprev)./G*100)) > 1e-3 && counter < 1000
            k = -k;
            l = l+2;
            d = d+4*l;
            n = n+(1+k)*l;
            A = d./(d - n*A.*x);
            B = (A-1).*B;
            Gprev = G;
            G = G + B;
            counter = counter+1;
        end
        if counter == 1000
            disp('Failed to converge on G');
        end
    end
end