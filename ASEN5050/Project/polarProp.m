function [ X_dot, alpha ] = polarProp( t,X,ctrl_law,beta )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

CelestialConstants; % import useful constants
% beta = 0.05;

r = X(1);
th = X(2);
rd = X(3);
thd = X(4);

alpha = pi/2; % No affect on 2-body acceleration

while th > 2*pi
    th = th - 2*pi;
end

if strcmp(ctrl_law,'HoldMax')
    alpha = 35*pi/180;
elseif strcmp(ctrl_law,'Zero')
    alpha = 0;
elseif strcmp(ctrl_law,'Optimal')
    alpha_calc = pi/2 - atan2(rd,r*thd);
    % The optimal angle.
    alpha = atan2(-3+sqrt(9+8*tan(alpha_calc)^2),4*tan(alpha_calc));
    if alpha > atan(1/sqrt(2))
        alpha = atan(1/sqrt(2));
    end
elseif strcmp(ctrl_law,'RaisePeri')
    r_i = Euler2DCM('3',X(2))*[X(1);0;0];
    v_i = Euler2DCM('3',X(2))*[X(3); ...
                            (X(1)*X(4));...
                             0];
    [a,e,~,~,~,f] = cart2OE(r_i,v_i,Sun.mu);
    p = a*(1-e^2);
    alpha_calc = atan2((2*e/(1+e)-1)*sin(f), 2*p/r/(1+e)-r/p*(1+cos(f))-1);
    % The optimal angle.
    alpha = atan2(-3+sqrt(9+8*tan(alpha_calc)^2),4*tan(alpha_calc));
    if alpha > atan(1/sqrt(2))
        alpha = atan(1/sqrt(2));
    end
elseif strcmp(ctrl_law,'Perpendicular')
    alpha = pi/2;
end

X_dot = [rd;
        thd;
        r*thd*thd + Sun.mu/r/r*(-1+beta*cos(alpha)*cos(alpha)*cos(alpha));
        (-2*rd*thd + beta*Sun.mu/r/r*cos(alpha)*cos(alpha)*sin(alpha))/r];

end

