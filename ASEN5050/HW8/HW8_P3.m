%% HW8 Problem 3
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

% Orbit normal vector
normal_vec = @(X) cross(X(1:3),X(4:6))/norm(cross(X(1:3),X(4:6)));

% Anonymous function to calculate 2-body accel
two_body = @(t,X) [X(4);X(5);X(6);...
    -Earth.mu*X(1)/norm(X(1:3))^3;...
    -Earth.mu*X(2)/norm(X(1:3))^3;...
    -Earth.mu*X(3)/norm(X(1:3))^3] + 1e-6*[0;0;0;normal_vec(X)];

hp = 400;
ha = 1000;
a = (ha+hp+2*Earth.R)/2;
e = (ha+Earth.R-a)/a;
i = 51.5*pi/180;
RAAN = 0;
w = 60*pi/180;
f = 0;
P = 2*pi*sqrt(a^3/Earth.mu);

[r,v] = OE2cart(a,e,i,RAAN,w,f,Earth.mu);
X0 = [r;v];

tol=1e-12;
options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);

[t_array,X_array]=ode45(two_body,[0 P],X0,options);

OE_array = zeros(6,length(t_array));
di_dt = zeros(length(t_array),1);
dRAAN_dt = zeros(length(t_array),1);
dw_dt = zeros(length(t_array),1);
for ii = 1:length(t_array)
    [OE_array(1,ii),...
        OE_array(2,ii),...
        OE_array(3,ii),...
        OE_array(4,ii),...
        OE_array(5,ii),...
        OE_array(6,ii)] = cart2OE(X_array(ii,1:3)',X_array(ii,4:6)',Earth.mu);
    
    r = norm(X_array(ii,1:3));
    h = norm(cross(X_array(ii,1:3),X_array(ii,4:6)));
    di_dt(ii) = r*cos(OE_array(5,ii)+OE_array(6,ii))/...
        (sqrt(Earth.mu/OE_array(1,ii)^3)...
        *OE_array(1,ii)^2*sqrt(1-OE_array(2,ii)^2))*1e-6;
    dRAAN_dt(ii) = r*sin(OE_array(5,ii)+OE_array(6,ii))/...
        (sqrt(Earth.mu/OE_array(1,ii)^3)...
        *OE_array(1,ii)^2*sqrt(1-OE_array(2,ii)^2)*sin(OE_array(3,ii)))...
        *1e-6;
    dw_dt(ii) = r*cot(OE_array(3,ii))*sin(OE_array(5,ii)+OE_array(6,ii))...
        /(h)*1e-6;
        
end
% plot(t_array,OE_array(3,:));
% figure
% plot(t_array,[OE_array(4,OE_array(4,:)<pi), OE_array(4,OE_array(4,:)>=pi)-2*pi]);
% figure
% plot(t_array,OE_array(5,:));

fprintf(['a) The inclination, RAAN, and argument of periapse are\n'...
         'directly affected by this force, according to Gaussian VOP.\n'])
fprintf(['b) The inclination, RAAN, and argument of periapse will\n'...
         'all experience secular drift because energy is constantly\n'...
         'added to the system.\n'])
fprintf(['c) When the orbit is exagerated, you can see evidence of\n'...
         'secular drift. The rate of change is biased toward either side\n'...
         'of zero over the course of an orbit.\n'])
figure('Position',[0 0 hw_pub.figWidth hw_pub.figHeight])
subplot(3,1,1)
plot(t_array,di_dt*180/pi*day2sec,'LineWidth',2)
ylabel('$\dot{i}$ (deg/day)','interpreter','latex')
subplot(3,1,2)
plot(t_array,dRAAN_dt*180/pi*day2sec,'LineWidth',2)
ylabel('$\dot{\Omega}$ (deg/day)','interpreter','latex')
subplot(3,1,3)
plot(t_array,dw_dt*180/pi*day2sec,'LineWidth',2)
ylabel('$\dot{\omega}$ (deg/day)','interpreter','latex')
xlabel('Time (sec)')

figure('Position',[0 0 hw_pub.figWidth hw_pub.figHeight])
subplot(3,1,1)
plot(OE_array(6,:),di_dt*180/pi*day2sec,'LineWidth',2)
ylabel('$\dot{i}$ (deg/day)','interpreter','latex')
subplot(3,1,2)
plot(OE_array(6,:),dRAAN_dt*180/pi*day2sec,'LineWidth',2)
ylabel('$\dot{\Omega}$ (deg/day)','interpreter','latex')
subplot(3,1,3)
plot(OE_array(6,:),dw_dt*180/pi*day2sec,'LineWidth',2)
ylabel('$\dot{\omega}$ (deg/day)','interpreter','latex')
xlabel('True Anomaly (deg)')