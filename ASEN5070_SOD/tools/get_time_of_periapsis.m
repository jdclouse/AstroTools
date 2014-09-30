function Tp = get_time_of_periapsis(a,e,f,t)
%get_time_of_periapsis   Return time of periapsis
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

mu = 3.986004e5;
% a = oe_vec(1,:);
% e = oe_vec(2,:);
% f = oe_vec(6,:);
mean_motion = sqrt(mu./(a.*a.*a));
E = get_eccentric_anomaly(f,e);
Tp = t - get_mean_anomaly(E, e)./mean_motion;
% plot(Tp)
% plot(mod(Tp, (2*pi)./mean_motion))
time_periapse_angle = mod(t, (2*pi)./mean_motion).*mean_motion ...
    - (E - e.*sin(E));

time_periapse_angle = unwrap(time_periapse_angle);
% plot(time_periapse_angle)

Tp = time_periapse_angle./mean_motion;

% plot(times/3600, Tp)