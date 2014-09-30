function plot_oe(oe_vec, times, prefix, P, vec2)
%plot_oe   Plot OEs
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

mu = 3.986004e5;
a = oe_vec(1,:);
e = oe_vec(2,:);
f = oe_vec(6,:);
if nargin < 4
    P = 2*pi/sqrt(mu/(a(1)*a(1)*a(1)));
end
elem_vec = oe_vec;
Tp = get_time_of_periapsis(a,e,f, times);
if nargin == 5
    elem_vec = vec2 - oe_vec;
    Tp = get_time_of_periapsis(vec2(1,:),vec2(2,:),vec2(6,:), times)- Tp;
end
fighandle = figure;
set(fighandle, 'Position', [100, 100, 600, 850])
subplot(6,1,1)
plot(times/3600, elem_vec(1,:))
period_lines(elem_vec(1,:), P, times)
ylabel(sprintf('%sa (km)', prefix))
xlabel('time (hours)')
subplot(6,1,2)
plot(times/3600, elem_vec(2,:))
period_lines(elem_vec(2,:), P, times)
ylabel(sprintf('%se', prefix))
xlabel('time (hours)')
subplot(6,1,3)
plot(times/3600, elem_vec(3,:)*180/pi)
period_lines(elem_vec(3,:)*180/pi, P, times)
ylabel(sprintf('%si (deg)', prefix))
xlabel('time (hours)')
subplot(6,1,4)
plot(times/3600, elem_vec(4,:)*180/pi)
period_lines(elem_vec(4,:)*180/pi, P, times)
ylabel(sprintf('%s\\Omega (deg)', prefix))
xlabel('time (hours)')
subplot(6,1,5)
plot(times/3600, elem_vec(5,:)*180/pi)
period_lines(elem_vec(5,:)*180/pi, P, times)
ylabel(sprintf('%s\\omega (deg)', prefix))
xlabel('time (hours)')
subplot(6,1,6)
plot(times/3600, Tp)
period_lines(Tp, P, times)
ylabel(sprintf('%sT_p (s)', prefix))
xlabel('time (hours)')

end

function period_lines(vec, P, times)
hold on
ymax = max(vec);
ymin = min(vec);
num_p = floor(times(end)/P);

for ii = 1:num_p
    [~, x] = min(abs(times - P*ii));
    plot([times(x)/3600, times(x)/3600], [ymax, ymin], 'r--')
end

hold off
end