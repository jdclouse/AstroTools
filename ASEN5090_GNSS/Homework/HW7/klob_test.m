% Display what the issue is with Klobuchar and higher altitudes. 
lats = [0 10 20 25 40];
colors = {'b', 'r', 'g', 'k', 'c'};
times = 1:1:24*3600;
I = zeros(1,length(times));
cnt = 1;
figure
hold on
for lat = lats
    for t = times
        I(t) = ...
        klobuchar( 0, pi/2, ...
        lat*pi/180, nist.lon*(pi/180), ...
        t, solar_min.a, solar_min.b );
    end
    adj_t = times/3600-(105/180*12);
    I = [I(adj_t>=0) I(adj_t<0)];
    adj_t = [adj_t(adj_t>=0) adj_t(adj_t<0)+24];
    plot(adj_t, I, colors{cnt})
    cnt = cnt+1;
end

figure
test = [];
cnt = 1;
is = 0:.01:.5;
for i = is
    test(cnt) = solar_min.a(1) + solar_min.a(2)*i ...
        + solar_min.a(3)*i*i + solar_min.a(4)*i*i*i;
    cnt = cnt + 1;
end
plot(is*180, test)