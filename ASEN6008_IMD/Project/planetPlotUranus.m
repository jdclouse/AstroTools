d_vec = MGA_arr(1):1:JGA_arr(end);
num_dates = length(d_vec);
r_mars_store = zeros(3,num_dates);
r_jup_store = zeros(3,num_dates);
r_ura_store = zeros(3,num_dates);
jup_color = 2;
sat_color = 3;
nep_color = 1;
for ii = 1:num_dates
    [r_mars_store(:,ii), ~] = MeeusEphemeris(Mars, d_vec(ii), Sun);
    [r_jup_store(:,ii), ~] = MeeusEphemeris(Jupiter, d_vec(ii), Sun);
    [r_ura_store(:,ii), ~] = MeeusEphemeris(Uranus, d_vec(ii), Sun);
end

figure
hold on
plot(r_mars_store(1,:), r_mars_store(2,:), 'Color',color_order(sat_color,:))
plot(r_jup_store(1,:), r_jup_store(2,:), 'Color',color_order(jup_color,:))
plot(r_ura_store(1,:), r_ura_store(2,:), 'Color',color_order(nep_color,:))
axis equal