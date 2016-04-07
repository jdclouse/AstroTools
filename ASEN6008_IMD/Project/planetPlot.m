d_vec = JGA_arr(525):1:SGA_arr(1452);
num_dates = length(d_vec);
r_sat_store = zeros(3,num_dates);
r_jup_store = zeros(3,num_dates);
r_nep_store = zeros(3,num_dates);
jup_color = 2;
sat_color = 3;
nep_color = 1;
for ii = 1:num_dates
    [r_sat_store(:,ii), ~] = MeeusEphemeris(Saturn, d_vec(ii), Sun);
    [r_jup_store(:,ii), ~] = MeeusEphemeris(Jupiter, d_vec(ii), Sun);
    [r_nep_store(:,ii), ~] = MeeusEphemeris(Neptune, d_vec(ii), Sun);
end

figure
hold on
plot(r_sat_store(1,:), r_sat_store(2,:), 'Color',color_order(sat_color,:))
plot(r_jup_store(1,:), r_jup_store(2,:), 'Color',color_order(jup_color,:))
plot(r_nep_store(1,:), r_nep_store(2,:), 'Color',color_order(nep_color,:))
axis equal