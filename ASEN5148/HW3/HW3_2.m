%% Team Homework
%% Initial conditions
% Constraints, assumptions, constants
m_payload = 2.7e-3; % kg (Wikipedia)
r_earth = 149598261; %km (Wikipedia)
r_mars = 227939100; %km (Wikipedia)
Re = 6378.; %km (Wikipedia)
mu = 3.986e5; % km3/s2 (Wikipedia)
Rm = 3396.2; % km (Wikipedia)
mu_mars = 4.2828e4; % km3/s2 (Wikipedia)
soi_mars = 5.77e5; % km, Laplace (Brown)
mu_sun = 1.32712440018e11; % km3/s2 (Wikipedia)
alt_final = 200; %km

% LV assumptions
% LV takes SC to circular parking orbit
lv_apogee_alt = 300; %km

% Isp
Isp = [333; % s, UDMH
       450; % s, LOX + LH
       250]; % s, shuttle SRB
prop_types = {'UDMH';
       'LOX + LH';
       'Shuttle SRB'};
gc = 9.80665; % m/s2

for ii = 1:3
Ve = Isp(ii)*gc;

fs = 0.1; % ratio of prop struct mass to prop

%% Calculate dV's
% All velocities in km/s
% Using patched conics

% Required Excess velocity
% Hohmann transfer
V_earth = sqrt(mu_sun/r_earth);
V_mars = sqrt(mu_sun/r_mars);
a_hoh = (r_earth + r_mars)/2;
Vp = sqrt(2*mu_sun/r_earth - mu_sun/a_hoh);
Va = sqrt(2*mu_sun/r_mars - mu_sun/a_hoh);
V_he = Vp - V_earth;


% Escape velocity from equator
dV_lo = sqrt(2*mu/Re - mu/(Re+lv_apogee_alt/2));
V_la = sqrt(2*mu/(Re+lv_apogee_alt) - mu/(Re+lv_apogee_alt/2));
V_leo = sqrt(mu/(Re+lv_apogee_alt));
dV_loc = V_leo - V_la;
dV_loss = 1.5; % aero loss
dV_leo = dV_loc + dV_loss + dV_lo;
V_pe = sqrt(2*mu/(Re+300) + (V_he^2));
dV_esc = V_pe - V_leo;

% Mars orbit injection
V_inf = abs(Va - V_mars);
a_hyp_mars = mu_mars/(V_inf*V_inf);
V_mars_peri = sqrt(2*mu_mars/(Rm+alt_final) + mu_mars/a_hyp_mars);
V_circ = sqrt(mu_mars/(Rm+alt_final));
dV_final = V_mars_peri - V_circ;

%% Calculate mass at important stages
% Work backwards from Mars
% Anon fcn to calculate the propellant mass
% mf = m_payload + fs*m_prop
calc_mp = @(m_payload, dV, Ve, fs) m_payload*(exp(dV/Ve)-1)/...
    (1-fs*(exp(dV/Ve)-1));

% Prop required for Mars orbit injection 
prop_at_mars = calc_mp(m_payload, dV_final*1e3, Ve, fs);
mass_mars_arrival = m_payload + prop_at_mars + fs*prop_at_mars;

% Prop required for Hohmann xfer
prop_at_x = calc_mp(mass_mars_arrival, dV_esc*1e3, Ve, fs);
mass_x = mass_mars_arrival + prop_at_x + fs*prop_at_x;

% Prop for single-stage LV to V_esc
prop_lv_ss = calc_mp(mass_x, dV_leo*1e3, Ve, fs);
mass_lv_ss = mass_x + prop_lv_ss + fs*prop_lv_ss;

% Prop for 2-stage LV to V_esc, equal dV
prop_lv_2s_2 = calc_mp(mass_x, dV_leo/2*1e3, Ve, fs);
mass_lv_2s_2 = mass_x + prop_lv_2s_2 + fs*prop_lv_2s_2;

prop_lv_2s_1 = calc_mp(mass_lv_2s_2, dV_leo/2*1e3, Ve, fs);
mass_lv_2s_1 = mass_lv_2s_2 + prop_lv_2s_1 + fs*prop_lv_2s_1;

% Prop for 3-stage LV to V_esc, equal dV
prop_lv_3s_3 = calc_mp(mass_x, dV_leo/3*1e3, Ve, fs);
mass_lv_3s_3 = mass_x + prop_lv_3s_3 + fs*prop_lv_3s_3;

prop_lv_3s_2 = calc_mp(mass_lv_3s_3, dV_leo/3*1e3, Ve, fs);
mass_lv_3s_2 = mass_lv_3s_3 + prop_lv_3s_2 + fs*prop_lv_3s_2;

prop_lv_3s_1 = calc_mp(mass_lv_3s_2, dV_leo/3*1e3, Ve, fs);
mass_lv_3s_1 = mass_lv_3s_2 + prop_lv_3s_1 + fs*prop_lv_3s_1;

fprintf(strcat(prop_types{ii}, ', Isp = ', num2str(Isp(ii)), '\n'));
if prop_lv_ss < 0
    fprintf('\tSingle-stage LV will not work for this propellant.\n')
    fprintf('\tNeed to increase Isp or decrease inert mass.\n\n')
else
    p_tot = prop_at_mars+prop_at_x+prop_lv_ss;
    fprintf('\tTotal prop: Single-stage LV: %f kg\n', p_tot)
    fprintf('\tTotal launch mass: Single-stage LV: %f kg\n\n', mass_lv_ss)
end

p_tot = prop_at_mars+prop_at_x+prop_lv_2s_2+prop_lv_2s_1;
fprintf('\tTotal prop: 2-stage LV: %f kg\n', p_tot)
fprintf('\tTotal launch mass: 2-stage LV: %f kg\n\n', mass_lv_2s_1)

p_tot = prop_at_mars+prop_at_x+prop_lv_3s_3+prop_lv_3s_2+prop_lv_3s_1;
fprintf('\tTotal prop: 3-stage LV: %f kg\n', p_tot)
fprintf('\tTotal launch mass: 3-stage LV: %f kg\n\n', mass_lv_3s_1)
end