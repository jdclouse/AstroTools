theta0 = 0;
theta_dot = 7.29211585275553e-005;
JD_init = 2456296.25;
mu_sun = 132712440017.987; % km3/s2
solar_flux = 1357; %W/m2 @ 1 AU
A_m_ratio =  0.01; % m2/kg
Re = 6378.1363;

c = 299792.458; %km/s
mu_earth = 3.98600432896939e5; %km3/s2

EMO2EME_theta = 23.4393; %deg
EMO2EME = [1 0 0; 0 cosd(EMO2EME_theta) -sind(EMO2EME_theta); ...
    0 sind(EMO2EME_theta) cosd(EMO2EME_theta)];

au2km = 149597870.700; % km/AU

state_ap = [-274096790.0 ;
            -92859240.0 ;
            -40199490.0 ;
            32.67 ;
            -8.94 ;
            -3.88 ;
            1.2];

P = diag([100 100 100 0.1 0.1 0.1 0.1]);

filename = 'Interplanetary_Obs.txt';
delimiterIn = ',';
headerlinesIn = 2;
y = importdata(filename,delimiterIn,headerlinesIn);
Obs=y.data;

ObsMassaged = zeros(length(Obs), 5);
ObsMassaged(:,2) = Obs(:,1);
for ii = 1:length(Obs)
    for jj = 1:3
        if ~isnan(Obs(ii,jj+1))
            ObsMassaged(ii,1) = jj;
            ObsMassaged(ii,3) = Obs(ii,jj+1);
            ObsMassaged(ii,4) = Obs(ii,jj+4);
            break;
        end
    end
end

% i. Canberra Station (DSS 34) at latitude = -35.398333?
% , longitude = 148.981944?
% and altitude = 0.691750 km
% ii. Madrid Station (DSS 65) at latitude = 40.427222?
% , longitude = 4.250556?
% and altitude = 0.834539 km
% 1
% iii. Goldstone Station (DSS 13) at latitude = 35.247164?
% , longitude = 243.205?
% and altitude = 1.07114904 km

site(1).name = 'Canberra Station (DSS 34)';
site(1).id = 1;
site(1).lat_lon_alt = [-35.398333; 148.981944;0.691750]; % geodetic, deg, km
site(1).r = lla2ecef(site(1).lat_lon_alt', 0, Re); % m
site(2).name = 'Madrid Station (DSS 65)';
site(2).id = 2;
site(2).lat_lon_alt = [40.427222; 4.250556; 0.834539]; % geodetic, deg, km
site(2).r = lla2ecef(site(2).lat_lon_alt', 0, Re); % m
site(3).name = 'Goldstone Station (DSS 13)';
site(3).id = 3;
site(3).lat_lon_alt = [35.247164; 243.205; 1.07114904]; % geodetic, deg, km
site(3).r = lla2ecef(site(3).lat_lon_alt', 0, Re); % m