%% CelestialConstants
%% Description
% All sorts of constants for orbital mechanics purposes
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

%% Earth
Earth.name = 'Earth';
Earth.mu = 3.986004415e5; %km3/s2
Earth.R = 6378; %km
Earth.a = 149598023; %km
Earth.spin_rate = 7.2921158553e-05; %rad/s
Earth.flattening = 1/298.25722; %WGS-84
Earth.J2 = 0.0010826267;
Earth.P_days = 365.2421897; %days
Earth.P_years = 0.99997862; %days
Earth.m = 5.9742e24; %kg

%% Moon
Moon.name = 'Moon';
Moon.R = 1738.0; %km
Moon.J2 = 0.0002027;
Moon.P_days = 27.321582; %days
Moon.mu = 4902.799; %km3/s2
Moon.m = 7.3483e22; %kg
Moon.a = 384400; %km

%% Sun
Sun.mu = 1.32712428e11; %km3/s2
Sun.m = 1.9891e30; %kg

%% Mercury
Mercury.name = 'Mercury';
Mercury.R = 2439.0; %km
Mercury.J2 = 0.00006;
Mercury.P_days = 87.9666; %days
Mercury.mu = 2.2032e4; %km3/s2

%% Venus
Venus.name = 'Venus';
Venus.a = 108208601; %km
Venus.R = 6052.0; %km
Venus.J2 = 0.000027;
Venus.P_days = 224.6906; %days
Venus.mu = 3.257e5; %km3/s2

%% Mars
Mars.name = 'Mars';
Mars.a = 227939186; %km
Mars.R = 3397.2; %km
Mars.J2 = 0.001964;
Mars.P_days = 686.9150; %days
Mars.mu = 4.305e4; %km3/s2
Mars.m = 6.4191e23; %kg

%% Jupiter
Jupiter.name = 'Jupiter';
Jupiter.a = 778298361; %km
Jupiter.R = 71492; %km
Jupiter.J2 = 0.01475;
Jupiter.P_years = 11.856525; %days
Jupiter.P_days = Jupiter.P_years/Earth.P_years*Earth.P_days; %days
Jupiter.mu = 1.268e8; %km3/s2
Jupiter.m = 1.8988e27; %kg

%% Saturn
Saturn.name = 'Saturn';
Saturn.a = 1429394133; %km
Saturn.R = 60268; %km
Saturn.J2 = 0.01645;
Saturn.P_years = 29.423519; %days
Saturn.P_days = Saturn.P_years/Earth.P_years*Earth.P_days; %days
Saturn.mu = 3.794e7; %km3/s2
Saturn.m = 5.685e26; %kg

%% Uranus
Uranus.name = 'Uranus';
Uranus.R = 25559; %km
Uranus.J2 = 0.012;
Uranus.P_years = 83.747406; %days
Uranus.P_days = Uranus.P_years/Earth.P_years*Earth.P_days; %days
Uranus.mu = 5.794e6; %km3/s2

%% Neptune
Neptune.name = 'Neptune';
Neptune.R = 24764; %km
Neptune.J2 = 0.004;
Neptune.P_years = 163.7232045; %days
Neptune.P_days = Neptune.P_years/Earth.P_years*Earth.P_days; %days
Neptune.mu = 6.809e6; %km3/s2

%% Celestial units
au2km = 149597870.7;

%% Physical constants
day2sec = 86400; % sec/day
speed_of_light = 299792458; %m/s