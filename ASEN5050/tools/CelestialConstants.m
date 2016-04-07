%% CelestialConstants
%% Description
% All sorts of constants for orbital mechanics purposes
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 


%% Celestial units
au2km = 149597870.7;

%% Physical constants
day2sec = 86400; % sec/day
speed_of_light = 299792458; %m/s

%% Earth
Earth.name = 'Earth';
Earth.mu = 3.986004415e5; %km3/s2
Earth.R = 6378; %km
Earth.a = 149598023; %km
Earth.spin_rate = 7.2921158553e-05; %rad/s
Earth.flattening = 1/298.25722; %WGS-84
Earth.oblate_ecc = 0.081819221456; %WGS-84
Earth.J2 = 0.0010826267;
Earth.P_days = 365.2421897; %days
Earth.P_years = 0.99997862; %days
Earth.m = 5.9742e24; %kg
% Meeus ephemeris parameters
Earth.Meeus.J200.L = [100.466449 35999.3728519 -0.00000568 0.0]; %deg
Earth.Meeus.J200.a = 1.000001018*au2km; %km
Earth.Meeus.J200.e = [0.01670862 -0.000042037 -0.0000001236 0.00000000004];
Earth.Meeus.J200.i = [0 0.0130546 -0.00000931 -0.000000034]; % deg
Earth.Meeus.J200.RAAN = [174.873174 -0.2410908 0.00004067 -0.000001327]; %deg
Earth.Meeus.J200.Pi = [102.937348 0.3225557 0.00015026 0.000000478]; %deg

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
Venus.m = 4.869e24; %km
Venus.Meeus.J200.L = [181.979801 58517.8156760 0.00000165 -0.000000002];%deg
Venus.Meeus.J200.a = 0.72332982*au2km; % km
Venus.Meeus.J200.e = [0.00677188 -0.000047766 0.0000000975 0.00000000044];
Venus.Meeus.J200.i = [3.394662 -0.0008568 -0.00003244 0.000000010];%deg
Venus.Meeus.J200.RAAN = [76.679920 -0.2780080 -0.00014256 -0.000000198];%deg
Venus.Meeus.J200.Pi = [131.563707 0.0048646 -0.00138232 -0.000005332];%deg

%% Mars
Mars.name = 'Mars';
Mars.a = 227939186; %km
Mars.R = 3397.2; %km
Mars.J2 = 0.001964;
Mars.P_days = 686.9150; %days
Mars.mu = 4.305e4; %km3/s2
Mars.m = 6.4191e23; %kg

% Meeus ephemeris parameters
Mars.Meeus.J200.L = [355.433275 19140.2993313 0.00000261 -0.000000003]; %deg
Mars.Meeus.J200.a = 1.523679342*au2km; %km
Mars.Meeus.J200.e = [0.09340062 0.000090483 -0.0000000806 -0.00000000035];
Mars.Meeus.J200.i = [1.849726 -0.0081479 -0.00002255 -0.000000027]; %deg
Mars.Meeus.J200.RAAN = [49.558093 -0.2949846 -0.00063993 -0.000002143]; %deg
Mars.Meeus.J200.Pi = [336.060234 0.4438898 -0.00017321 0.000000300]; %deg

%% Jupiter
Jupiter.name = 'Jupiter';
Jupiter.a = 778298361; %km
Jupiter.R = 71492; %km
Jupiter.J2 = 0.01475;
Jupiter.P_years = 11.856525; %days
Jupiter.P_days = Jupiter.P_years/Earth.P_years*Earth.P_days; %days
Jupiter.mu = 1.268e8; %km3/s2
Jupiter.m = 1.8988e27; %kg
Jupiter.Meeus.J200.L = [34.351484 3034.9056746 -0.00008501 0.000000004 ];
Jupiter.Meeus.J200.a = [5.202603191 0.0000001913 ]*au2km;
Jupiter.Meeus.J200.e = [0.04849485 0.000163244 -0.0000004719 -0.00000000197 ];
Jupiter.Meeus.J200.i = [1.303270 -0.0019872 0.00003318 0.000000092 ];
Jupiter.Meeus.J200.RAAN = [100.464441 0.1766828 0.00090387 -0.000007032 ];
Jupiter.Meeus.J200.Pi = [14.331309 0.2155525 0.00072252 -0.000004590 ];


%% Saturn
Saturn.name = 'Saturn';
Saturn.a = 1429394133; %km
Saturn.R = 60268; %km
Saturn.J2 = 0.01645;
Saturn.P_years = 29.423519; %days
Saturn.P_days = Saturn.P_years/Earth.P_years*Earth.P_days; %days
Saturn.mu = 3.794e7; %km3/s2
Saturn.m = 5.685e26; %kg
Saturn.Meeus.J200.L = [50.077471 1222.1137943 0.00021004 -0.000000019 ];
Saturn.Meeus.J200.a = [9.554909596 -0.0000021389 ]*au2km;
Saturn.Meeus.J200.e = [0.05550862 -0.000346818 -0.0000006456 0.00000000338 ];
Saturn.Meeus.J200.i = [2.488878 -0.0025515 -0.00004903 0.000000018 ];
Saturn.Meeus.J200.RAAN = [113.665524 -0.2566649 -0.00018345 0.000000357 ];
Saturn.Meeus.J200.Pi = [93.056787 0.5665496 0.00052809 0.000004882 ];

%% Uranus
Uranus.name = 'Uranus';
Uranus.R = 25559; %km
Uranus.J2 = 0.012;
Uranus.P_years = 83.747406; %days
Uranus.P_days = Uranus.P_years/Earth.P_years*Earth.P_days; %days
Uranus.mu = 5.794e6; %km3/s2
Uranus.Meeus.J200.L = [314.055005 428.4669983 -.00000486 0.000000006 ];
Uranus.Meeus.J200.a = [19.218446062 -0.0000000372 0.00000000098 ]*au2km;
Uranus.Meeus.J200.e = [0.04629590 -0.000027337 0.0000000790 0.00000000025 ];
Uranus.Meeus.J200.i = [0.773196 -0.0016869 0.00000349 0.000000016 ];
Uranus.Meeus.J200.RAAN = [74.005947 0.0741461 0.00040540 0.000000104 ];
Uranus.Meeus.J200.Pi = [173.005159 0.0893206 -0.00009470 0.000000413 ];

%% Neptune
Neptune.name = 'Neptune';
Neptune.R = 24764; %km
Neptune.J2 = 0.004;
Neptune.P_years = 163.7232045; %days
Neptune.P_days = Neptune.P_years/Earth.P_years*Earth.P_days; %days
Neptune.mu = 6.809e6; %km3/s2
Neptune.Meeus.J200.L = [304.348665 218.4862002 0.00000059 -0.000000002 ];
Neptune.Meeus.J200.a = [30.110386869 -0.000000166 0.00000000069 ]*au2km;
Neptune.Meeus.J200.e = [0.00898809 0.000006408 -0.0000000008 0 ];
Neptune.Meeus.J200.i = [1.769952 0.0002257 0.00000023 -0 ];
Neptune.Meeus.J200.RAAN = [131.784057 -0.0061651 -0.00000219 -0.000000078 ];
Neptune.Meeus.J200.Pi = [48.123691 0.0291587 0.00007051 -0 ];

%% Pluto
Pluto.name = 'Pluto';
Pluto.Meeus.J200.L = [238.92903833 145.20780515 0 0]; %deg
Pluto.Meeus.J200.a = [39.48211675 -0.00031596 0 0]*au2km; %km
Pluto.Meeus.J200.e = [0.24882730 0.00005170 0 0];
Pluto.Meeus.J200.i = [17.14001206 0.00004818 0 0]; %deg
Pluto.Meeus.J200.RAAN = [110.30393684 -0.01183482 0 0];%deg
Pluto.Meeus.J200.Pi = [224.06891629 -0.04062942 0 0];%deg
