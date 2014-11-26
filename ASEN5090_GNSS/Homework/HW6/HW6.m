%% HW 6
% read_GPSbroadcast, broadcast2xv, adjust year, and read_rinex functions 
% provided from class are used in this homework.

%% Initialize
clearvars -except function_list pub_opt 
close all
c = 2.99792458e8; %m/s

%% Geometric Range
fprintf('1)\n')
fprintf('Broadcast file for 2014-09-12 is brdc2550.14n.\n\n')
eph = read_GPSbroadcast('brdc2550.14n');
%Week is col 19
%TOE is col 20

fprintf('2)\n')
GPS_Week = eph(1,19);
fprintf('GPS week is %d.\n', GPS_Week)
GPS_TOD = [1 03 00];
TOW = eph(1,20)+GPS_TOD(1)*3600 + GPS_TOD(2)*60 + GPS_TOD(3);
fprintf('GPS time of week is %d s.\n\n', TOW)

fprintf('3)\n')
[health,x,v,relcorr,satClkCorr] = broadcast2xv(eph,[GPS_Week TOW],21);
fprintf('PRN 21 position, in meters ECEF:\n')
fprintf('  %.2f\n',x(1))
fprintf('  %.2f\n',x(2))
fprintf('  %.2f\n',x(3))
fprintf('\n')
% disp(x)


fprintf('5)\n')
% User position is at Schriever AFB, home of the control segment
fprintf('User position is at Schriever AFB, home of the control segment\n\n')
userpos = [-1248596.2520 -4819428.2840 3976506.0340]'; %m

fprintf('6-7)\n')
prn_list = [21 22 26];
gps_vel = sqrt(3.986e5/26600);
range_store = zeros(length(prn_list),1);
clock_store = zeros(length(prn_list),1);
cntr = 1;
for ii = prn_list
[range0, range1] = compute_range(eph, ii, TOW, userpos);
range_diff = abs(range0 - range1);
fprintf('Uncorrected range for PRN %d:\n', ii)
fprintf('\t%.2f m\n', range0)
fprintf('TOF corrected range for PRN %d:\n', ii)
fprintf('\t%.2f m\n', range1)
fprintf('Difference in ranges for PRN %d:\n', ii)
fprintf('\t%.2f m\n', range_diff)
tof = range1/c;

fprintf('Approx time of flight is %.1e seconds.\n', tof)
fprintf('Approx in-track displacement of the satellite\n')
fprintf(' during that time is %.0f m.\n', gps_vel*1e3*tof)

range_store(cntr) = range1;
cntr = cntr+1;
fprintf('\n')
end
fprintf('All of the range differences make sense due to Tr-Tt\n')
fprintf(' and how far a satellite would move in that time.\n')
fprintf(' Satellites with longer range (like 22 and 26), have\n')
fprintf(' larger differences because the signal takes longer to travel.\n\n')

%% RINEX
%There are 7 observations in the file. P1 is the third one.
fprintf('8)\n')
fprintf('There are 7 observations in the file. P1 is the third one.\n\n');
fprintf('9)\n')
[ rinexv3 ] = read_rinex_obs('test.14o', prn_list);
fprintf('\n\n');

fprintf('10)\n')
C1_store = zeros(length(prn_list),1);
P1_store = zeros(length(prn_list),1);
TOW_measurements = rinexv3.data(rinexv3.data(:,2)==TOW,:);

for ii = 1:3
    [health,x,v,relcorr,satClkCorr] = broadcast2xv(eph,[GPS_Week TOW],prn_list(ii));
    C1_store(ii) = TOW_measurements(TOW_measurements(:,3) == prn_list(ii),8);
    P1_store(ii) = TOW_measurements(TOW_measurements(:,3) == prn_list(ii),6);
    clock_store(ii) = satClkCorr;
    fprintf('PRN %d C1 - R = %.2f m\n', ii, C1_store(ii) - range_store(ii))
    fprintf('PRN %d P1 - R = %.2f m\n', ii, P1_store(ii) - range_store(ii))
    fprintf('PRN %d clock error = %.2f m\n\n', ii, clock_store(ii))
end
    
fprintf('The code measurements vs the computed geometric range are\n')
fprintf(' pretty huge, up to 100s of km. But the error in the clocks\n')
fprintf(' appear to account for the discrepency, down to 10s of m.\n\n')