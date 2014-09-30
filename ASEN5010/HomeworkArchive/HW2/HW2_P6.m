%% Problem 6: S&J, Problem 3.14
fprintf('Problem 6: S&J, Problem 3.14\n')
clearvars -except function_list pub_opt

phi = 45 * pi/180; %radians
e = [1;1;1;]/sqrt(3); %unit vectorPlot
PRV=phi*e;
euler_ypr = DCM2Euler('321', PRV2DCM(PRV)) * 180/pi; %degrees
fprintf('Euler 3-2-1 angles: ')
printVector(euler_ypr, 'deg')