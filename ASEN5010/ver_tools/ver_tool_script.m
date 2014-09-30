%% Init
addpath('C:\Users\John\Documents\ASEN5010\AIAA_Software\Matlab Toolbox\DOS-Windows\')

%% Euler angles to DCM
clear all
euler_set = rand(3,1) * pi;
tol = 1e-10;
fprintf('Euler Angles to DCM:\n')
myDCM = Euler2DCM('123', euler_set);
trueDCM = Euler1232C(euler_set);
matrixCompare('123', trueDCM - myDCM, tol)

myDCM = Euler2DCM('121', euler_set);
trueDCM = Euler1212C(euler_set);
matrixCompare('121', trueDCM - myDCM, tol)

myDCM = Euler2DCM('132', euler_set);
trueDCM = Euler1322C(euler_set);
matrixCompare('132', trueDCM - myDCM, tol)

myDCM = Euler2DCM('131', euler_set);
trueDCM = Euler1312C(euler_set);
matrixCompare('131', trueDCM - myDCM, tol)

myDCM = Euler2DCM('213', euler_set);
trueDCM = Euler2132C(euler_set);
matrixCompare('213', trueDCM - myDCM, tol)

myDCM = Euler2DCM('212', euler_set);
trueDCM = Euler2122C(euler_set);
matrixCompare('212', trueDCM - myDCM, tol)

myDCM = Euler2DCM('231', euler_set);
trueDCM = Euler2312C(euler_set);
matrixCompare('231', trueDCM - myDCM, tol)

myDCM = Euler2DCM('232', euler_set);
trueDCM = Euler2322C(euler_set);
matrixCompare('232', trueDCM - myDCM, tol)

myDCM = Euler2DCM('312', euler_set);
trueDCM = Euler3122C(euler_set);
matrixCompare('312', trueDCM - myDCM, tol)

myDCM = Euler2DCM('313', euler_set);
trueDCM = Euler3132C(euler_set);
matrixCompare('313', trueDCM - myDCM, tol)

myDCM = Euler2DCM('321', euler_set);
trueDCM = Euler3212C(euler_set);
matrixCompare('321', trueDCM - myDCM, tol)

myDCM = Euler2DCM('323', euler_set);
trueDCM = Euler3232C(euler_set);
matrixCompare('323', trueDCM - myDCM, tol)

%% Euler angles to inv(B)
clear all
euler_set = rand(3,1) * pi;
tol = 1e-10;
fprintf('Euler Angles to inv(B):\n')

myBinv = BinvEuler('321', euler_set);
trueBinv = BinvEuler321(euler_set);
matrixCompare('321', trueBinv - myBinv, tol)

%% Euler angles to B
clear all
euler_set = rand(3,1) * pi;
tol = 1e-10;
fprintf('Euler Angles to B:\n')

myBmat = BmatEuler('321', euler_set);
trueBmat = BmatEuler321(euler_set);
matrixCompare('321', trueBmat - myBmat, tol)

%% Euler angle derivative
clear all
euler_set = rand(3,1) * pi;
w_body_set = rand(3,1) * 5 * pi/180; % up to 5 deg/s
tol = 1e-10;
fprintf('Euler Angles and body rates to Euler Angle derivative:\n')

mydE = dEuler('321', euler_set, w_body_set);
truedE = dEuler321(euler_set, w_body_set);
matrixCompare('321', truedE - mydE, tol)

%% Euler angles to PRV
clear all
euler_set = rand(3,1) * pi;
tol = 1e-10;
fprintf('Euler Angles to principle rotation vector:\n')
fprintf('Euler Angles to Euler Params:\n')
fprintf('Euler angles: ')
printVector(euler_set*180/pi, 'deg')

myPRV = Euler2PRV('121', euler_set);
truePRV = Euler1212PRV(euler_set);
matrixCompare('121', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('123', euler_set);
truePRV = Euler1232PRV(euler_set);
matrixCompare('123', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('131', euler_set);
truePRV = Euler1312PRV(euler_set);
matrixCompare('131', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('132', euler_set);
truePRV = Euler1322PRV(euler_set);
matrixCompare('132', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('212', euler_set);
truePRV = Euler2122PRV(euler_set);
matrixCompare('212', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('213', euler_set);
truePRV = Euler2132PRV(euler_set);
matrixCompare('213', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('231', euler_set);
truePRV = Euler2312PRV(euler_set);
matrixCompare('231', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('231', euler_set);
truePRV = Euler2312PRV(euler_set);
matrixCompare('231', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('312', euler_set);
truePRV = Euler3122PRV(euler_set);
matrixCompare('312', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('313', euler_set);
truePRV = Euler3132PRV(euler_set);
matrixCompare('313', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('321', euler_set);
truePRV = Euler3212PRV(euler_set);
matrixCompare('321', abs(truePRV) - abs(myPRV), tol)

myPRV = Euler2PRV('323', euler_set);
truePRV = Euler3232PRV(euler_set);
matrixCompare('323', abs(truePRV) - abs(myPRV), tol)

%% Euler angles to Euler parameters
clear all
euler_set = rand(3,1) * pi;
tol = 1e-10;
fprintf('Euler Angles to Euler Params:\n')
fprintf('Euler angles: ')
printVector(euler_set*180/pi, 'deg')

myEP = Euler2EP('121', euler_set);
trueEP = Euler1212EP(euler_set);
matrixCompare('121', trueEP - myEP, tol)

myEP = Euler2EP('123', euler_set);
trueEP = Euler1232EP(euler_set);
matrixCompare('123', trueEP - myEP, tol)

myEP = Euler2EP('131', euler_set);
trueEP = Euler1312EP(euler_set);
matrixCompare('131', trueEP - myEP, tol)

myEP = Euler2EP('132', euler_set);
trueEP = Euler1322EP(euler_set);
matrixCompare('132', trueEP - myEP, tol)

myEP = Euler2EP('212', euler_set);
trueEP = Euler2122EP(euler_set);
matrixCompare('212', trueEP - myEP, tol)

myEP = Euler2EP('213', euler_set);
trueEP = Euler2132EP(euler_set);
matrixCompare('213', trueEP - myEP, tol)

myEP = Euler2EP('231', euler_set);
trueEP = Euler2312EP(euler_set);
matrixCompare('231', trueEP - myEP, tol)

myEP = Euler2EP('232', euler_set);
trueEP = Euler2322EP(euler_set);
matrixCompare('232', trueEP - myEP, tol)

myEP = Euler2EP('313', euler_set);
trueEP = Euler3132EP(euler_set);
matrixCompare('313', trueEP - myEP, tol)

myEP = Euler2EP('312', euler_set);
trueEP = Euler3122EP(euler_set);
matrixCompare('312', trueEP - myEP, tol)

myEP = Euler2EP('321', euler_set);
trueEP = Euler3212EP(euler_set);
matrixCompare('321', trueEP - myEP, tol)

myEP = Euler2EP('323', euler_set);
trueEP = Euler3232EP(euler_set);
matrixCompare('323', trueEP - myEP, tol)

%% DCM to EA
clear all
euler_set = rand(3,1) * pi;
tol = 1e-10;
fprintf('DCM to Euler Angles:\n')
fprintf('Euler angles: ')
printVector(euler_set*180/pi, 'deg')
% seq_string = {'121','123','131','132','213','212','231','232','312','313','321','323'};
seq_string = '321';
DCM = Euler2DCM(seq_string, euler_set);
euler_angles = DCM2Euler( seq_string, DCM );
truth_EA = C2Euler321(DCM);
matrixCompare(seq_string, truth_EA - euler_angles, tol)


%% Cleanup
rmpath('C:\Users\John\Documents\ASEN5010\AIAA_Software\Matlab Toolbox\DOS-Windows\')
rmpath('C:\Users\John\Documents\ASEN5010\ver_tools\')