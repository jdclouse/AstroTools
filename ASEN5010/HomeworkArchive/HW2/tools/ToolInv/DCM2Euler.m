function euler_angles = DCM2Euler( seq_string, DCM ) 
%DCM2Euler Turn a DCM set into an Euler angle set in the specified sequence

fcnPrintQueue(mfilename('fullpath'))
euler_angles = zeros(3,1);


if strcmp(seq_string,'321')
    euler_angles(2) = asin(-DCM(1,3));
    euler_angles(1) = atan2(DCM(1,2),DCM(1,1));
    euler_angles(3) = atan2(DCM(2,3),DCM(3,3));
else
    fprintf('this rotation sequence is not supported');
end
end
