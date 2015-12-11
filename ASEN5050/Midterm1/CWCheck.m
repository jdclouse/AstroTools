w = sqrt(398600.4415/(6378.1363+380)^3);
X_A1 = [0;0;0;0;0;0.4];
X_B1 = [0;0;0;0;0.2;0];
X_B2 = CWHillSTM(w,8*60)*X_B1;
X_A2 = CWHillSTM(w,8*60)*X_A1;

X_B3 = CWHillSTM(w,20*60)*X_B1;
STM_12_min = CWHillSTM(w,12*60);
V_Aimp_plus = inv(STM_12_min(1:3,4:6))*(X_B3(1:3) - STM_12_min(1:3,1:3)*X_A2(1:3));
V_burn1 = V_Aimp_plus - X_A2(4:6);
V_burn1_wrtB = V_burn1 - X_B2(4:6);
Vf = STM_12_min(4:6,:)*[X_A2(1:3); V_Aimp_plus]