%% Problem 2: S&J, Problem 3.23
fprintf('Problem 2: S&J, Problem 3.23\n');
clearvars -except function_list pub_opt
close all
fprintf('CRP vector:\n');
CRP = [0.5; -0.2; 0.8];
printVector(CRP, '');
Q = vecSkew(CRP);
I = eye(3);
C_Cayley = (I-Q)*inv((I+Q))
C_Eqn_3p120 = ((1-CRP'*CRP)*I + 2*(CRP*CRP') - 2*Q)/(1+CRP'*CRP)
Diff = C_Cayley - C_Eqn_3p120