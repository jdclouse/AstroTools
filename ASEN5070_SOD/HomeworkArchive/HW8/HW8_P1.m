%% HW8 Problem 1
% John Clouse
%% Initialize
close all
logx_fig = figure;
loglog_fig = figure;
% Epsilon values
len = 10;
eps = ones(len,1)*1e-5; % start at 1e-6 after the loop iterates
for ii = 1:len
    eps(ii:end) = eps(ii:end)/10;
end
std_dev = 1./eps;
R = 1;


Beta = @(e) 1-2*e+2*e*e*(2+e*e);
P2_exact_trace = zeros(len,1);
for ii = 1:len
    e = eps(ii);
    P2_exact_trace(ii) = trace([1+2*e*e, -(1+e);-(1+e), 2+e*e]/Beta(e));
end
clear e

%% CKF P2
P2_ckf_trace = zeros(len,1);
for ii = 1:len
    H1 = [1 eps(ii)];
    P1_ap = eye(2)*std_dev(ii)*std_dev(ii);
    K1 = P1_ap*H1'*inv(H1*P1_ap*H1' + R);
    P1 = (eye(2) - K1*H1)*P1_ap;
    
    %P1 is now P2_ap
    P2_ap = P1;
    H2 = [1 1];
    K2 = P2_ap*H2'*inv(H2*P2_ap*H2' + R);
    P2 = (eye(2) - K2*H2)*P2_ap;
    P2_ckf_trace(ii) = trace(P2);
end

diff_ckf = P2_exact_trace - P2_ckf_trace;
HW8_plot(eps, diff_ckf, 1, 'Kalman', logx_fig, loglog_fig);

%% Joseph P2
I = eye(2);
P2_joseph_trace = zeros(len,1);
for ii = 1:len
    H = [1 eps(ii); 1 1];
    P_ap = I*std_dev(ii)*std_dev(ii);
    for jj = 1:2
        K = P_ap*H(jj,:)'*inv(H(jj,:)*P_ap*H(jj,:)' + R);
        P = (I-K*H(jj,:))*P_ap*(I-K*H(jj,:))' + K*R*K';
        P_ap = P; % a priori for next measurement
    end
    P2_joseph_trace(ii) = trace(P);    
end

diff_joseph = P2_exact_trace - P2_joseph_trace;
HW8_plot(eps, diff_joseph, 2, 'Joseph', logx_fig, loglog_fig);

%% Potter P2
P2_potter_trace = zeros(len,1);
for ii = 1:len    
    H = [1 eps(ii); 1 1];
    P_ap = I*std_dev(ii)*std_dev(ii);
    W_bar = sqrt(P_ap);
    for jj = 1:2
        F = W_bar*H(jj,:)';
        alpha = inv(F'*F + R);
        gamma = 1/(1+sqrt(R*alpha));
        K = alpha*W_bar*F;
        W = W_bar-gamma*K*F';
        W_bar = W; % sequential update
    end
    P2_potter = W*W';
    P2_potter_trace(ii) = trace(P2_potter);
end

diff_potter = P2_exact_trace - P2_potter_trace;
HW8_plot(eps, diff_potter, 3, 'Potter', logx_fig, loglog_fig);

%% Batch P2
P2_batch_trace = zeros(len,1);
for ii = 1:len    
    H = [1 eps(ii); 1 1];
    info_mat = inv(I*std_dev(ii)*std_dev(ii));
    info_mat = info_mat + H(1,:)'*inv(R)*H(1,:);
    info_mat = info_mat + H(2,:)'*inv(R)*H(2,:);
    P2_batch = inv(info_mat);
    P2_batch_trace(ii) = trace(P2_batch);
end

diff_batch = P2_exact_trace - P2_batch_trace;
HW8_plot(eps, diff_batch, 4, 'Batch', logx_fig, loglog_fig);

