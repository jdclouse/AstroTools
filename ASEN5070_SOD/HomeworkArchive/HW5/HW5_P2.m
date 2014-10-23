%% HW5 Problem 2
%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all
% Bring in answers to compare
hw5_p1_answers

%% Find x_est and observation error with least squares
y = [1 2 1]';
W = [2 0 0; 0 1 0; 0 0 1];
H = [1 1 1]';
x_bar = 2;
W_bar = 2;

% BLS algorithm, just one pass
lam = W_bar;
N=lam*x_bar;

lam = lam + H'*W*H;
N = N + H'*W*y;

x_est = lam\N

error = y - H*x_est