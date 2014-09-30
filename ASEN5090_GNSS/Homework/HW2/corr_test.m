%% HW2 Problem 2: Code correlations
%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all

prn19=PRNCode(19);
prn25=PRNCode(25);
prn5=PRNCode(5);
make_delay = @(code, delay) [code(delay+1:end), code(1:delay)];
prn19_delay=350;
prn25_delay=905;
prn5_delay=75;

for i = 1:1023
    prn19.update()
    prn25.update()
    prn5.update()
end


offsets_correlation = zeros(1,1024);
%% a) auto-correlation of PRN19
for ii = 0:1023
    offsets_correlation(ii+1) = ...
        normalized_correlation(prn19.CA_code, prn19.CA_code, ii);
end
figure
plot(offsets_correlation)

%% b) cross-correlation of PRN19 and PRN19 delayed by 200
% The peak value's location makes sense because the 2nd code was shifted by
% 200, making the correlation happen 200 bits from the end.
prn19_200delay = make_delay(prn19.CA_code, 200);
for ii = 0:1023
    offsets_correlation(ii+1) = ...
        normalized_correlation(prn19.CA_code, prn19_200delay, ii);
end
figure
plot(offsets_correlation)

%% c) cross-correlation of PRN19 and PRN25
% No good correlations between PRNs 19 and 25
for ii = 0:1023
    offsets_correlation(ii+1) = ...
        normalized_correlation(prn19.CA_code, prn25.CA_code, ii);
end
figure
plot(offsets_correlation)

%% d) cross-correlation of PRN19 and PRN5
% No good correlations between PRNs 19 and 5
for ii = 0:1023
    offsets_correlation(ii+1) = ...
        normalized_correlation(prn19.CA_code, prn5.CA_code, ii);
end
figure
plot(offsets_correlation)

%% e) cross-correlation of PRN19 and summed PRNs
% The peak value's location makes sense because the 2nd code was shifted by
% 350, making the correlation happen 350 bits from the end.  The codes
% correlate despite the other codes on top of PRN 19.
x1 = make_delay(prn19.CA_code, prn19_delay);
x2 = make_delay(prn25.CA_code, prn25_delay);
x3 = make_delay(prn5.CA_code, prn5_delay);
summed_prns = x1+x2+x3;
for ii = 0:1023
    offsets_correlation(ii+1) = ...
        normalized_correlation(prn19.CA_code, summed_prns, ii);
end
figure
plot(offsets_correlation)

%% f) Noise
noise = 4*randn(1,1023);
figure
subplot(2, 2, 1)
plot(x1)
title('x_1')
xlabel('chip')
ylabel('bit value')
ylim([min(noise)*1.1, max(noise)*1.1])
xlim([1, 1023])
subplot(2, 2, 2)
plot(x2)
title('x_2')
xlabel('chip')
ylabel('bit value')
ylim([min(noise)*1.1, max(noise)*1.1])
xlim([1, 1023])
subplot(2, 2, 3)
plot(x3)
title('x_3')
xlabel('chip')
ylabel('bit value')
ylim([min(noise)*1.1, max(noise)*1.1])
xlim([1, 1023])
subplot(2, 2, 4)
plot(noise)
title('Noise')
xlabel('chip')
ylabel('Noise value')
ylim([min(noise)*1.1, max(noise)*1.1])
xlim([1, 1023])

%% g) cross-correlation of PRN19 and summed PRNs + noise
% The peak is where I expect it to be, however it's much closer in scale to
% other local peaks. I ran the noise function multiple times and sometimes
% the correlation really stood out, sometimes it didn't. 
summed_prns = x1+x2+x3+noise;
for ii = 0:1023
    offsets_correlation(ii+1) = ...
        normalized_correlation(prn19.CA_code, summed_prns, ii);
end
figure
plot(offsets_correlation)