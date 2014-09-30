%% HW2 Problem 1: PRN generation
fprintf('\n');
clearvars -except function_list pub_opt
close all

%% a) Plot PRN 19 chips
prn19=PRNCode(19);
for i = 1:1023
prn19.update()
end

binaryVectorToHex((prn19.CA_code(1:16)))
hexToBinaryVector('E6D6');
hw2_code_plot(prn19)

%% b) PRN 19 chips 1024:2046
% These chips are the same as 1:1023
for i = 1:1023
prn19.update()
end

sum(prn19.CA_code(1:1023)-prn19.CA_code(1024:end))

%% c) Plot PRN 25 chips
prn25=PRNCode(25);
for i = 1:1023
prn25.update()
end

hw2_code_plot(prn25)

%% d) Plot PRN 5 chips
prn5=PRNCode(5);
for i = 1:1023
prn5.update()
end

hw2_code_plot(prn5)
