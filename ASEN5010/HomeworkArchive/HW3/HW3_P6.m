%% Problem 6
fprintf('Problem 6\n');
clearvars -except function_list pub_opt
close all
body_meas(:,1) = [0.8273;0.5541;-0.0920];
body_meas(:,2) = [-0.8285;0.5522;-0.0955];
body_meas(:,3) = [0.2155;0.5522;0.8022];
body_meas(:,4) = [0.5570;-0.7442;-0.2884];

inrtl_meas(:,1) = [-0.1517;-0.9669;0.2050];
inrtl_meas(:,2) = [-0.8393;0.4494;-0.3044];
inrtl_meas(:,3) = [-0.0886;-0.5856;-0.8000];
inrtl_meas(:,4) = [0.8814;-0.0303;0.5202];

% enforce unit vectors
for i = 1:4
    body_meas(:,i) = body_meas(:,i)/norm(body_meas(:,i));
    inrtl_meas(:,i) = inrtl_meas(:,i)/norm(inrtl_meas(:,i));
end

for i = 2:4
    fprintf('DCM for measurements 1 primary and %d secondary:\n', i);
    DCM = triadBI(body_meas(:,1), body_meas(:,i),...
        inrtl_meas(:,1), inrtl_meas(:,i))
end
