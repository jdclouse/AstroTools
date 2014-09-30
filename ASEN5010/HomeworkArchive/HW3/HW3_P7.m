%% Problem 7
fprintf('Problem 7\n');
clearvars -except function_list pub_opt
close all
body_meas(:,1) = [0.8273;0.5541;-0.0920];
body_meas(:,2) = [-0.8285;0.5522;-0.0955];

inrtl_meas(:,1) = [-0.1517;-0.9669;0.2050];
inrtl_meas(:,2) = [-0.8393;0.4494;-0.3044];

% enforce unit vectors
for i = 1:2
    body_meas(:,i) = body_meas(:,i)/norm(body_meas(:,i));
    inrtl_meas(:,i) = inrtl_meas(:,i)/norm(inrtl_meas(:,i));
end

DCM_triad = triadBI(body_meas(:,1), body_meas(:,2),...
    inrtl_meas(:,1), inrtl_meas(:,2));

w = [1;1];
DCM_davenport = davenportDCM(body_meas, inrtl_meas, w);

% Find the phi about the PRV for each, show error.
% From eqn 3.73 in S&J:
phi_triad = acos(0.5*(trace(DCM_triad) - 1));
phi_davenport = acos(0.5*(trace(DCM_davenport) - 1));
phi_error = abs(phi_triad - phi_davenport);
fprintf('Principle rotation angle between Triad and Davenport solutions:\n');
fprintf('%f degrees\n', phi_error*180/pi);