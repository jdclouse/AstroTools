%% SNC process noise
function Q = compute_SNC_Q(fo, X)
%compute_SNC_Q SNC process noise
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

if fo.SNC_use_RIC
    % Find the rotation matrix RIC->Inertial
    r_inrtl = X(1:3)/norm(X(1:3));
    h = cross(X(1:3),X(4:6));
    cross_track_inrtl = h/norm(h);
    in_track = cross(cross_track_inrtl, r_inrtl);

    R_eci_ric = [r_inrtl in_track cross_track_inrtl];
    Q = R_eci_ric'*fo.SNC_Q*R_eci_ric;
else
    Q = fo.SNC_Q;
end

end