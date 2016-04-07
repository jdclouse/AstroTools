incoming_v = eval(inc_vel(Venus_Earth));
outgoing_v = eval(out_vel(Earth_Jupiter));
for ii = 1:num_VGA_window
    % 
    for jj = 1:num_JGA_window;
        ResoOrb_vel_err = abs(outgoing_v(:,jj) - incoming_v(ii,:)');
        valid_RO = ...
            ResoOrb_vel_err <= max_GA_diff;
        ResoOrb_vel_err_3d(ii,:,jj) = ResoOrb_vel_err;
        ResoOrb_valid(ii,:,jj) = valid_RO;
    end
end