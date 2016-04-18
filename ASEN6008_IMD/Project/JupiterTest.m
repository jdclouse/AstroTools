max_err = 0.3;
best_traj = zeros(6,1);
best_err = max_err;
lowest_vf_traj = zeros(6,1);
lowest_vf = V_inf_final_max;
lowest_vf_err = max_err;
lowest_C3_traj = zeros(6,1);
lowest_C3 = C3_max;
lowest_C3_err = max_err;

get_val_c3 = @(leg) ['valid_c3_' leg.route(1) 'w'];
get_c3_store = @(leg) ['lambert_out(' num2str(leg.idx) ').' ...
    leg.route(1) 'w_c3_store'];


incoming_v = eval(inc_vel(Earth_Venus));
outgoing_v = eval(out_vel(Venus_Earth));

x = 0;
get_val_c3(Earth_Venus)
valid_c3 = eval(get_val_c3(Earth_Venus));
% c3_store = lambert_out(1).lw_c3_store;
c3_store =  eval(get_c3_store(Earth_Venus));
dv_final_store = lambert_out(3).long_way_dv2_store;
dv_final_store = eval(inc_vel(Earth_Jupiter));
all_traj = zeros(4,467485);
for launch_idx = launch_idx1:launch_idx2
%     launch_idx
    for VGA_idx = 1:num_VGA_window
%         VGA_idx
%         if VGA_idx == 172
%             return
%         end
        if valid_c3(launch_idx, VGA_idx) == 0
            continue
        end
        for EGA_idx = 1:num_EGA
%             EGA_idx
            if VGA_valid(launch_idx, VGA_idx, EGA_idx) == 0
                continue
            end
            for JGA_idx = 1:num_JGA_window
                if ResoOrb_valid(VGA_idx, EGA_idx, JGA_idx) == 0
                    continue
                end
%                 use_traj = [launch_idx; VGA_idx; EGA_idx; JGA_idx];
%                 ResonantOrbit;
%                 if ~good_reso 
%                     continue
%                 end
                % Congrats. you made it!
                traj_error = ...
                    +VGA_vel_err_3d(launch_idx, VGA_idx, EGA_idx)...
                    +ResoOrb_vel_err_3d(VGA_idx, EGA_idx, JGA_idx);

                if traj_error < max_err &&...
                        c3_store(launch_idx,VGA_idx) < C3_max &&...
                        dv_final_store(EGA_idx,JGA_idx) < V_inf_final_max
                x = x+1;
                all_traj(:,x) = [launch_idx; VGA_idx; EGA_idx; JGA_idx];
                    if traj_error < best_err
                        best_err = traj_error;
                        best_traj = [launch_idx;
                            VGA_idx; EGA_idx; JGA_idx];
                    end
                    if c3_store(launch_idx,VGA_idx) < lowest_C3
                        lowest_C3 = c3_store(launch_idx,VGA_idx);
                        lowest_C3_err = traj_error;
                        lowest_C3_traj = [launch_idx;
                            VGA_idx; EGA_idx; JGA_idx];
                    end
                    if dv_final_store(EGA_idx,JGA_idx) < lowest_vf
                        lowest_vf = dv_final_store(EGA_idx,JGA_idx);
                        lowest_vf_err = traj_error;
                        lowest_vf_traj = [launch_idx;
                            VGA_idx; EGA_idx; JGA_idx];
                    end
%                             return

                end
            end
        end
    end
end