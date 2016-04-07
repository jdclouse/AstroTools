max_err = 1.5;
best_traj = zeros(6,1);
best_err = max_err;
lowest_vf_traj = zeros(6,1);
lowest_vf = V_inf_final_max;
lowest_vf_err = max_err;
lowest_C3_traj = zeros(6,1);
lowest_C3 = C3_max;
lowest_C3_err = max_err;
x = 0;
valid_c3 = valid_c3_lw;
c3_store = lambert_out(1).lw_c3_store;
dv_final_store = lambert_out(5).short_way_dv2_store;
for launch_idx = launch_idx1:launch_idx2
    launch_idx
    for VGA_idx = 1:num_VGA_window
        VGA_idx
%         if VGA_idx == 172
%             return
%         end
        if valid_c3(launch_idx, VGA_idx) == 0
            continue
        end
        for EGA_idx = 1:num_EGA
            if VGA_valid(launch_idx, VGA_idx, EGA_idx) == 0
                continue
            end
            for JGA_idx = 1:num_JGA_window
                if ResoOrb_valid(VGA_idx, EGA_idx, JGA_idx) == 0 || ...
                        sum(JGA_valid(EGA_idx, JGA_idx, :)) == 0
                    continue
                end                            
                for SGA_idx = 1:num_SGA
                    if JGA_valid(EGA_idx, JGA_idx, SGA_idx) == 0 || ...
                            sum(SGA_valid(JGA_idx, SGA_idx, :)) == 0
                        continue
                    end                   
                    for NOI_idx = 1:num_NOI
                        if SGA_valid(JGA_idx, SGA_idx, NOI_idx) == 0
                            continue
                        end 
%                         NOI_idx
                    x = x+1;
%                     launch_idx
                        % Congrats. you made it!
                        traj_error = ...
                            +VGA_vel_err_3d(launch_idx, VGA_idx, EGA_idx)...
                            +ResoOrb_vel_err_3d(VGA_idx, EGA_idx, JGA_idx)...
                            +JGA_vel_err_3d(EGA_idx, JGA_idx, SGA_idx)...
                            +abs(lambert_out(5).short_way_dv1_store(SGA_idx,NOI_idx) ...
                                - lambert_out(4).short_way_dv2_store(JGA_idx,SGA_idx)');
                        
                        if traj_error < max_err &&...
                                c3_store(launch_idx,VGA_idx) < C3_max &&...
                                lambert_out(5).short_way_dv1_store(SGA_idx,NOI_idx) < V_inf_final_max
                            if traj_error < best_err
                                best_err = traj_error;
                                best_traj = [launch_idx;
                                    VGA_idx; EGA_idx; JGA_idx; SGA_idx; NOI_idx];
                            end
                            if c3_store(launch_idx,VGA_idx) < lowest_C3
                                lowest_C3 = c3_store(launch_idx,VGA_idx);
                                lowest_C3_err = traj_error;
                                lowest_C3_traj = [launch_idx;
                                    VGA_idx; EGA_idx; JGA_idx; SGA_idx; NOI_idx];
                            end
                            if lambert_out(5).short_way_dv2_store(SGA_idx,NOI_idx) < lowest_vf
                                lowest_vf = lambert_out(5).short_way_dv2_store(SGA_idx,NOI_idx);
                                lowest_vf_err = traj_error;
                                lowest_vf_traj = [launch_idx;
                                    VGA_idx; EGA_idx; JGA_idx; SGA_idx; NOI_idx];
                            end
%                             return
                        end
                    end
                end
            end
        end
    end
end