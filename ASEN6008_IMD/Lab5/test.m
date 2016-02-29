% for each launch date, see what has a valid GA
for ii = 1:1
    % 
    for jj = 1:392%fb_idx1:fb_idx2
    valid_GA = ...
        abs(out2_refined.short_way_dv1_store(:,jj) - out1_refined.short_way_dv2_store(ii,:)') <= .4;
    if sum(valid_GA) > 0
        return
    end
    end
end
sum(valid_GA)