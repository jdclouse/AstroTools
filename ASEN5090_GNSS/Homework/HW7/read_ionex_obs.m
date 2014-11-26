function [ ionex ] = read_ionex_obs(fname)
% Uses some techniques found in read_rinex_obs, obtained from P. Axelrad
fid=fopen(fname);

if fid == -1
    fprintf('File does not exist!')
    return
end

current_line = fgetl(fid);
header_done = 0;
tec_maps = 1;
exponent = 0;
    
% Process all lines
while (current_line ~= -1)
    num_lats = 0;
    num_lons = 0;
    num_maps = 0;
    % Find the end of the header
    % TODO: Header Parser would be nice
    while ~header_done
        if isempty(strfind(current_line, 'END OF HEADER'))
            if ~isempty(strfind(current_line, '# OF MAPS IN FILE'))
                num_maps = str2num(current_line(1:6));
            end
            if ~isempty(strfind(current_line, 'LAT1 / LAT2 / DLAT'))
                temp = str2num(current_line(1:20));
                ionex.lats = (temp(1):temp(3):temp(2))'*pi/180;
                num_lats = length(ionex.lats);
            end
            if ~isempty(strfind(current_line, 'LON1 / LON2 / DLON'))
                temp = str2num(current_line(1:20));
                ionex.lons = (temp(1):temp(3):temp(2))'*pi/180;
                num_lons = length(ionex.lons);
            end
            if ~isempty(strfind(current_line, 'EXPONENT'))
                temp = str2num(current_line(1:6));
                if temp == -1
                    exponent = 0.1;
                end
            end
            %Find lats
            %
            current_line = fgetl(fid);
        else
            header_done = 1;
            ionex.maps = zeros(num_maps, num_lats, num_lons);
            ionex.epochs = zeros(num_maps,1);
            current_line = fgetl(fid);
        end
    end
    
    
    % Only looking for TEC MAP entries for now
    if isempty(strfind(current_line, 'START OF TEC MAP'))
        current_line = fgetl(fid);
        continue
    end
    current_line = fgetl(fid);
    hr = str2num(current_line(19:24));
    min = str2num(current_line(25:30));
    sec = str2num(current_line(31:36));
    ionex.epochs(tec_maps)=3600*hr + 60 * min + sec;
    
    current_line = fgetl(fid);
    
    collected_tecs = [];
    lat_idx = 0;
    in_map = 1;
    while in_map 
        if ~isempty(strfind(current_line, 'END OF TEC MAP'))
            in_map = 0;
        end
        
        if isempty(strfind(current_line, 'LAT/LON1/LON2/DLON/H')) && in_map
            collected_tecs = [collected_tecs current_line];
            current_line = fgetl(fid);
            continue
        end
        
        if lat_idx > 0
            ionex.maps(tec_maps, lat_idx, :) = str2num(collected_tecs)*exponent;
            collected_tecs = [];
        end
        lat_idx = lat_idx + 1;
        current_line = fgetl(fid);
        
        
    end
    
    tec_maps = tec_maps + 1;
    
    if ~isempty(strfind(current_line, 'START OF RMS MAP'))
        break
    end
    
    1;
end

end