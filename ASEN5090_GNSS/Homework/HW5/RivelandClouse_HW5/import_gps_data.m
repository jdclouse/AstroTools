%==========================================================================
%==========================================================================
% import_gps_data.m
%
% Author: Johnathan Clouse
%
% Function to import GPS track data
%
%==========================================================================
%==========================================================================

function [struct1, struct2] = import_gps_data(file1, file2)
% import_gps_data Import two data that have similar time series to compare
% against eachother
% lat/lon output in rad

% csv data is expected to be lat, lon, alt, ISO 8601 time string
[struct1.lat, struct1.lon, struct1.alt, struct1.time] = ...
    textread(file1, '%f%f%f%s', 'delimiter', ',');
struct1.time = convert_times(struct1.time);

[struct2.lat, struct2.lon, struct2.alt, struct2.time] = ...
    textread(file2, '%f%f%f%s', 'delimiter', ',');
struct2.time = convert_times(struct2.time); 

begin_time = max([struct1.time(1), struct2.time(1)]);
end_time = min([struct1.time(end), struct2.time(end)]);

% Trim both times series to common beginning/end times
trim = @(x,y,t1,t2) x(y >= t1 & y <= t2); %anon fcn for consistency
struct1.lat = trim(struct1.lat,struct1.time, begin_time, end_time);
struct1.lon = trim(struct1.lon,struct1.time, begin_time, end_time);
struct1.alt = trim(struct1.alt,struct1.time, begin_time, end_time);
struct1.time = trim(struct1.time,struct1.time, begin_time, end_time);
struct2.lat = trim(struct2.lat,struct2.time, begin_time, end_time);
struct2.lon = trim(struct2.lon,struct2.time, begin_time, end_time);
struct2.alt = trim(struct2.alt,struct2.time, begin_time, end_time);
struct2.time = trim(struct2.time,struct2.time, begin_time, end_time);

% Convert lat/lon to rads
struct1.lat = struct1.lat*pi/180;
struct1.lon = struct1.lon*pi/180;
struct2.lat = struct2.lat*pi/180;
struct2.lon = struct2.lon*pi/180;

% Deal with lack of data for some times
% Find the common data and stick them in an array
cnt = 1;
series1=zeros(1,4);
series2=zeros(1,4);
for ii = 1:length(struct2.time)
    cond = struct1.time == struct2.time(ii);
    if struct1.time(cond) == struct2.time(ii)
        series1(cnt,:) = [struct1.lat(cond) struct1.lon(cond) struct1.alt(cond) struct1.time(cond)];
        series2(cnt,:) = [struct2.lat(ii) struct2.lon(ii) struct2.alt(ii) struct2.time(ii)];
        cnt = cnt + 1;
    end
end
for ii = 1:length(struct1.time)
    cond = struct2.time == struct1.time(ii);
    if struct2.time(cond) == struct1.time(ii)
        series1(cnt,:) = [struct1.lat(ii) struct1.lon(ii) struct1.alt(ii) struct1.time(ii)];
        series2(cnt,:) = [struct2.lat(cond) struct2.lon(cond) struct2.alt(cond) struct2.time(cond)];
        cnt = cnt + 1;
    end
end

% Remove duplicates, sort by time
series1 = unique(series1, 'rows', 'stable');
series2 = unique(series2, 'rows', 'stable');
series1 = sortrows(series1, 4);
series2 = sortrows(series2, 4);

% Reassign data structs.
struct1.lat = series1(:,1);
struct1.lon = series1(:,2);
struct1.alt = series1(:,3);
struct1.time = series1(:,4);
struct2.lat = series2(:,1);
struct2.lon = series2(:,2);
struct2.alt = series2(:,3);
struct2.time = series2(:,4);
end

function time_nums = convert_times(input)
len = length(input);
time_nums = zeros(len,1);
for ii = 1:len
    time_nums(ii) = datenum8601(input{ii},'*ymdHMS');
end
end