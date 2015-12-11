function JD = computeJD(yr,mo,day,hr,mn,sec, is_leap_sec_day)
%computeJD return Julian date for UT1
% ONLY VALID BETWEEN 1900-2100!!!!
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

sec_denom = 60;
if nargin < 6
    fprintf('ERROR - not enough arguments for computeJD.\n')
    JD = -1;
    return;
elseif nargin == 7
    if is_leap_sec_day == 1
        sec_denom = 61;
    end
end
JD = 367*yr -floor(7/4*(yr+floor((mo+9)/12)))...
    + floor(275*mo/9) + day + 1721013.5 ...
    + ((sec/sec_denom + mn)/60 + hr)/24;