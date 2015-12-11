function R = R_ECEF2ENU(lat, lon)
fcnPrintQueue(mfilename('fullpath')); % Add this code to code app

s_lat = sin(lat);
s_lon = sin(lon);
c_lat = cos(lat);
c_lon = cos(lon);
R = [-s_lon c_lon 0;
    -s_lat*c_lon -s_lat*s_lon c_lat;
     c_lat*c_lon  c_lat*s_lon s_lat];