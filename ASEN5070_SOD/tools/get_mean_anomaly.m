function M = get_mean_anomaly(E, e)
%get_mean_anomaly Calculate mean anomaly
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

M = E - e.*sin(E);